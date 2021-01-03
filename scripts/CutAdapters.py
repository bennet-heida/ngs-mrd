#Cut NGS-MRD reads to genomic region and look for random tag for read family analysis
#Reads are valid if genomic sequence of forward primer is found, RF-tag is valid if common sequence is found and
#number of bases between common sequence and forward primer is 8 or 16
#process_files() function should be run with following arguments:
#list_PF_path_input_files: list of paths to input fastq-files
#str_PF_prefix_output_files: path and prefix for output_files
#str_PF_CS_f/r: base sequence of common sequences
#str_PF_gPrimer_f/r: genomic primer sequences
#int_PF_length_gReg: length of genomic region between primers
#int_PF_length_indel: length of indel target, 0 for SNP-targets
#int_PF_maxPen_x: max number of penalties allowed in sequence
#int_PF_buffer_x: number of bases around expected position to be included in search for sequence, higher number increases chance of finding sequence,
#but impacts speed of script
#list_PF_RF_lengths: list of valid RF-tag lengths
#int_PF_max_read_length: max read length of illumina sequencing run
#int_PF_max_primer_length: max length of used primers
#int_PF_read_chunks: size of read chunks to be processed per thread
#int_PF_number_processes: number of threads, leave at 0 to default to CPU-cores

import sys
import os
from time import monotonic
import contextlib
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import concurrent.futures
from LocalNeedleAlnClass import LocalNeedleAln
from NGSMRDFuncs import format_time, format_fastq, fastq_zip_equal, fastq_iterator, iter_chunked, iter_double_chunked, RevSeq, fastq_get_read_number
import argparse



def aln_speed(obj_AS_LocalNeedle,str_AS_read,str_AS_query,bool_allowEnd=False):
	int_AS_exact_match = str_AS_read.find(str_AS_query)
	#check for exact match
	if int_AS_exact_match != -1:
		list_AS_LN_output = [int_AS_exact_match,int_AS_exact_match+len(str_AS_query),0]
	#if no exact match was found, perform error tolerant search
	else:
		list_AS_LN_output = obj_AS_LocalNeedle.alnQuery(str_AS_read,str_AS_query)
	return list_AS_LN_output


def mark_read(obj_MR_LocalNeedle,str_MR_read,dict_MR_seqs,dict_MR_slice_points,dict_MR_maxPenalties):
	list_MR_search_order = ['gPrimer','CS','gPrimer_rev']
	list_MR_found = [{'CS':False,'gPrimer':False,'gPrimer_rev':False},{'CS':[0,0],'gPrimer':[0,0],'gPrimer_rev':[0,0]}]
	#check read length, what points can be found on the read
	int_MR_len_read = len(str_MR_read)

	for str_MR_point in list_MR_search_order:
		#only perform search if target is expected within read
		if dict_MR_slice_points[str_MR_point][0] <= int_MR_len_read - (len(dict_MR_seqs[str_MR_point]) - dict_MR_maxPenalties[str_MR_point]):
			str_MR_temp_read_slice = str_MR_read[dict_MR_slice_points[str_MR_point][0]:(dict_MR_slice_points[str_MR_point][1] if dict_MR_slice_points[str_MR_point][1] <= int_MR_len_read else int_MR_len_read)]
			list_MR_temp_LN_output = aln_speed(obj_MR_LocalNeedle,str_MR_temp_read_slice,dict_MR_seqs[str_MR_point])

			if list_MR_temp_LN_output[2] <= dict_MR_maxPenalties[str_MR_point]:
				list_MR_found[0][str_MR_point] = True
				list_MR_found[1][str_MR_point] = [dict_MR_slice_points[str_MR_point][0] + list_MR_temp_LN_output[0],dict_MR_slice_points[str_MR_point][0] + list_MR_temp_LN_output[1]]
			else:
				#if gPrimer pos. is not found, break loop
				if str_MR_point == 'gPrimer':
					break

	return list_MR_found

def calc_slice_points(dict_CSP_seqs,int_CSP_length_gReg,int_CSP_length_indel,dict_CSP_buffer,list_CSP_RF_lengths):
	dict_CSP_seq_lengths = {'CS':len(dict_CSP_seqs['CS']),'gPrimer':len(dict_CSP_seqs['gPrimer']),'gPrimer_rev':len(dict_CSP_seqs['gPrimer_rev'])}
	dict_CSP_slice_points = {}
	dict_CSP_slice_points['CS'] = [0,dict_CSP_seq_lengths['CS']+dict_CSP_buffer['CS']]
	dict_CSP_slice_points['gPrimer'] = [dict_CSP_seq_lengths['CS'] + min(list_CSP_RF_lengths) - dict_CSP_buffer['gPrimer'], \
		dict_CSP_seq_lengths['CS'] + max(list_CSP_RF_lengths) + dict_CSP_seq_lengths['gPrimer'] + dict_CSP_buffer['gPrimer']]
	#check first base (in case large buffer was used)
	if dict_CSP_slice_points['gPrimer'][0] < 0: dict_CSP_slice_points['gPrimer'][0] = 0

	dict_CSP_slice_points['gPrimer_rev'] = [dict_CSP_seq_lengths['CS'] + min(list_CSP_RF_lengths) + dict_CSP_seq_lengths['gPrimer'] + int_CSP_length_gReg \
		- dict_CSP_buffer['gPrimer_rev'] + (int_CSP_length_indel if int_CSP_length_indel < 0 else 0), \
		dict_CSP_seq_lengths['CS'] + max(list_CSP_RF_lengths) + dict_CSP_seq_lengths['gPrimer'] + int_CSP_length_gReg + dict_CSP_seq_lengths['gPrimer_rev'] \
		+ dict_CSP_buffer['gPrimer_rev'] + (int_CSP_length_indel if int_CSP_length_indel > 0 else 0)]
	#check first base (in case large buffer was used)
	if dict_CSP_slice_points['gPrimer_rev'][0] < dict_CSP_seq_lengths['CS'] + min(list_CSP_RF_lengths) + dict_CSP_seq_lengths['gPrimer']:
		dict_CSP_slice_points['gPrimer_rev'][0] = dict_CSP_seq_lengths['CS'] + min(list_CSP_RF_lengths) + dict_CSP_seq_lengths['gPrimer']

	return dict_CSP_slice_points

def process_read_pair(obj_PRP_LocalNeedle,list_PRP_reads,list_PRP_seqs,list_PRP_slice_points,dict_PRP_maxPenalties,list_PRP_RF_lengths):
	list_PRP_cut_seqs = ['read','quals']

	list_PRP_output = [False,[{},{}]]#status for read cut, and read itself
	str_PRP_RF = ""

	list_PRP_marks = [mark_read(obj_PRP_LocalNeedle,dict_PRP_temp_read['read'],dict_PRP_temp_seqs,dict_PRP_temp_slice_points,dict_PRP_maxPenalties) \
		for dict_PRP_temp_read,dict_PRP_temp_seqs,dict_PRP_temp_slice_points in zip(list_PRP_reads,list_PRP_seqs,list_PRP_slice_points)]

	#check if gPrimer was found on both reads
	if all([list_PRP_temp_mark[0]['gPrimer'] is True for list_PRP_temp_mark in list_PRP_marks]):
		list_PRP_output[0] = True
		#if CS was found on both reads, extract RF-sequence
		if all([list_PRP_temp_mark[0]['CS'] is True for list_PRP_temp_mark in list_PRP_marks]):
			#get RF-sequences
			list_PRP_RF_seqs = [dict_PRP_temp_read['read'][dict_PRP_temp_mark[1]['CS'][1]:dict_PRP_temp_mark[1]['gPrimer'][0]] \
				for dict_PRP_temp_read, dict_PRP_temp_mark in zip(list_PRP_reads,list_PRP_marks)]
			#if RF-length matches pattern RF is found
			int_PRP_RF_length = len(list_PRP_RF_seqs[0]) if [len(int_PRP_temp_RF) for int_PRP_temp_RF in list_PRP_RF_seqs].count(len(list_PRP_RF_seqs[0])) \
				== len(list_PRP_RF_seqs) else 0
			#concatenate RF-seqs
			if int_PRP_RF_length in list_PRP_RF_lengths:
				str_PRP_RF = "-".join(list_PRP_RF_seqs)

	for int_PRP_n, dict_PRP_temp_read_input in enumerate(list_PRP_reads):
		list_PRP_output[1][int_PRP_n]['name'] = str_PRP_RF + '/' + dict_PRP_temp_read_input['name']
		list_PRP_output[1][int_PRP_n]['read'] = dict_PRP_temp_read_input['read'][list_PRP_marks[int_PRP_n][1]['gPrimer'][1] if list_PRP_output[0] else 0: \
			list_PRP_marks[int_PRP_n][1]['gPrimer_rev'][0] if list_PRP_output[0] and list_PRP_marks[int_PRP_n][0]['gPrimer_rev'] else len(dict_PRP_temp_read_input['read'])]
		list_PRP_output[1][int_PRP_n]['quals'] = dict_PRP_temp_read_input['quals'][list_PRP_marks[int_PRP_n][1]['gPrimer'][1] if list_PRP_output[0] else 0: \
			list_PRP_marks[int_PRP_n][1]['gPrimer_rev'][0] if list_PRP_output[0] and list_PRP_marks[int_PRP_n][0]['gPrimer_rev'] else len(dict_PRP_temp_read_input['quals'])]

	return list_PRP_output


def process_read_chunk(list_PRC_read_chunk,list_PRC_seqs,list_PRC_slice_points,dict_PRC_maxPenalties,list_PRC_RF_lengths,int_PRC_max_read_length,int_PRC_max_primer_length):
	obj_PRC_LocalNeedle = LocalNeedleAln(int_PRC_max_read_length,int_PRC_max_primer_length,1,-1,-1)

	list_PRC_read_chunk_output = []

	for list_PRC_read_item in list_PRC_read_chunk:
		list_PRC_read_chunk_output.append(process_read_pair(obj_PRC_LocalNeedle,list_PRC_read_item,list_PRC_seqs, \
			list_PRC_slice_points,dict_PRC_maxPenalties,list_PRC_RF_lengths))

	return list_PRC_read_chunk_output


def process_files(list_PF_path_input_files,str_PF_prefix_output_files,str_PF_CS_f,str_PF_CS_r,str_PF_gPrimer_f,str_PF_gPrimer_r,int_PF_length_gReg,int_PF_length_indel, \
	int_PF_maxPen_CS,int_PF_maxPen_gPrimer,int_PF_maxPen_gPrimer_rev,int_PF_buffer_CS,int_PF_buffer_gPrimer,int_PF_buffer_gPrimer_rev,list_PF_RF_lengths, \
	int_PF_max_read_length=241,int_PF_max_primer_length=30,int_PF_read_chunks=10000,int_PF_number_processes=0):
	#number processes defaults to cpu count
	int_PF_number_processes = int_PF_number_processes if int_PF_number_processes > 0 else os.cpu_count()
	#get start time
	float_PF_time_start = monotonic()
	#create dictionaries
	list_PF_seqs = [{'CS':str_PF_CS_f,'gPrimer':str_PF_gPrimer_f,'gPrimer_rev':RevSeq(str_PF_gPrimer_r)}, \
		{'CS':str_PF_CS_r,'gPrimer':str_PF_gPrimer_r,'gPrimer_rev':RevSeq(str_PF_gPrimer_f)}]
	dict_PF_maxPenalties = {'CS':int_PF_maxPen_CS,'gPrimer':int_PF_maxPen_gPrimer,'gPrimer_rev':int_PF_maxPen_gPrimer_rev}
	dict_PF_buffer = {'CS':int_PF_buffer_CS,'gPrimer':int_PF_buffer_gPrimer,'gPrimer_rev':int_PF_buffer_gPrimer_rev}
	#get read count
	print(f'Open input file: {list_PF_path_input_files[0]}')
	print(f'Open input file: {list_PF_path_input_files[1]}')
	int_PF_read_count_input_files = fastq_get_read_number(*list_PF_path_input_files)

	list_PF_slice_points = [calc_slice_points(dict_PF_seqs_temp,int_PF_length_gReg,int_PF_length_indel,dict_PF_buffer,list_PF_RF_lengths) \
		for dict_PF_seqs_temp in list_PF_seqs]

	list_PF_temp_output_endings, list_PF_temp_output_discarded_endings = ['_cut_F.fq','_cut_R.fq'], ['_cut_discarded_F.fq','_cut_discarded_R.fq']

	list_PF_path_output_files = [str_PF_prefix_output_files + str_PF_temp_output_end \
		for str_PF_temp_output_end in list_PF_temp_output_endings]
	list_PF_path_output_discarded_files = [str_PF_prefix_output_files + str_PF_temp_output_end \
		for str_PF_temp_output_end in list_PF_temp_output_discarded_endings]

	with contextlib.ExitStack() as stack:
		#Input files
		list_PF_input_files = [stack.enter_context(open(str_PF_temp_path_input_file, 'r')) for str_PF_temp_path_input_file \
			in list_PF_path_input_files]
		it_PF_input_fastq = fastq_iterator(fastq_zip_equal(*[FastqGeneralIterator(obj_PF_temp_input_file) \
			for obj_PF_temp_input_file in list_PF_input_files]))

		#Output files and output files for discarded reads
		list_PF_output_files = [stack.enter_context(open(str_PF_temp_path_output_file, 'w')) \
			for str_PF_temp_path_output_file in list_PF_path_output_files]
		list_PF_output_discarded_files = [stack.enter_context(open(str_PF_temp_path_output_discarded_file, 'w')) \
			for str_PF_temp_path_output_discarded_file in list_PF_path_output_discarded_files]

		dict_PF_output_switch = {True:list_PF_output_files,False:list_PF_output_discarded_files}

		int_PF_reads_processed = 0

		for list_PF_reads_for_processing in iter_double_chunked(it_PF_input_fastq,int_PF_read_chunks,int_PF_number_processes):
			with concurrent.futures.ProcessPoolExecutor() as executor:
				list_PF_reads_processed = executor.map(process_read_chunk,list_PF_reads_for_processing,*(itertools.repeat(obj_PF_temp_arg) \
					for obj_PF_temp_arg in (list_PF_seqs,list_PF_slice_points,dict_PF_maxPenalties,list_PF_RF_lengths,int_PF_max_read_length,int_PF_max_primer_length)))
				
				for list_PF_temp_read_output in itertools.chain.from_iterable(list_PF_reads_processed):
					for int_PF_temp_read_FR, dict_PF_temp_read in enumerate(list_PF_temp_read_output[1]):
						dict_PF_output_switch[list_PF_temp_read_output[0]][int_PF_temp_read_FR].write(format_fastq(dict_PF_temp_read))

					int_PF_reads_processed += 1

			str_PF_percent_processed = f'{round(int_PF_reads_processed*100/int_PF_read_count_input_files,2):.2f}%' \
				if int_PF_read_count_input_files > 0 else '100.00%'

			sys.stdout.write(f'\rProcessing Reads: {int_PF_reads_processed} / {int_PF_read_count_input_files} ({str_PF_percent_processed})')
			sys.stdout.flush()

	float_PF_time_total = monotonic() - float_PF_time_start
	print(f'\nEnd of files: processed {int_PF_reads_processed} reads in {format_time(float_PF_time_total)}')
			


def main():
	cut_argparse = argparse.ArgumentParser(description='Cut NGS-Reads')
	cut_argparse.add_argument('--input_file1', type=str, required=True, help='Fastq file forward reads')
	cut_argparse.add_argument('--input_file2', type=str, required=True, help='Fastq file reverse reads')
	cut_argparse.add_argument('--output_prefix', type=str, required=True, help='Path and file prefix for output files')
	cut_argparse.add_argument('--seq_CS_F', type=str, required=True, help='forward common sequence')
	cut_argparse.add_argument('--seq_CS_R', type=str, required=True, help='reverse common sequence')
	cut_argparse.add_argument('--seq_gPrimer_F', type=str, required=True, help='forward genomic primer')
	cut_argparse.add_argument('--seq_gPrimer_R', type=str, required=True, help='reverse genomic primer')
	cut_argparse.add_argument('--length_gReg', type=int, required=True, help='length of covered genomic region between primers')
	cut_argparse.add_argument('--length_indel', type=int, required=True, help='length of target insertion or deletion, positive value means insertion, negative value deletion, 0 for SNPs')
	cut_argparse.add_argument('--maxPen_CS', type=int, required=True, help='number of penalties allowed in CS match')
	cut_argparse.add_argument('--maxPen_gPrimer', type=int, required=True, help='number of penalties allowed in genomic primer match')
	cut_argparse.add_argument('--maxPen_gPrimer_rev', type=int, required=True, help='number of penalties allowed in reverse genomic primer match')
	cut_argparse.add_argument('--buffer_CS', type=int, required=True, help='buffer for search position in CS match')
	cut_argparse.add_argument('--buffer_gPrimer', type=int, required=True, help='buffer for search position in genomic primer match')
	cut_argparse.add_argument('--buffer_gPrimer_rev', type=int, required=True, help='buffer for search position in reverse genomic primer match')
	cut_argparse.add_argument('--max_read_length', type=int, required=True, help='max length of reads')
	cut_argparse.add_argument('--max_primer_length', type=int, required=True, help='max length of primers')
	cut_argparse.add_argument('--read_chunks', type=int, required=True, help='reads to be processed at a time')
	cut_argparse.add_argument('--number_processes', type=int, required=True, help='number of processes to be run in parallel')

	cut_args = cut_argparse.parse_args()

	list_cut_input_files = [cut_args.input_file1,cut_args.input_file2]
	list_cut_RF_lengths = [8,16]

	process_files(list_cut_input_files,cut_args.output_prefix,cut_args.seq_CS_F,cut_args.seq_CS_R,cut_args.seq_gPrimer_F,cut_args.seq_gPrimer_R, \
		cut_args.length_gReg,cut_args.length_indel,cut_args.maxPen_CS,cut_args.maxPen_gPrimer,cut_args.maxPen_gPrimer_rev,cut_args.buffer_CS, \
		cut_args.buffer_gPrimer,cut_args.buffer_gPrimer_rev,list_cut_RF_lengths,cut_args.max_read_length,cut_args.max_primer_length, \
		cut_args.read_chunks,cut_args.number_processes,)


if __name__ == '__main__':
	main()