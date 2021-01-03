#Script finds consensus sequence from R1 and R2 reads by using the Gotoh-algorithm
#process_files() function should be run with the following arguments:
#list_PF_path_input_files: path to fastq-files
#str_PF_path_output_file: path to output-file
#int_PF_minLength: minimum length of consensus sequence
#float_PF_minScoreFraction: minimum alignment score between sequences as fraction of length of consensus sequence
#bool_PF_cut_reads: if genomic region between primers is long, reads can be cut to overlapping region
#int_PF_gReg_length: length of genomic region, only needed if bool_PF_cut_reads=True
#list_PF_max_read_lengths: max position of reads to be used, only needed if bool_PF_cut_reads=True
#int_PF_matchScore: score for matching bases
#int_PF_arraySize: max length of illumina reads
#int_PF_read_chunk: number of reads to be processed per thread
#int_PF_number_processes: number of threads to be run in parallel, 0 defaults to number of CPU-cores

import sys
import os
from time import monotonic
import contextlib
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import concurrent.futures
from LocalGotohMergeClass import LocalGotohMerge
from NGSMRDFuncs import RevSeq,Chr2Phred,Phred2Chr,fastq_zip_equal,fastq_iterator,iter_double_chunked,format_time,format_fastq,fastq_get_read_number
import argparse



def merge_read_pair(instance_MRP_LocalGotoh,list_MRP_reads,int_MRP_matchScore,int_MRP_minLength,float_RMP_minScoreFraction):
	#check if read sequences are equal
	if list_MRP_reads[0]['read'] == list_MRP_reads[1]['read']:
		int_MPR_len_reads = len(list_MRP_reads[0]['read'])
		str_MRP_output_quals = "".join([Phred2Chr(max(Chr2Phred(list_MRP_reads[0]['quals'][int_MRP_i]), \
			Chr2Phred(list_MRP_reads[1]['quals'][int_MRP_i]))) for int_MRP_i in range(int_MPR_len_reads)])
		list_MRP_merge_output_read = [int_MPR_len_reads*int_MRP_matchScore, [list_MRP_reads[0]['read'],str_MRP_output_quals]]

	else:
		list_MRP_merge_output_read = instance_MRP_LocalGotoh.alnMerge([list_MRP_reads[0]['read'],list_MRP_reads[0]['quals']], \
			[list_MRP_reads[1]['read'],list_MRP_reads[1]['quals']])

	str_MRP_output_read_name = list_MRP_reads[0]['name'].split(' ',1)[0] + ' R1R2'
	int_MRP_output_read_len = len(list_MRP_merge_output_read[1][0])

	if int_MRP_output_read_len >= int_MRP_minLength and list_MRP_merge_output_read[0] >= int_MRP_output_read_len*float_RMP_minScoreFraction:
		bool_MRP_match_filter = True
	else:
		bool_MRP_match_filter = False

	return [bool_MRP_match_filter,{'name':str_MRP_output_read_name,'read':list_MRP_merge_output_read[1][0],'quals':list_MRP_merge_output_read[1][1]}]


def cut_reads(list_CR_reads,int_CR_gReg_length,list_CR_max_read_lengths):
	list_CR_read_lengths = [len(str_CR_temp_read['read']) if len(str_CR_temp_read['read']) <= int_CR_temp_max_length else int_CR_temp_max_length \
		for int_CR_temp_max_length,str_CR_temp_read in zip(list_CR_max_read_lengths,list_CR_reads)]

	list_CR_cut_points = [int_CR_gReg_length - int_CR_temp_read_length for int_CR_temp_read_length in list_CR_read_lengths[::-1]]

	list_CR_output_reads = []

	for int_CR_num, int_CR_temp_cut_point in enumerate(list_CR_cut_points):
		if 0 <= int_CR_temp_cut_point <= list_CR_read_lengths[int_CR_num]:
			str_CR_temp_output_read = list_CR_reads[int_CR_num]['read'][int_CR_temp_cut_point:list_CR_read_lengths[int_CR_num]]
			str_CR_temp_output_quals = list_CR_reads[int_CR_num]['quals'][int_CR_temp_cut_point:list_CR_read_lengths[int_CR_num]]
		else:
			str_CR_temp_output_read = list_CR_reads[int_CR_num]['read']
			str_CR_temp_output_quals = list_CR_reads[int_CR_num]['quals']

		list_CR_output_reads.append({'name':list_CR_reads[int_CR_num]['name'],'read':str_CR_temp_output_read,'quals':str_CR_temp_output_quals})

	return list_CR_output_reads


def process_read_chunk(list_PRC_reads,bool_PRC_cut_reads,int_PRC_gReg_length,list_PRC_max_read_lengths,int_PRC_matchScore,int_PRC_minLength,float_PRC_minScoreFraction,int_PRC_arraySize):
	instance_PRC_LocalGotoh = LocalGotohMerge(int_PRC_arraySize,int_PRC_arraySize)

	list_PRC_output = []

	for list_PRC_temp_read in list_PRC_reads:
		if bool_PRC_cut_reads is True:
			list_PRC_temp_read = cut_reads(list_PRC_temp_read,int_PRC_gReg_length,list_PRC_max_read_lengths)

		list_PRC_temp_read[1]['read'] = RevSeq(list_PRC_temp_read[1]['read'])
		list_PRC_temp_read[1]['quals'] = list_PRC_temp_read[1]['quals'][::-1]

		list_PRC_output.append(merge_read_pair(instance_PRC_LocalGotoh,list_PRC_temp_read,int_PRC_matchScore,int_PRC_minLength,float_PRC_minScoreFraction))

	return list_PRC_output


def process_files(list_PF_path_input_files,str_PF_path_output_file,int_PF_minLength,float_PF_minScoreFraction,bool_PF_cut_reads,int_PF_gReg_length,list_PF_max_read_lengths, \
	int_PF_matchScore=1,int_PF_arraySize=241,int_PF_read_chunk=10000,int_PF_number_processes=0):
	#number processes defaults to cpu count
	int_PF_number_processes = int_PF_number_processes if int_PF_number_processes > 0 else os.cpu_count()

	float_PF_time_start = monotonic()
	print('Open input file: ' + list_PF_path_input_files[0])
	print('Open input file: ' + list_PF_path_input_files[1])
	int_PF_total_reads = fastq_get_read_number(*list_PF_path_input_files)

	#open files
	with contextlib.ExitStack() as stack:
		list_PF_input_files = [stack.enter_context(open(str_PF_temp_path_input_file, 'r')) \
			for str_PF_temp_path_input_file in list_PF_path_input_files]

		it_PF_input_fastq = fastq_iterator(fastq_zip_equal(*[FastqGeneralIterator(obj_PF_temp_input_file) \
			for obj_PF_temp_input_file in list_PF_input_files]))

		obj_PF_output_file = stack.enter_context(open(str_PF_path_output_file, 'w'))

		int_PF_reads_processed = 0
		for list_PF_reads_for_processing in iter_double_chunked(it_PF_input_fastq,int_PF_read_chunk,int_PF_number_processes):
			with concurrent.futures.ProcessPoolExecutor() as executor:
				list_PF_reads_processed = executor.map(process_read_chunk,list_PF_reads_for_processing,*(itertools.repeat(obj_PF_temp_arg) \
					for obj_PF_temp_arg in (bool_PF_cut_reads,int_PF_gReg_length,list_PF_max_read_lengths,int_PF_matchScore,int_PF_minLength,float_PF_minScoreFraction,int_PF_arraySize)))

				for list_PF_temp_read in itertools.chain.from_iterable(list_PF_reads_processed):
					if list_PF_temp_read[0] is True:
						#write read to file
						obj_PF_output_file.write(format_fastq(list_PF_temp_read[1]))
					int_PF_reads_processed += 1
				

			str_PF_percent_processed = f'{round(int_PF_reads_processed*100/int_PF_total_reads,2):.2f}%' \
				if int_PF_total_reads > 0 else '100.00%'

			sys.stdout.write(f'\rProcessing Reads: {int_PF_reads_processed} / {int_PF_total_reads} ({str_PF_percent_processed})')
			sys.stdout.flush()

	float_PF_time_total = monotonic() - float_PF_time_start
	print(f'\nEnd of files: processed {int_PF_reads_processed} reads in {format_time(float_PF_time_total)}')	

def main():

	merge_argparse = argparse.ArgumentParser(description='Merge R1 and R2 NGS-Reads')
	merge_argparse.add_argument('--input_file1', type=str, required=True, help='Fastq file forward reads')
	merge_argparse.add_argument('--input_file2', type=str, required=True, help='Fastq file reverse reads')
	merge_argparse.add_argument('--output_file', type=str, required=True, help='Path and file for output file')
	merge_argparse.add_argument('--min_length', type=int, required=True, help='minimum length of consensus sequence')
	merge_argparse.add_argument('--score_fraction_percent', type=int, required=True, help='minimum score of consensus in relation to the length of the consensus')
	merge_argparse.add_argument('--cut_reads', type=int, required=True, choices=[0,1], help='cut reads to overlapping region')
	merge_argparse.add_argument('--gReg_length', type=int, required=True, help='length of covered genomic region')
	merge_argparse.add_argument('--cap_read_f_length', type=int, required=True, help='cap read length, only active if --cut_reads=1')
	merge_argparse.add_argument('--cap_read_r_length', type=int, required=True, help='cap read length, only active if --cut_reads=1')

	merge_args = merge_argparse.parse_args()

	float_merge_score_fraction = merge_args.score_fraction_percent / 100

	bool_merge_cut_reads = True if merge_args.cut_reads == 1 else False

	process_files([merge_args.input_file1,merge_args.input_file2],merge_args.output_file,merge_args.min_length,float_merge_score_fraction, \
		bool_merge_cut_reads,merge_args.gReg_length,[merge_args.cap_read_f_length,merge_args.cap_read_r_length])

if __name__ == '__main__':
	main()