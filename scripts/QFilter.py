#Quality filtering for NGS-Reads
#Script calculates percentage of bases under phred quality threshold. If percentage is below cut-off read is discarded.
#End-trimming looks for ends of reads, that drop off in quality and removes end from reads.
#minimum amount of reads to be kept can be set, this results in lowered cut-off to increase the amount of reads, that pass filter
#process_files() function should be run with the following arguments:
#list_PF_path_input_files: list of input fastq-files
#str_PF_path_output_prefix: path and prefix for output files
#int_PF_phred_cutoff: phred-score threshold, bases below this score will be considered low-quality
#int_PF_input_max_low_q_percent: cut-off for maximum percentage of bases per read to be below phred-score threshold
#int_PF_max_read_filter_percent: max percentage of reads to be filtered out, if not enough reads pass the cut-off, cut-off will be lowered
#to reach at least this percentage of reads to pass the filter
#int_PF_start_end_trim: no end-trimming before this position
#int_PF_end_end_trim: end trimming can't start after this position
#int_PF_end_trim_max_low_q_percent: cut-off for maximum percentage of bases to be below phred-score threshold for end-trimming
#int_PF_read_chunks: number of reads to be processed per thread
#int_PF_number_processes: number of threads to be run in parallel, leave at 0 to default to number of CPU-cores

import sys
import os
import contextlib
import concurrent.futures
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from NGSMRDFuncs import Chr2Phred,format_fastq,format_time,fastq_get_read_number,iter_double_chunked,fastq_zip_equal,fastq_iterator,division_zero_tolerant
from time import monotonic
import itertools
from math import ceil
import argparse

def get_quality(dict_GQ_read,int_GQ_phred_cutoff):
	list_GQ_low_phreds = [1 if Chr2Phred(str_GQ_temp_chr) < int_GQ_phred_cutoff else 0 for str_GQ_temp_chr in dict_GQ_read['quals']]

	int_GQ_len_read = len(dict_GQ_read['read'])

	float_GQ_low_q_fraction = division_zero_tolerant(sum(list_GQ_low_phreds),int_GQ_len_read)

	int_GQ_low_q_fraction_ceil = ceil(float_GQ_low_q_fraction*100)

	return int_GQ_low_q_fraction_ceil

def quality_read_pair(list_QRP_reads,int_QRP_phred_cutoff):
	list_QRP_low_q_fraction_ceil = [get_quality(dict_QRP_temp_read,int_QRP_phred_cutoff) for dict_QRP_temp_read in list_QRP_reads]

	int_QRP_low_q_fraction_ceil_max = max(list_QRP_low_q_fraction_ceil)

	return int_QRP_low_q_fraction_ceil_max

def quality_read_chunk(list_QRC_read_chunk,int_QRC_phred_cutoff):
	list_QRC_quality_distr_output = [0 for _ in range(101)]

	for list_QRC_temp_read in list_QRC_read_chunk:
		int_QRC_temp_low_q_fraction_ceil_max = quality_read_pair(list_QRC_temp_read,int_QRC_phred_cutoff)

		list_QRC_quality_distr_output[int_QRC_temp_low_q_fraction_ceil_max] += 1

	return list_QRC_quality_distr_output


def adjust_quality_distr(list_AQD_path_input_files,int_AQD_phred_cutoff,int_AQD_input_max_low_q_percent,int_AQD_max_read_filter_percent,int_AQD_read_chunks,int_AQD_number_processes):

	list_AQD_quality_distr = [0 for _ in range(101)]

	with contextlib.ExitStack() as stack:
		list_AQD_input_files = [stack.enter_context(open(str_AQD_temp_path_input_file, 'r')) for str_AQD_temp_path_input_file \
			in list_AQD_path_input_files]
		it_AQD_input_fastq = fastq_iterator(fastq_zip_equal(*[FastqGeneralIterator(obj_AQD_temp_input_file) \
			for obj_AQD_temp_input_file in list_AQD_input_files]))


		for list_AQD_reads_for_processing in iter_double_chunked(it_AQD_input_fastq,int_AQD_read_chunks,int_AQD_number_processes):
			with concurrent.futures.ProcessPoolExecutor() as executor:

				list_AQD_read_chunk_distr = executor.map(quality_read_chunk,list_AQD_reads_for_processing,itertools.repeat(int_AQD_phred_cutoff))


				for list_AQD_temp_read_chunk_distr in list_AQD_read_chunk_distr:
					list_AQD_quality_distr = [int_AQD_temp_this_score+int_AQD_temp_total for int_AQD_temp_this_score,int_AQD_temp_total in zip(list_AQD_temp_read_chunk_distr,list_AQD_quality_distr)]

	int_AQD_total_reads = sum(list_AQD_quality_distr)

	int_AQD_reads_pot_filter = 0
	for int_AQD_i in range(101):
		int_AQD_reads_pot_filter += list_AQD_quality_distr[int_AQD_i]
		int_AQD_min_low_q = int_AQD_i
		if division_zero_tolerant(int_AQD_reads_pot_filter,int_AQD_total_reads)*100 >= 100 - int_AQD_max_read_filter_percent:
			break

	if int_AQD_min_low_q > int_AQD_input_max_low_q_percent:
		print('Adjusted filter criterion to ' + str(int_AQD_min_low_q) + '%')
		return int_AQD_min_low_q
	else:
		return int_AQD_input_max_low_q_percent



def quality_filter(dict_QF_read,int_QF_phred_cutoff,float_QF_max_low_q_fraction,int_QF_start_end_trim,int_QF_end_end_trim,float_QF_end_trim_max_low_q_fraction):
	list_QF_low_phreds = [1 if Chr2Phred(str_QF_temp_chr) < int_QF_phred_cutoff else 0 for str_QF_temp_chr in dict_QF_read['quals']]

	int_QF_len_read = len(dict_QF_read['read'])

	bool_QF_read_end_cut = False
	#check if read meets basic quality criteria
	float_QF_low_q_fraction = division_zero_tolerant(sum(list_QF_low_phreds),int_QF_len_read)
	if float_QF_low_q_fraction <= float_QF_max_low_q_fraction:
		bool_QF_read_passed = True
		int_QF_end_cut_point = int_QF_len_read
		#check if quality at read end is low
		int_QF_end_sum_low_q_bases = sum(list_QF_low_phreds[int_QF_start_end_trim:])
		for int_QF_i in range(int_QF_start_end_trim,int_QF_len_read-int_QF_end_end_trim):
			if int_QF_end_sum_low_q_bases/(int_QF_len_read-int_QF_i) > float_QF_end_trim_max_low_q_fraction:
				int_QF_end_cut_point = int_QF_i
				bool_QF_read_end_cut = True
				break
			int_QF_end_sum_low_q_bases -= list_QF_low_phreds[int_QF_i]
		dict_QF_output_read = {str_QF_temp_read_item:dict_QF_read[str_QF_temp_read_item][:int_QF_end_cut_point] for str_QF_temp_read_item in ['read','quals']}
		dict_QF_output_read['name'] = dict_QF_read['name']

	else:
		bool_QF_read_passed = False
		dict_QF_output_read = dict_QF_read

	return [bool_QF_read_passed,bool_QF_read_end_cut,dict_QF_output_read]


def filter_read_pair(list_FRP_reads,int_FRP_phred_cutoff,float_FRP_max_low_q_fraction,int_FRP_start_end_trim,int_FRP_end_end_trim,float_FRP_end_trim_max_low_q_fraction):
	list_FRP_filtered_reads = [quality_filter(dict_FRP_temp_read,int_FRP_phred_cutoff,float_FRP_max_low_q_fraction, \
		int_FRP_start_end_trim,int_FRP_end_end_trim,float_FRP_end_trim_max_low_q_fraction) for dict_FRP_temp_read in list_FRP_reads]

	bool_FRP_reads_passed = all(list_FRP_temp_read[0] is True for list_FRP_temp_read in list_FRP_filtered_reads)
	bool_FRP_reads_cut = any(list_FRP_temp_read[1] is True for list_FRP_temp_read in list_FRP_filtered_reads)

	return [bool_FRP_reads_passed,bool_FRP_reads_cut,[list_FRP_temp_read[2] for list_FRP_temp_read in list_FRP_filtered_reads]]

def process_read_chunk(list_PRC_read_chunk,int_PRC_phred_cutoff,float_PRC_max_low_q_fraction,int_PRC_start_end_trim,int_PRC_end_end_trim,float_PRC_end_trim_max_low_q_fraction):
	list_PRC_processed_reads = []

	for list_PRC_reads in list_PRC_read_chunk:
		list_PRC_processed_reads.append(filter_read_pair(list_PRC_reads,int_PRC_phred_cutoff,float_PRC_max_low_q_fraction,int_PRC_start_end_trim, \
			int_PRC_end_end_trim,float_PRC_end_trim_max_low_q_fraction))

	return list_PRC_processed_reads




def process_files(list_PF_path_input_files,str_PF_path_output_prefix,int_PF_phred_cutoff,int_PF_input_max_low_q_percent,int_PF_max_read_filter_percent,int_PF_start_end_trim, \
	int_PF_end_end_trim,int_PF_end_trim_max_low_q_percent,int_PF_read_chunks=10000,int_PF_number_processes=0):
	#number processes defaults to cpu count
	int_PF_number_processes = int_PF_number_processes if int_PF_number_processes > 0 else os.cpu_count()

	#get start time
	float_PF_time_start = monotonic()
	#get read count
	print(f'Open input file: {list_PF_path_input_files[0]}')
	print(f'Open input file: {list_PF_path_input_files[1]}')
	int_PF_read_count_input_files = fastq_get_read_number(*list_PF_path_input_files)

	print('Analyze quality distribution of input files')		

	int_PF_adjusted_max_low_q_percent = adjust_quality_distr(list_PF_path_input_files,int_PF_phred_cutoff,int_PF_input_max_low_q_percent,int_PF_max_read_filter_percent,int_PF_read_chunks,int_PF_number_processes)

	float_PF_max_low_q_fraction = int_PF_adjusted_max_low_q_percent/100

	if int_PF_adjusted_max_low_q_percent > int_PF_input_max_low_q_percent:
		float_PF_end_trim_max_low_q_fraction = 1
		print('End trimming disabled')
	else:
		float_PF_end_trim_max_low_q_fraction = int_PF_end_trim_max_low_q_percent/100

	list_PF_temp_output_endings, list_PF_temp_output_discarded_endings = ['_qcut_F.fq','_qcut_R.fq'], ['_qcut_discarded_F.fq','_qcut_discarded_R.fq']

	list_PF_path_output_files = [str_PF_path_output_prefix + str_PF_temp_ouput_file_ending \
		for str_PF_temp_ouput_file_ending in list_PF_temp_output_endings]
	list_PF_path_output_discarded_files = [str_PF_path_output_prefix + str_PF_temp_ouput_file_ending \
		for str_PF_temp_ouput_file_ending in list_PF_temp_output_discarded_endings]

	with contextlib.ExitStack() as stack:
		#Input files
		list_PF_input_files = [stack.enter_context(open(str_PF_temp_path_input_file, 'r')) for str_PF_temp_path_input_file \
			in list_PF_path_input_files]
		it_PF_input_fastq = fastq_iterator(fastq_zip_equal(*[FastqGeneralIterator(obj_PF_temp_input_file) \
			for obj_PF_temp_input_file in list_PF_input_files]))

		#output files
		list_PF_output_files = [stack.enter_context(open(str_PF_temp_path_output_file, 'w')) for str_PF_temp_path_output_file \
			in list_PF_path_output_files]
		list_PF_output_discarded_files = [stack.enter_context(open(str_PF_temp_path_output_discarded_file, 'w')) for str_PF_temp_path_output_discarded_file \
			in list_PF_path_output_discarded_files]

		dict_PF_output_switch = {True:list_PF_output_files,False:list_PF_output_discarded_files}

		int_PF_reads_processed = 0
		int_PF_reads_passed_quality_filter = 0
		int_PF_reads_end_trimmed = 0


		for list_PF_reads_for_processing in iter_double_chunked(it_PF_input_fastq,int_PF_read_chunks,int_PF_number_processes):

			with concurrent.futures.ProcessPoolExecutor() as executor:
				list_PF_reads_processed = executor.map(process_read_chunk,list_PF_reads_for_processing,*(itertools.repeat(obj_PF_temp_arg) \
					for obj_PF_temp_arg in (int_PF_phred_cutoff,float_PF_max_low_q_fraction,int_PF_start_end_trim, int_PF_end_end_trim,float_PF_end_trim_max_low_q_fraction)))


				for list_PF_temp_read_output in itertools.chain.from_iterable(list_PF_reads_processed):
					for int_PF_temp_read_FR, dict_PF_temp_read in enumerate(list_PF_temp_read_output[2]):
						dict_PF_output_switch[list_PF_temp_read_output[0]][int_PF_temp_read_FR].write(format_fastq(dict_PF_temp_read))

					if list_PF_temp_read_output[0] is True:
						int_PF_reads_passed_quality_filter += 1
						if list_PF_temp_read_output[1] is True:
							int_PF_reads_end_trimmed += 1
					int_PF_reads_processed += 1

			str_PF_percent_processed = f'{round(int_PF_reads_processed*100/int_PF_read_count_input_files,2):.2f}%' \
				if int_PF_read_count_input_files > 0 else '100.00%'

			sys.stdout.write(f'\rProcessing Reads: {int_PF_reads_processed} / {int_PF_read_count_input_files} ({str_PF_percent_processed})')
			sys.stdout.flush()


	float_PF_time_total = monotonic() - float_PF_time_start

	str_PF_percent_passed_filter = f'{round(int_PF_reads_passed_quality_filter*100/int_PF_read_count_input_files,2):.1f}%' \
		if int_PF_read_count_input_files > 0 else '100.00%'
	str_PF_percent_end_trimmed = f'{round(int_PF_reads_end_trimmed*100/int_PF_reads_passed_quality_filter,2):.1f}%' \
		if int_PF_reads_passed_quality_filter > 0 else '100.00%'

	print(f'\nRead Counts: {int_PF_reads_passed_quality_filter} reads passed quality criterion ({str_PF_percent_passed_filter})')
	print(f'Read Counts: {int_PF_reads_end_trimmed} reads were end trimmed ({str_PF_percent_end_trimmed})')
	print(f'End of files: processed {int_PF_reads_processed} reads in {format_time(float_PF_time_total)}')



def main():

	QFF_argparse = argparse.ArgumentParser(description='Filter NGS-reads by phred quality scores')
	QFF_argparse.add_argument('--input_file1', type=str, required=True, help='Fastq file forward reads')
	QFF_argparse.add_argument('--input_file2', type=str, required=True, help='Fastq file reverse reads')
	QFF_argparse.add_argument('--output_prefix', type=str, required=True, help='Path and file prefix for output files')
	QFF_argparse.add_argument('--phred_cutoff', type=int, required=True, help='lowest phred score not to be considered low quality')
	QFF_argparse.add_argument('--percent_low_q', type=int, required=True, help='percent of bases with low phred score allowed')
	QFF_argparse.add_argument('--percent_max_read_filter', type=int, required=True, help='max percent of reads to be filtered, regardless of quality')
	QFF_argparse.add_argument('--start_end_trim', type=int, required=True, help='min length of read for end trimming')
	QFF_argparse.add_argument('--min_end_trim_length', type=int, required=True, help='min length of read end to be considered for end trimming')
	QFF_argparse.add_argument('--end_trim_percent_low_q', type=int, required=True, help='percent of bases with low phred score allowed in read end')
	QFF_argparse.add_argument('--read_chunks', type=int, required=True, help='reads to be processed at a time')
	QFF_argparse.add_argument('--number_processes', type=int, required=True, help='number of processes to be run in parallel')

	QFF_args = QFF_argparse.parse_args()

	list_QFF_input_files = [QFF_args.input_file1,QFF_args.input_file2]



	process_files(list_QFF_input_files,QFF_args.output_prefix,QFF_args.phred_cutoff,QFF_args.percent_low_q,QFF_args.percent_max_read_filter,QFF_args.start_end_trim, \
		QFF_args.min_end_trim_length,QFF_args.end_trim_percent_low_q,QFF_args.read_chunks,QFF_args.number_processes)


if __name__ == '__main__':
	main()
