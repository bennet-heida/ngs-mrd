#Filter reads from NGS-MRD run by phred quality scores
#this script collects data from csv-files and command line, prepares the data and runs 'process_files()' function in QFilter.py
import csv
import QFilter
import argparse
from time import monotonic
from NGSMRDFuncs import format_time


def main():

	QFR_argparse = argparse.ArgumentParser(description='Filter NGS-reads by phred quality scores')
	QFR_argparse.add_argument('--run_name', type=str, required=True, help='Fastq file forward reads')
	QFR_argparse.add_argument('--run_info_file', type=str, required=True, help='Fastq file forward reads')
	QFR_argparse.add_argument('--input_folder', type=str, required=True, help='Path and file prefix for output files')
	QFR_argparse.add_argument('--output_folder', type=str, required=True, help='Path and file prefix for output files')
	QFR_argparse.add_argument('--phred_cutoff', type=int, required=True, help='lowest phred score not to be considered low quality')
	QFR_argparse.add_argument('--percent_low_q', type=int, required=True, help='percent of bases with low phred score allowed')
	QFR_argparse.add_argument('--percent_max_read_filter', type=int, required=True, help='max percent of reads to be filtered, regardless of quality')
	QFR_argparse.add_argument('--start_end_trim', type=int, required=True, help='min length of read for end trimming')
	QFR_argparse.add_argument('--min_end_trim_length', type=int, required=True, help='min length of read end to be considered for end trimming')
	QFR_argparse.add_argument('--end_trim_percent_low_q', type=int, required=True, help='percent of bases with low phred score allowed in read end')
	QFR_argparse.add_argument('--read_chunks', type=int, required=True, help='reads to be processed at a time')
	QFR_argparse.add_argument('--number_processes', type=int, default=0, help='number of processes to be run in parallel')

	QFR_args = QFR_argparse.parse_args()

	float_QFR_start_time = monotonic()

	print(f'Process run: {QFR_args.run_name}')


	str_QFR_input_folder = QFR_args.input_folder + ('/' if QFR_args.input_folder[-1] != '/' else '')
	str_QFR_output_folder = QFR_args.output_folder + ('/' if QFR_args.output_folder[-1] != '/' else '')

	#read sample info
	list_QFR_samples = []
	with open(QFR_args.run_info_file,'r') as obj_QFR_run_info_file:
		it_QFR_run_info_file_csv = csv.DictReader(obj_QFR_run_info_file,delimiter='\t')
		for dict_QFR_temp_sample in it_QFR_run_info_file_csv:
			str_QFR_sample_name_prefix = f'{QFR_args.run_name}_{dict_QFR_temp_sample["Number"].zfill(2)}' \
				f'_{dict_QFR_temp_sample["Primer"]}_{dict_QFR_temp_sample["Sample"]}'

			str_QFR_temp_input_sample_prefix = str_QFR_input_folder + str_QFR_sample_name_prefix + '_cut'
			str_QFR_temp_output_sample_prefix = str_QFR_output_folder + str_QFR_sample_name_prefix

			list_QFR_samples.append([[str_QFR_temp_input_sample_prefix + str_QFR_temp_input_file_ending \
				for str_QFR_temp_input_file_ending in ['_F.fq', '_R.fq']],str_QFR_temp_output_sample_prefix])


	for list_QFR_temp_sample in list_QFR_samples:

		QFilter.process_files(list_QFR_temp_sample[0],list_QFR_temp_sample[1],QFR_args.phred_cutoff,QFR_args.percent_low_q,QFR_args.percent_max_read_filter,QFR_args.start_end_trim, \
			QFR_args.min_end_trim_length,QFR_args.end_trim_percent_low_q,QFR_args.read_chunks,QFR_args.number_processes)


	print(f'Run processed in: {format_time(monotonic()-float_QFR_start_time)}')

if __name__ == '__main__':
	main()