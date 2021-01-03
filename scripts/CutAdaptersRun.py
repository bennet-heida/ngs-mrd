#Cut NGS-MRD reads to genomic region and look for random tag for read family analysis
#this script collects data from csv-files and command line, prepares the data and runs 'process_files()' function in CutAdapters.py
import CutAdapters
import csv
from NGSMRDFuncs import format_time
from time import monotonic
import argparse

def main():
	CR_argparse = argparse.ArgumentParser(description='Cut NGS-Reads')
	CR_argparse.add_argument('--run_info_file', type=str, required=True, help='path to run info file')
	CR_argparse.add_argument('--run_name', type=str, required=True, help='name of run')
	CR_argparse.add_argument('--primer_info_file', type=str, required=True, help='path to genomic primer info file')
	CR_argparse.add_argument('--input_folder', type=str, required=True, help='folder with demultiplexed input files')
	CR_argparse.add_argument('--output_folder', type=str, required=True, help='path where output files should be written')
	CR_argparse.add_argument('--maxPen_CS', type=int, required=True, help='number of penalties allowed in CS match')
	CR_argparse.add_argument('--maxPen_gPrimer', type=int, required=True, help='number of penalties allowed in genomic primer match')
	CR_argparse.add_argument('--maxPen_gPrimer_rev', type=int, required=True, help='number of penalties allowed in reverse genomic primer match')
	CR_argparse.add_argument('--buffer_CS', type=int, required=True, help='buffer for search position in CS match')
	CR_argparse.add_argument('--buffer_gPrimer', type=int, required=True, help='buffer for search position in genomic primer match')
	CR_argparse.add_argument('--buffer_gPrimer_rev', type=int, required=True, help='buffer for search position in reverse genomic primer match')
	CR_argparse.add_argument('--max_read_length', type=int, default=241, help='max length of reads')
	CR_argparse.add_argument('--max_primer_length', type=int, default=40, help='max length of primers')
	CR_argparse.add_argument('--read_chunks', type=int, default=10000, help='reads to be processed at a time')
	CR_argparse.add_argument('--number_processes', type=int, default=0, help='number of processes to be run in parallel')

	CR_args = CR_argparse.parse_args()

	str_Common_Sequence_F = "GGTAAACACAAGGGCACTGG"
	str_Common_Sequence_R = "CGGACTACAGCTCCCATCAT"
	list_CR_RF_lenghts = [8,16]


	#read primer info
	dict_CR_gPrimer = {}
	with open(CR_args.primer_info_file,'r') as obj_CR_gPrimer_file:
		it_CR_gPrimer_csv = csv.DictReader(obj_CR_gPrimer_file, delimiter='\t')
		for dict_CR_gPrimer_line in it_CR_gPrimer_csv:
			dict_CR_gPrimer[dict_CR_gPrimer_line['Primer'].strip()] = {'seq_f':dict_CR_gPrimer_line['Seq_1'].strip(), \
				'seq_r':dict_CR_gPrimer_line['Seq_2'].strip(),'gReg_length':abs(int(dict_CR_gPrimer_line['FirstBase_1'])-int(dict_CR_gPrimer_line['FirstBase_2']))+1}

	#get list of samples
	list_CR_samples = []
	with open(CR_args.run_info_file,'r') as obj_CR_run_info_file:
		it_CR_run_info_csv = csv.DictReader(obj_CR_run_info_file, delimiter='\t')

		for dict_CR_run_info_line in it_CR_run_info_csv:
			str_CR_temp_gPrimer_seq_f = dict_CR_gPrimer[dict_CR_run_info_line['Primer']]['seq_f']
			str_CR_temp_gPrimer_seq_r = dict_CR_gPrimer[dict_CR_run_info_line['Primer']]['seq_r']
			int_CR_temp_len_gReg = dict_CR_gPrimer[dict_CR_run_info_line['Primer']]['gReg_length']
			list_CR_variant_change = dict_CR_run_info_line['Variant'].split('>')
			if len(list_CR_variant_change) != 2:
				raise ValueError('Variant not recognized')
			int_CR_length_indel = len(list_CR_variant_change[1]) - len(list_CR_variant_change[0])

			list_CR_samples.append({'sample_file':f'{CR_args.run_name}_{dict_CR_run_info_line["Number"].zfill(2)}' \
				+ f'_{dict_CR_run_info_line["Primer"]}_{dict_CR_run_info_line["Sample"]}', \
				'gPrimer_seq_f': str_CR_temp_gPrimer_seq_f, 'gPrimer_seq_r': str_CR_temp_gPrimer_seq_r, \
				'gReg_length':int_CR_temp_len_gReg,'indel_length':int_CR_length_indel})

	print(f'Process Run: {CR_args.run_name}')

	float_CR_start_time = monotonic()

	#give sample info to 'process_files()' function
	for dict_CR_sample in list_CR_samples:
		list_CR_input_files = [CR_args.input_folder + dict_CR_sample['sample_file'] + '_F.fq', CR_args.input_folder + dict_CR_sample['sample_file'] + '_R.fq']
		str_CR_output_prefix = CR_args.output_folder + dict_CR_sample['sample_file']

		CutAdapters.process_files(list_CR_input_files,str_CR_output_prefix,str_Common_Sequence_F,str_Common_Sequence_R,dict_CR_sample['gPrimer_seq_f'],dict_CR_sample['gPrimer_seq_r'], \
			dict_CR_sample['gReg_length'],dict_CR_sample['indel_length'],CR_args.maxPen_CS,CR_args.maxPen_gPrimer,CR_args.maxPen_gPrimer_rev,CR_args.buffer_CS, \
			CR_args.buffer_gPrimer,CR_args.buffer_gPrimer_rev,list_CR_RF_lenghts,CR_args.max_read_length,CR_args.max_primer_length, \
			CR_args.read_chunks,CR_args.number_processes,)

	print(f'Run processed in: {format_time(monotonic()-float_CR_start_time)}')

if __name__ == '__main__':
	main()