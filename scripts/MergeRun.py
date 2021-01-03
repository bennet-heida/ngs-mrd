#Script for NGS-MRD run read merging
#this script collects data from csv-files and command line, prepares the data and runs 'process_files()' function in MergeReads.py
import csv
import MergeReads
import argparse
from time import monotonic
from NGSMRDFuncs import format_time


def main():

	MR_argparse = argparse.ArgumentParser(description='Merge NGS-Reads')
	MR_argparse.add_argument('--run_info_file', type=str, required=True, help='path to run info file')
	MR_argparse.add_argument('--run_name', type=str, required=True, help='name of run')
	MR_argparse.add_argument('--MID_seq_file', type=str, required=True, help='path to MID primer info file')
	MR_argparse.add_argument('--primer_info_file', type=str, required=True, help='path to genomic primer info file')
	MR_argparse.add_argument('--input_folder', type=str, required=True, help='folder with demultiplexed and cut input files')
	MR_argparse.add_argument('--output_folder', type=str, required=True, help='path where output files should be written')
	MR_argparse.add_argument('--min_consensus_length', type=int, required=True, help='minimum length of consensus read')
	MR_argparse.add_argument('--min_score_fraction', type=int, required=True, help='score of alignment must reach this percentage of the reads length')
	MR_argparse.add_argument('--read_chunks', type=int, default=10000, help='reads to be processed at a time')
	MR_argparse.add_argument('--number_processes', type=int, default=0, help='number of processes to be run in parallel')

	MR_args = MR_argparse.parse_args()

	int_MR_min_length = MR_args.min_consensus_length #min length of output consensus string
	float_MR_min_ScoreFraction = MR_args.min_score_fraction/100 #alignment score must be higher than length of output consensus string times ScoreFraction
	int_MR_length_illumina_read = 250
	int_MR_length_CS = 20
	int_MR_length_RF = 16
	int_MR_buffer_to_read_end = 5
	int_MR_matchScore = 1
	int_MR_cutoff_long_primer = 170
	int_MR_arraySize = 240 #max length of consensus alignment input

	#read MID-info from file
	list_MR_MIDs = [{},{}]
	dict_MR_MID_FR = {'F':0,'R':1}
	with open(MR_args.MID_seq_file,'r') as obj_MR_MID_file:
		it_MR_MID_file_csv = csv.DictReader(obj_MR_MID_file,delimiter='\t')
		for dict_MR_temp_MID in it_MR_MID_file_csv:
			int_MR_temp_FR = dict_MR_MID_FR[dict_MR_temp_MID['Stockname'][-1]]
			list_MR_MIDs[int_MR_temp_FR][dict_MR_temp_MID['Primername']] = int(dict_MR_temp_MID['Spacer_Length']) + len(dict_MR_temp_MID['MID_Sequence'])

	#read primer info
	dict_MR_gPrimer = {}
	with open(MR_args.primer_info_file,'r') as obj_MR_primer_file:
		it_MR_primer_file_csv = csv.DictReader(obj_MR_primer_file,delimiter='\t')
		for dict_MR_temp_primer in it_MR_primer_file_csv:
			dict_MR_gPrimer[dict_MR_temp_primer['Primer'].strip()] = {'gLength':abs(int(dict_MR_temp_primer['FirstBase_1']) - int(dict_MR_temp_primer['FirstBase_2']))+1, \
				'Seq_Lengths':[len(dict_MR_temp_primer['Seq_1']),len(dict_MR_temp_primer['Seq_2'])]}
			

	#read sample info
	list_MR_samples = []
	list_MR_temp_input_file_endings = ['_qcut_F.fq','_qcut_R.fq']
	with open(MR_args.run_info_file,'r') as obj_MR_run_info_file:
		it_MR_run_info_file_csv = csv.DictReader(obj_MR_run_info_file,delimiter='\t')
		for dict_MR_temp_sample in it_MR_run_info_file_csv:

			list_MR_temp_path_input_files = [f"{MR_args.input_folder}{MR_args.run_name}_{dict_MR_temp_sample['Number'].zfill(2)}_{dict_MR_temp_sample['Primer']}" \
				+ f"_{dict_MR_temp_sample['Sample']}" + str_MR_temp_end for str_MR_temp_end in list_MR_temp_input_file_endings]
			str_MR_temp_path_output_file = f"{MR_args.output_folder}{MR_args.run_name}_{dict_MR_temp_sample['Number'].zfill(2)}_{dict_MR_temp_sample['Primer']}" \
				+ f"_{dict_MR_temp_sample['Sample']}_R1R2.fq"

			bool_MR_temp_cut_reads = True if dict_MR_gPrimer[dict_MR_temp_sample['Primer']]['gLength'] + sum(dict_MR_gPrimer[dict_MR_temp_sample['Primer']]['Seq_Lengths']) > int_MR_cutoff_long_primer else False

			if bool_MR_temp_cut_reads is True:
				int_MR_temp_gReg_Length = dict_MR_gPrimer[dict_MR_temp_sample['Primer']]['gLength']
				list_MR_temp_max_read_lengths = [int_MR_length_illumina_read-int_MR_temp_MID_length-int_MR_length_CS-int_MR_length_RF-int_MR_temp_gPrimer_seq_length-int_MR_buffer_to_read_end \
					for int_MR_temp_MID_length,int_MR_temp_gPrimer_seq_length in zip([list_MR_MIDs[0][dict_MR_temp_sample['MID_F']],list_MR_MIDs[1][dict_MR_temp_sample['MID_R']]], \
					dict_MR_gPrimer[dict_MR_temp_sample['Primer']]['Seq_Lengths'])]
			else:
				int_MR_temp_gReg_Length = 0
				list_MR_temp_max_read_lengths = [0,0]

			list_MR_samples.append({'input_files':list_MR_temp_path_input_files,'output_file':str_MR_temp_path_output_file, \
				'cut_reads':bool_MR_temp_cut_reads,'gReg_Length':int_MR_temp_gReg_Length,'max_read_length':list_MR_temp_max_read_lengths})


	float_MR_start_time = monotonic()
	print(f'Process run: {MR_args.run_name}')

	for dict_MR_sample_for_processing in list_MR_samples:
		MergeReads.process_files(dict_MR_sample_for_processing['input_files'],dict_MR_sample_for_processing['output_file'], \
			int_MR_min_length,float_MR_min_ScoreFraction,dict_MR_sample_for_processing['cut_reads'],dict_MR_sample_for_processing['gReg_Length'], \
			dict_MR_sample_for_processing['max_read_length'],int_MR_matchScore,int_MR_arraySize,MR_args.read_chunks,MR_args.number_processes)


	print(f'Run processed in: {format_time(monotonic()-float_MR_start_time)}')

if __name__ == '__main__':
	main()