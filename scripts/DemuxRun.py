#Script for NGS-MRD run demultiplexing
#this script collects data from csv-files and command line, prepares the data and runs 'process_files()' function in Demux.py
import csv
import Demux
import argparse

def main():
	DR_argparse = argparse.ArgumentParser(description='Demultiplex NGS-Reads')
	DR_argparse.add_argument('--input_file1', type=str, required=True, help='path to forward read fastq file')
	DR_argparse.add_argument('--input_file2', type=str, required=True, help='path to reverse read fastq file')
	DR_argparse.add_argument('--output_folder', type=str, required=True, help='path where output files should be written')
	DR_argparse.add_argument('--run_name', type=str, required=True, help='name of run')
	DR_argparse.add_argument('--run_info_file', type=str, required=True, help='path to run info file')
	DR_argparse.add_argument('--MID_seq_file', type=str, required=True, help='path to MID sequences file')
	DR_argparse.add_argument('--maxPen_MID', type=int, required=True, help='number of penalties allowed in MID match')
	DR_argparse.add_argument('--check_MID_pos', type=int, choices=[0,1], required=True, help='number of penalties allowed in MID match')
	DR_argparse.add_argument('--read_chunks', type=int, default=10000, help='reads to be processed at a time')
	DR_argparse.add_argument('--number_processes', type=int, default=0, help='number of processes to be run in parallel')

	DR_args = DR_argparse.parse_args()

	list_DR_input_files = [DR_args.input_file1,DR_args.input_file2]
	str_DR_output_folder = DR_args.output_folder + ('' if DR_args.output_folder[-1] == '/' else '/')

	bool_DR_check_MID_pos = False if DR_args.check_MID_pos == 0 else True

	#read MID-info from file
	list_DR_MIDs = [[],[]]
	dict_DR_MID_FR = {'F':0,'R':1}
	with open(DR_args.MID_seq_file,'r') as obj_DR_MID_file:
		it_DR_MID_file_csv = csv.DictReader(obj_DR_MID_file,delimiter='\t')

		for dict_DR_MID in it_DR_MID_file_csv:
			dict_DR_temp_output_MID = {'ID':dict_DR_MID['Primername'],'seq':dict_DR_MID['MID_Sequence'],'start':int(dict_DR_MID['Spacer_Length']),}
			list_DR_MIDs[dict_DR_MID_FR[dict_DR_MID['Stockname'][-1:]]].append(dict_DR_temp_output_MID)

	#get list of valid MID-IDs to check construction of MID-grid in next step
	list_DR_MID_IDs = [[dict_DR_temp_MID['ID'] for dict_DR_temp_MID in list_DR_temp_MID_list] for list_DR_temp_MID_list in list_DR_MIDs]
	
	#read run info file, make list of samples and construct MID-grid
	list_DR_samples = []
	dict_DR_MID_grid = {}
	with open(DR_args.run_info_file,'r') as obj_DR_run_file:
		it_DR_run_file_csv = csv.DictReader(obj_DR_run_file,delimiter='\t')

		for dict_DR_run_sample in it_DR_run_file_csv:
			#read sample names
			list_DR_samples.append([int(dict_DR_run_sample['Number']),f"{str_DR_output_folder}{DR_args.run_name}_{dict_DR_run_sample['Number'].zfill(2)}" \
				+ f"_{dict_DR_run_sample['Primer']}_{dict_DR_run_sample['Sample']}"])

			#create MID-grid
			if all(str_DR_temp_MID in list_DR_MID_IDs[int_DR_n] for int_DR_n,str_DR_temp_MID in enumerate([dict_DR_run_sample['MID_F'],dict_DR_run_sample['MID_R']])):
				if dict_DR_run_sample['MID_F'] not in dict_DR_MID_grid:
					dict_DR_MID_grid[dict_DR_run_sample['MID_F']] = {}
				dict_DR_MID_grid[dict_DR_run_sample['MID_F']][dict_DR_run_sample['MID_R']] = int(dict_DR_run_sample['Number'])
			else:
				raise ValueError(f"MID for sample {dict_DR_run_sample['Number'].zfill(2)}_{dict_DR_run_sample['Primer']}" \
					+ f"_{dict_DR_run_sample['Sample']} not recognized")

	#start demultiplexing
	Demux.process_files(list_DR_input_files,list_DR_samples,str_DR_output_folder+DR_args.run_name,list_DR_MIDs,dict_DR_MID_grid, \
		DR_args.maxPen_MID,bool_DR_check_MID_pos,DR_args.read_chunks,DR_args.number_processes)


if __name__ == '__main__':
	main()