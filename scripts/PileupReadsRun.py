#Read run info from file and run process_file() from PileupReads.py
import PileupReads
import csv
import argparse


def main():

	PRR_argparse = argparse.ArgumentParser(description='Pileup NGS-Reads')
	PRR_argparse.add_argument('--run_info_file', type=str, required=True, help='path to run info file')
	PRR_argparse.add_argument('--run_name', type=str, required=True, help='name of run')
	PRR_argparse.add_argument('--ref_fasta_file', type=str, required=True, help='name of run')
	PRR_argparse.add_argument('--primer_info_file', type=str, required=True, help='name of run')
	PRR_argparse.add_argument('--input_folder', type=str, required=True, help='folder with demultiplexed input files')
	PRR_argparse.add_argument('--output_folder', type=str, required=True, help='path where output files should be written')

	PRR_args = PRR_argparse.parse_args()

	bool_PRR_target_zero_based = False

	str_PRR_input_folder = PRR_args.input_folder + ('/' if PRR_args.input_folder[-1] != '/' else '')
	str_PRR_output_folder = PRR_args.output_folder + ('/' if PRR_args.output_folder[-1] != '/' else '')

	dict_PRR_gPrimer = {}

	with open(PRR_args.primer_info_file,'r') as obj_PRR_primer_info_file:
		it_PRR_primer_info_csv = csv.DictReader(obj_PRR_primer_info_file,delimiter='\t')

		for dict_PRR_temp_primer in it_PRR_primer_info_csv:

			int_PRR_temp_primer_start = min(int(dict_PRR_temp_primer['FirstBase_1']),int(dict_PRR_temp_primer['FirstBase_2']))
			int_PRR_temp_primer_end = max(int(dict_PRR_temp_primer['FirstBase_1']),int(dict_PRR_temp_primer['FirstBase_2']))


			dict_PRR_gPrimer[dict_PRR_temp_primer['Primer'].strip()] = (dict_PRR_temp_primer['Chr_1'],int_PRR_temp_primer_start, \
				int_PRR_temp_primer_end-int_PRR_temp_primer_start+1)


	list_PRR_samples = []

	with open(PRR_args.run_info_file,'r') as obj_PRR_run_info_file:
		it_PRR_run_info_csv = csv.DictReader(obj_PRR_run_info_file,delimiter='\t')

		for dict_PRR_sample in it_PRR_run_info_csv:

			str_PRR_base_sample_name = f'{PRR_args.run_name}_{dict_PRR_sample["Number"].zfill(2)}_{dict_PRR_sample["Primer"]}_{dict_PRR_sample["Sample"]}'

			for str_PRR_temp_input_file_ending in ['_R1','_R2','_R1R2']:

				list_PRR_samples.append({'input_file':str_PRR_input_folder+str_PRR_base_sample_name+str_PRR_temp_input_file_ending+'.bam','output_file':str_PRR_output_folder \
					+str_PRR_base_sample_name+str_PRR_temp_input_file_ending+'.tsv','gReg_chr':dict_PRR_gPrimer[dict_PRR_sample['Primer']][0],'gReg_start':dict_PRR_gPrimer[dict_PRR_sample['Primer']][1], \
					'gReg_length':dict_PRR_gPrimer[dict_PRR_sample['Primer']][2]})



	for dict_PRR_temp_sample in list_PRR_samples:
		PileupReads.process_file(dict_PRR_temp_sample['input_file'],dict_PRR_temp_sample['output_file'],PRR_args.ref_fasta_file,dict_PRR_temp_sample['gReg_chr'], \
			dict_PRR_temp_sample['gReg_start'],dict_PRR_temp_sample['gReg_length'],bool_PRR_target_zero_based)



if __name__ == '__main__':
	main()
