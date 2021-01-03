#Read info from csv files and run process_file() from PileupRF.py
import PileupRF
import csv
import argparse


def main():

	PRFR_argparse = argparse.ArgumentParser(description='Pileup NGS-Reads')
	PRFR_argparse.add_argument('--run_info_file', type=str, required=True, help='path to run info file')
	PRFR_argparse.add_argument('--run_name', type=str, required=True, help='name of run')
	PRFR_argparse.add_argument('--ref_fasta_file', type=str, required=True, help='name of run')
	PRFR_argparse.add_argument('--primer_info_file', type=str, required=True, help='name of run')
	PRFR_argparse.add_argument('--input_folder', type=str, required=True, help='folder with demultiplexed input files')
	PRFR_argparse.add_argument('--output_folder', type=str, required=True, help='path where output files should be written')
	PRFR_argparse.add_argument('--min_member_count', type=int, required=True, help='min number of reads to form read family')
	PRFR_argparse.add_argument('--min_consensus_percent', type=int, required=True, help='min percentage of bases to accept consensus')

	PRFR_args = PRFR_argparse.parse_args()

	bool_PRFR_target_zero_based = False

	str_PRFR_input_folder = PRFR_args.input_folder + ('/' if PRFR_args.input_folder[-1] != '/' else '')
	str_PRFR_output_folder = PRFR_args.output_folder + ('/' if PRFR_args.output_folder[-1] != '/' else '')

	float_PRFR_min_consensus_fraction = PRFR_args.min_consensus_percent / 100


	dict_PRFR_gPrimer = {}

	with open(PRFR_args.primer_info_file,'r') as obj_PRFR_primer_info_file:
		it_PRFR_primer_info_csv = csv.DictReader(obj_PRFR_primer_info_file,delimiter='\t')

		for dict_PRFR_temp_primer in it_PRFR_primer_info_csv:

			int_PRFR_temp_primer_start = min(int(dict_PRFR_temp_primer['FirstBase_1']),int(dict_PRFR_temp_primer['FirstBase_2']))
			int_PRFR_temp_primer_end = max(int(dict_PRFR_temp_primer['FirstBase_1']),int(dict_PRFR_temp_primer['FirstBase_2']))


			dict_PRFR_gPrimer[dict_PRFR_temp_primer['Primer'].strip()] = (dict_PRFR_temp_primer['Chr_1'],int_PRFR_temp_primer_start, \
				int_PRFR_temp_primer_end-int_PRFR_temp_primer_start+1)


	list_PRFR_samples = []

	with open(PRFR_args.run_info_file,'r') as obj_PRFR_run_info_file:
		it_PRFR_run_info_csv = csv.DictReader(obj_PRFR_run_info_file,delimiter='\t')

		for dict_PRFR_sample in it_PRFR_run_info_csv:

			str_PRFR_base_sample_name = f'{PRFR_args.run_name}_{dict_PRFR_sample["Number"].zfill(2)}_{dict_PRFR_sample["Primer"]}_{dict_PRFR_sample["Sample"]}'

			list_PRFR_samples.append({'input_file':str_PRFR_input_folder+str_PRFR_base_sample_name+'_RF.bam','output_file':str_PRFR_output_folder \
				+str_PRFR_base_sample_name+'_RF.tsv','gReg_chr':dict_PRFR_gPrimer[dict_PRFR_sample['Primer']][0],'gReg_start':dict_PRFR_gPrimer[dict_PRFR_sample['Primer']][1], \
				'gReg_length':dict_PRFR_gPrimer[dict_PRFR_sample['Primer']][2]})



	for dict_PRFR_temp_sample in list_PRFR_samples:
		PileupRF.process_file(dict_PRFR_temp_sample['input_file'],dict_PRFR_temp_sample['output_file'],PRFR_args.ref_fasta_file,dict_PRFR_temp_sample['gReg_chr'], \
			dict_PRFR_temp_sample['gReg_start'],dict_PRFR_temp_sample['gReg_length'],bool_PRFR_target_zero_based,PRFR_args.min_member_count,float_PRFR_min_consensus_fraction)

if __name__ == '__main__':
	main()