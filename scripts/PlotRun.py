#read run info from csv files and plot target-VAF and background error
#runs plot_pileup_SNP_LVAF, plot_pileup_SNP_TVAF or plot_pileup_indel
#from PlotRun.py based on mutation type

import csv
from NGSMRDFuncs import get_analysis_range
import PlotMRD
import argparse



def main():
	PR_argparse = argparse.ArgumentParser(description='Output MRD scores to Excel file')
	PR_argparse.add_argument('--run_name', type=str, required=True, help='run name')
	PR_argparse.add_argument('--run_info_file', type=str, required=True, help='Path to run file')
	PR_argparse.add_argument('--MID_file', type=str, required=True, help='Path to MID sequence file')
	PR_argparse.add_argument('--gPrimer_file', type=str, required=True, help='Path to gPrimer sequence file')
	PR_argparse.add_argument('--pileup_folder', type=str, required=True, help='Path to folder, where pileup files are stored')
	PR_argparse.add_argument('--output_folder', type=str, required=True, help='Path to output Excel file')

	PR_args = PR_argparse.parse_args()

	str_PR_pileup_folder = PR_args.pileup_folder + ('/' if PR_args.pileup_folder[-1] != '/' else '')
	str_PR_output_folder = PR_args.output_folder + ('/' if PR_args.output_folder[-1] != '/' else '')

	float_PR_LOD_stdev_test = 3
	float_PR_LOD_stdev_cutoff_BE = 2.5
	int_PR_NOP_target_range = 25
	int_PR_R1_Nsupp_cutoff = 75

	int_PR_length_illumina_read = 250
	int_PR_length_CS = 20
	int_PR_length_RF = 16
	int_PR_buffer_to_read_end = 15

	list_PR_SNP_pileup_file_ends = ['_R1R2.tsv','_RF.tsv']
	list_PR_SNP_output_file_ends = ['_R1R2','_RF']

	list_PR_indel_pileup_file_ends = ['_R1R2.tsv']
	list_PR_indel_output_file_ends = ['_R1R2']

	list_PR_samples = []

	with open(PR_args.run_info_file,'r') as obj_PR_run_file:
		it_PR_run_csv = csv.DictReader(obj_PR_run_file,delimiter='\t')

		for dict_PR_temp_sample in it_PR_run_csv:
			str_PR_temp_file_base = f'{PR_args.run_name}_{dict_PR_temp_sample["Number"].zfill(2)}_{dict_PR_temp_sample["Primer"]}_{dict_PR_temp_sample["Sample"]}'

			int_PR_target_position = int(dict_PR_temp_sample['Coordinate'])

			list_PR_analysis_range = get_analysis_range(dict_PR_temp_sample['Primer'],dict_PR_temp_sample['MID_F'],dict_PR_temp_sample['MID_R'],PR_args.MID_file, \
				PR_args.gPrimer_file,int_PR_length_illumina_read,int_PR_length_CS,int_PR_length_RF,int_PR_buffer_to_read_end)

			if list_PR_analysis_range[0] == dict_PR_temp_sample['Chr'] and (list_PR_analysis_range[1][0] <= int_PR_target_position <= list_PR_analysis_range[1][1]):

				list_PR_base_change = dict_PR_temp_sample['Variant'].split('>')

				if len(list_PR_base_change) != 2 or 0 in map(len,list_PR_base_change):
					raise ValueError('Base change not recognized')

				int_PR_indel_length = abs(len(list_PR_base_change[1])-len(list_PR_base_change[0]))

				if int_PR_indel_length == 0:

					list_PR_temp_pileup_files = [str_PR_pileup_folder + str_PR_temp_file_base + str_PR_temp_file_end for str_PR_temp_file_end in list_PR_SNP_pileup_file_ends]
					list_PR_temp_output_files = [str_PR_output_folder + str_PR_temp_file_base + str_PR_temp_file_end for str_PR_temp_file_end in list_PR_SNP_output_file_ends]

				else:
					list_PR_temp_pileup_files = [str_PR_pileup_folder + str_PR_temp_file_base + str_PR_temp_file_end for str_PR_temp_file_end in list_PR_indel_pileup_file_ends]
					list_PR_temp_output_files = [str_PR_output_folder + str_PR_temp_file_base + str_PR_temp_file_end for str_PR_temp_file_end in list_PR_indel_output_file_ends]

				list_PR_samples.append({'indel_length':int_PR_indel_length,'pileup_files':list_PR_temp_pileup_files,'output_files':list_PR_temp_output_files,'target':{'chr':dict_PR_temp_sample['Chr'], \
					'position':int_PR_target_position,'variant':list_PR_base_change},'bounds':list_PR_analysis_range[1]})

			else:
				print(f'Sample {dict_PR_temp_sample["Number"]} target: {dict_PR_temp_sample["Chr"]}:{str(int_PR_target_position)} not on primer region ' + \
					f'({list_PR_analysis_range[0]}:{str(list_PR_analysis_range[1][0])}-{str(list_PR_analysis_range[1][1])})')


	for dict_PR_plot_sample in list_PR_samples:
		if dict_PR_plot_sample['indel_length'] == 0:

			for int_PR_i in range(2):

				PlotMRD.plot_pileup_SNP_LVAF(dict_PR_plot_sample['pileup_files'][int_PR_i],dict_PR_plot_sample['output_files'][int_PR_i]+'_LVAF.svg',dict_PR_plot_sample['target']['chr'],dict_PR_plot_sample['bounds'][0], \
					dict_PR_plot_sample['bounds'][1]-dict_PR_plot_sample['bounds'][0]+1,dict_PR_plot_sample['target']['chr'],dict_PR_plot_sample['target']['position'],dict_PR_plot_sample['target']['variant'], \
					float_PR_LOD_stdev_test,float_PR_LOD_stdev_cutoff_BE)
			
				PlotMRD.plot_pileup_SNP_TVAF(dict_PR_plot_sample['pileup_files'][int_PR_i],dict_PR_plot_sample['output_files'][int_PR_i]+'_TVAF.svg',dict_PR_plot_sample['target']['chr'],dict_PR_plot_sample['bounds'][0], \
					dict_PR_plot_sample['bounds'][1]-dict_PR_plot_sample['bounds'][0]+1,dict_PR_plot_sample['target']['chr'],dict_PR_plot_sample['target']['position'],dict_PR_plot_sample['target']['variant'], \
					float_PR_LOD_stdev_test,float_PR_LOD_stdev_cutoff_BE)

		else:
			PlotMRD.plot_pileup_indel(dict_PR_plot_sample['pileup_files'][0],dict_PR_plot_sample['output_files'][0]+'_indel.svg',dict_PR_plot_sample['target']['chr'],dict_PR_plot_sample['bounds'][0], \
				dict_PR_plot_sample['bounds'][1]-dict_PR_plot_sample['bounds'][0]+1,dict_PR_plot_sample['target']['chr'],dict_PR_plot_sample['target']['position'],dict_PR_plot_sample['target']['variant'], \
				float_PR_LOD_stdev_test,float_PR_LOD_stdev_cutoff_BE)


if __name__ == '__main__':
	main()