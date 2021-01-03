#Create results table from MRD pileup data
#script applies logic as described in thesis:
#SNP: if number of RF > 10000, R1/R2 and RF is used, if RF < 10000 only R1/R2 is used
#if VAF > av + 3sd and no other peak > VAF - (av + 3sd) then SNP is called positive
#small indel (<=3bp): only R1/R2 analysis, if VAF > av + 3sd then indel is called positive
#large indel (>3bp): only R1 or R2 seperate reads are used, if number of supporting reads
#exceeds >75 (NPM1: >10) then indel is called positive

import csv
import CalcMRDScore
import argparse
import xlsxwriter
from NGSMRDFuncs import fastq_get_read_number,get_RF_distr,division_zero_tolerant,get_analysis_range



def format_mrd_output(list_FMO_MRD_calc):
	#output format: Nsupp, Rdepth, VAF, [BEpos, av, stdev, av+3*stdev, T1, VAF-3*stdev, NOP, Test]

	list_FMO_tests = [[dict_FMO_temp_test['be'][0][0],dict_FMO_temp_test['be'][1][0]*100,dict_FMO_temp_test['be'][1][1]*100,dict_FMO_temp_test['lod'][0]*100, \
		1 if dict_FMO_temp_test['t'][0] else 0,dict_FMO_temp_test['lod'][1]*100,dict_FMO_temp_test['nop'], \
		1 if dict_FMO_temp_test['t'][2] else 0] for dict_FMO_temp_test in list_FMO_MRD_calc[3]]

	if len(list_FMO_tests) < 2:
		list_FMO_tests.append([None for _ in range(8)])

	list_FMO_MRD_output = list_FMO_MRD_calc[0:2] + [list_FMO_MRD_calc[2]*100] + [obj_FMO_item for list_FMO_sublist in list_FMO_tests for obj_FMO_item in list_FMO_sublist]

	return list_FMO_MRD_output

def get_mrd_scores(list_GMS_samples,float_GMS_LOD_stdev_test,float_GMS_LOD_stdev_cutoff_BE,int_GMS_NOP_target_range,int_GMS_R1_Nsupp_cutoff):
	list_GMS_R1_tests = ['_R1.tsv','_R2.tsv']
	list_GMS_EC_tests = ['_R1R2.tsv','_RF.tsv']

	list_GMS_MRD_scores = []

	for dict_GMS_temp_sample in list_GMS_samples:

		list_GMS_temp_final_line = []

		#general sample info
		list_GMS_sample_ident = [dict_GMS_temp_sample['number'],dict_GMS_temp_sample['sample'],dict_GMS_temp_sample['primer'],f'{dict_GMS_temp_sample["target"]["chr"]}:' \
			+ f'{dict_GMS_temp_sample["target"]["position"]}:{">".join(dict_GMS_temp_sample["target"]["base_change"])}']

		#check if primer on target region
		if dict_GMS_temp_sample['target_region']['bounds'][0] <= dict_GMS_temp_sample['target']['position'] <= dict_GMS_temp_sample['target_region']['bounds'][1] \
			and dict_GMS_temp_sample['target']['chr'] == dict_GMS_temp_sample['target_region']['chr']:

			#tests to run
			#get MRD scores
			list_GMS_EC_results = [CalcMRDScore.calc_mrd_score(dict_GMS_temp_sample['base_filename']+str_GMS_temp_file_end,dict_GMS_temp_sample['target']['chr'], \
				dict_GMS_temp_sample['target_region']['bounds'][0],dict_GMS_temp_sample['target_region']['bounds'][1]-dict_GMS_temp_sample['target_region']['bounds'][0]+1,dict_GMS_temp_sample['target']['position'], \
				dict_GMS_temp_sample['target']['base_change'],float_GMS_LOD_stdev_test,float_GMS_LOD_stdev_cutoff_BE,int_GMS_NOP_target_range) for str_GMS_temp_file_end in list_GMS_EC_tests]

			list_GMS_R1_temp_results = [CalcMRDScore.calc_mrd_score(dict_GMS_temp_sample['base_filename']+str_GMS_temp_file_end,dict_GMS_temp_sample['target']['chr'], \
				dict_GMS_temp_sample['target_region']['bounds'][0],dict_GMS_temp_sample['target_region']['bounds'][1]-dict_GMS_temp_sample['target_region']['bounds'][0]+1,dict_GMS_temp_sample['target']['position'], \
				dict_GMS_temp_sample['target']['base_change'],float_GMS_LOD_stdev_test,float_GMS_LOD_stdev_cutoff_BE,int_GMS_NOP_target_range) for str_GMS_temp_file_end in list_GMS_R1_tests]
			#check if R1 or R2 has more Nsupp
			if list_GMS_R1_temp_results[0][0] >= list_GMS_R1_temp_results[1][0]:
				list_GMS_R1_result = list_GMS_R1_temp_results[0]
				str_GMS_large_indel_test = 'R1'
			else:
				list_GMS_R1_result = list_GMS_R1_temp_results[1]
				str_GMS_large_indel_test = 'R2'

			list_GMS_test_results = list_GMS_EC_results + [list_GMS_R1_result]

			list_GMS_test_results_formatted = [format_mrd_output(list_GMS_temp_test_result) for list_GMS_temp_test_result in list_GMS_test_results]

			#add R1/R2 and RF results
			for list_GMS_temp_result_full_test in list_GMS_test_results_formatted[:2]:
				list_GMS_temp_final_line += list_GMS_sample_ident + list_GMS_temp_result_full_test + ['#']
			#add R1 results
			list_GMS_temp_final_line += list_GMS_sample_ident + list_GMS_test_results_formatted[2][:3] + ['#',None]

			##add final test columns
			#indel-length
			int_GMS_target_indel_length = max(map(len,dict_GMS_temp_sample['target']['base_change']))-1

			list_GSM_final_decision = final_decision(list_GMS_test_results,int_GMS_target_indel_length,str_GMS_large_indel_test,int_GMS_R1_Nsupp_cutoff)

			list_GMS_temp_final_line += list_GSM_final_decision + ['#',None]

			#add additional info, RF distributions
			#get read number
			int_GMS_temp_fastq_read_number = fastq_get_read_number(dict_GMS_temp_sample['fastq_file'])

			list_GMS_temp_RF_distr = get_RF_distr(dict_GMS_temp_sample['bam_file_RF'],3,0,32,2)
			str_GMS_temp_RF_distr_join = " ".join([":".join(str(int_GMS_item) for int_GMS_item in tuple_GMS_item) for tuple_GMS_item in list_GMS_temp_RF_distr[3]])


			list_GMS_temp_additional = [f'{dict_GMS_temp_sample["target_region"]["chr"]}:{dict_GMS_temp_sample["target_region"]["bounds"][0]}-{dict_GMS_temp_sample["target_region"]["bounds"][1]}', int_GMS_temp_fastq_read_number, \
				division_zero_tolerant(list_GMS_test_results[0][1],int_GMS_temp_fastq_read_number),list_GMS_temp_RF_distr[0],list_GMS_temp_RF_distr[1],division_zero_tolerant(list_GMS_temp_RF_distr[0], \
				list_GMS_temp_RF_distr[1]),list_GMS_temp_RF_distr[2], division_zero_tolerant(list_GMS_temp_RF_distr[1],list_GMS_temp_RF_distr[2]),str_GMS_temp_RF_distr_join]

			list_GMS_temp_final_line += list_GMS_temp_additional


		#target not on primer region
		else:
			list_GMS_temp_final_line = list_GMS_sample_ident + ['Target not on primer region']

		list_GMS_MRD_scores.append(list_GMS_temp_final_line)

	return list_GMS_MRD_scores


def final_decision(list_FD_test_results,int_FD_indel_length,str_FD_large_indel_test,int_FD_R1_Nsupp_cutoff):
	dict_FD_test_order = {'R1R2':0,'RF':1,'R1':2}

	#SNP case
	if int_FD_indel_length == 0:
		list_FD_output = RF_decision(list_FD_test_results)
	elif 1 <= int_FD_indel_length <= 3:
		list_FD_small_indel_test_to_use = list_FD_test_results[dict_FD_test_order['R1R2']]
		int_FD_small_indel_algorithm = 3
		str_FD_small_indel_decision = 'R1R2: ' + ('VAF>3SD' if list_FD_small_indel_test_to_use[3][0]['t'][0] else 'VAF<3SD')
		list_FD_output = [list_FD_small_indel_test_to_use[2],1 if list_FD_small_indel_test_to_use[3][0]['t'][0] else 0,int_FD_small_indel_algorithm,str_FD_small_indel_decision]
	else:
		list_FD_large_indel_test_to_use = list_FD_test_results[dict_FD_test_order['R1']]
		int_FD_large_indel_algorithm = 4
		if list_FD_large_indel_test_to_use[0] >= int_FD_R1_Nsupp_cutoff:
			int_FD_large_indel_test = 1
			str_FD_large_indel_decision = str_FD_large_indel_test + ': Nsupp>75'
		else:
			int_FD_large_indel_test = 0
			str_FD_large_indel_decision = str_FD_large_indel_test + ': Nsupp<75'
		
		list_FD_output = [list_FD_large_indel_test_to_use[2],int_FD_large_indel_test,int_FD_large_indel_algorithm,str_FD_large_indel_decision]

	list_FD_output = [list_FD_output[0]*100] + list_FD_output[1:]

	return list_FD_output

def RF_decision(list_RD_formal_tests):
	dict_RD_test_order = {'R1R2':0,'RF':1,'R1':2}
	int_RD_RF_cutoff = 10000

	if list_RD_formal_tests[dict_RD_test_order['RF']][1] >= int_RD_RF_cutoff:

		list_RD_final_results = [list_RD_formal_tests[dict_RD_test_order[str_RD_temp_test]][3][0]['t'][2] for str_RD_temp_test in ['R1R2','RF']]

		#when both negative, use RF results
		if not any(list_RD_final_results):
			int_RD_algorithm = 2
			list_RD_results_to_use = list_RD_formal_tests[dict_RD_test_order['RF']]
			str_RD_decision = 'RF: '
		elif list_RD_final_results[1] is True:
			int_RD_algorithm = 2
			list_RD_results_to_use = list_RD_formal_tests[dict_RD_test_order['RF']]
			str_RD_decision = 'RF: '
		else:
			int_RD_algorithm = 1
			list_RD_results_to_use = list_RD_formal_tests[dict_RD_test_order['R1R2']]
			str_RD_decision = 'R1R2 (RFneg): '

	else:
		int_RD_algorithm = 1
		list_RD_results_to_use = list_RD_formal_tests[dict_RD_test_order['R1R2']]
		str_RD_decision = 'R1R2 (RF<10000): '

	str_RD_decision += snp_decision_str(list_RD_results_to_use[3][0]['t'])

	list_RF_output = [list_RD_results_to_use[2],1 if list_RD_results_to_use[3][0]['t'][2] else 0,int_RD_algorithm,str_RD_decision]

	return list_RF_output

def snp_decision_str(list_SDS_formal_tests):
	str_SDS_output = ''
	if list_SDS_formal_tests[0] is True:
		str_SDS_output += 'VAF > 3SD'
		if list_SDS_formal_tests[1] is True:
			str_SDS_output += ' and NOP=0'
		else:
			str_SDS_output += ' but NOP>0'
	else:
		str_SDS_output += 'VAF < 3SD'

	return str_SDS_output

def write_output_file(str_WOF_output_file,list_WOF_MRD_scores):

	list_WOF_head2_sample = ['Number', 'Sample', 'Primer', 'Target']
	list_WOF_head2_tests1 = ['Nsupp', 'Rdepth', 'VAF']
	list_WOF_head2_tests2 = ['BEpos', 'av', 'stdev', 'av+3*stdev', 'T1', 'VAF-3*stdev', 'NOP', 'Test']
	list_WOF_head2_final_decision = ['Final VAF', 'Final Test', 'Algorithm', 'Decision']
	list_WOF_head2_addition = ['Pileup_region','fastq','aln','RF_incl','RF_reads','RF_incl / RF_reads','RF_N','RF_reads / RF_N','RF_distr']

	list_WOF_head2_complete = (list_WOF_head2_sample + list_WOF_head2_tests1 + 2*list_WOF_head2_tests2 + [None])*2 + list_WOF_head2_sample + list_WOF_head2_tests1 + [None,None] \
		+ list_WOF_head2_final_decision + [None,None] + list_WOF_head2_addition

	list_WOF_head1_complete = []

	list_WOF_head1_LVAF_TVAF = []
	for str_WOF_VAF in ['LVAF','TVAF']:
		list_WOF_head1_LVAF_TVAF += [str_WOF_VAF] + [None]*(len(list_WOF_head2_tests2)-1)

	for str_WOF_temp_full_test in ['R1/R2','RF']:
		list_WOF_head1_complete += [None]*len(list_WOF_head2_sample) + [str_WOF_temp_full_test] + [None]*(len(list_WOF_head2_tests1)-1) + list_WOF_head1_LVAF_TVAF + [None]

	list_WOF_head1_complete += [None]*len(list_WOF_head2_sample) + ['R1'] + [None]*(len(list_WOF_head2_tests1)-1) + [None]*2


	with xlsxwriter.Workbook(str_WOF_output_file) as obj_WOF_workbook:
		obj_WOF_worksheet = obj_WOF_workbook.add_worksheet()

		obj_WOF_worksheet.write_row(0,0,list_WOF_head1_complete)
		obj_WOF_worksheet.write_row(1,0,list_WOF_head2_complete)

		int_WOF_i = 2

		for list_WOF_temp_sample in list_WOF_MRD_scores:
			
			obj_WOF_worksheet.write_row(int_WOF_i,0,list_WOF_temp_sample)

			int_WOF_i += 1


def main():

	WX_argparse = argparse.ArgumentParser(description='Output MRD scores to Excel file')
	WX_argparse.add_argument('--run_name', type=str, required=True, help='run name')
	WX_argparse.add_argument('--run_info_file', type=str, required=True, help='Path to run file')
	WX_argparse.add_argument('--MID_file', type=str, required=True, help='Path to MID sequence file')
	WX_argparse.add_argument('--gPrimer_file', type=str, required=True, help='Path to gPrimer sequence file')
	WX_argparse.add_argument('--pileup_folder', type=str, required=True, help='Path to folder, where pileup files are stored')
	WX_argparse.add_argument('--fastq_folder', type=str, required=True, help='Path to folder, where demultiplexed fastq files are stored')
	WX_argparse.add_argument('--aln_folder', type=str, required=True, help='Path to folder, where aligned bam files are stored')
	WX_argparse.add_argument('--output_file', type=str, required=True, help='Path to output Excel file')

	WX_args = WX_argparse.parse_args()

	str_CMSR_pileup_folder = WX_args.pileup_folder + ('/' if WX_args.pileup_folder[-1] != '/' else '')
	str_CMSR_fastq_folder = WX_args.fastq_folder + ('/' if WX_args.fastq_folder[-1] != '/' else '')
	str_CMSR_aln_folder = WX_args.aln_folder + ('/' if WX_args.aln_folder[-1] != '/' else '')

	float_CMSR_LOD_stdev_test = 3
	float_CMSR_LOD_stdev_cutoff_BE = 2.5
	int_CMSR_NOP_target_range = 25
	int_CMSR_R1_Nsupp_cutoff = 75

	int_CMSR_length_illumina_read = 250
	int_CMSR_length_CS = 20
	int_CMSR_length_RF = 16
	int_CMSR_buffer_to_read_end = 15


	list_CMSR_samples = []

	with open(WX_args.run_info_file,'r') as obj_CMSR_run_file:
		it_CMSR_run_samples_csv = csv.DictReader(obj_CMSR_run_file,delimiter='\t')

		for dict_CMSR_temp_run_sample in it_CMSR_run_samples_csv:

			str_CMSR_temp_base_filename = f'{WX_args.run_name}_{dict_CMSR_temp_run_sample["Number"].zfill(2)}_{dict_CMSR_temp_run_sample["Primer"]}_{dict_CMSR_temp_run_sample["Sample"]}'

			str_CMSR_target_region_chr, list_CMSR_target_region_bounds = get_analysis_range(dict_CMSR_temp_run_sample['Primer'],dict_CMSR_temp_run_sample['MID_F'],dict_CMSR_temp_run_sample['MID_R'], \
				WX_args.MID_file,WX_args.gPrimer_file,int_CMSR_length_illumina_read,int_CMSR_length_CS,int_CMSR_length_RF,int_CMSR_buffer_to_read_end)


			list_CMSR_temp_target_base_change = dict_CMSR_temp_run_sample['Variant'].split('>')
			if len(list_CMSR_temp_target_base_change) != 2:
				raise ValueError('Variant not recognized')


			list_CMSR_samples.append({'number':int(dict_CMSR_temp_run_sample["Number"]),'sample':dict_CMSR_temp_run_sample['Sample'],'primer':dict_CMSR_temp_run_sample['Primer'], \
				'base_filename':str_CMSR_pileup_folder+str_CMSR_temp_base_filename, 'fastq_file':str_CMSR_fastq_folder + str_CMSR_temp_base_filename + '_F.fq', 'bam_file_RF':str_CMSR_aln_folder+str_CMSR_temp_base_filename + '_RF.bam', \
				'target':{'chr':dict_CMSR_temp_run_sample['Chr'], 'position':int(dict_CMSR_temp_run_sample['Coordinate']),'base_change':list_CMSR_temp_target_base_change}, \
				'target_region':{'chr':str_CMSR_target_region_chr,'bounds':list_CMSR_target_region_bounds}})

	

	write_output_file(WX_args.output_file,get_mrd_scores(list_CMSR_samples,float_CMSR_LOD_stdev_test,float_CMSR_LOD_stdev_cutoff_BE,int_CMSR_NOP_target_range,int_CMSR_R1_Nsupp_cutoff))

if __name__ == '__main__':
	main()