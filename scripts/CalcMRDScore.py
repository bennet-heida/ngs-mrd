#Calculate MRD scores from pileup files
import csv
import statistics
from NGSMRDFuncs import division_zero_tolerant
import argparse

def calc_be_rate(list_CBR_target_range_VAF,int_CBR_target_position,float_CBR_LOD_stdev_cutoff_BE):
	#exclude target position
	list_CBR_BE_range_VAF = [float_CBR_BE_VAF for int_CBR_BE_position,float_CBR_BE_VAF in list_CBR_target_range_VAF if int_CBR_BE_position != int_CBR_target_position]

	float_CBR_BE_range_mean = statistics.mean(list_CBR_BE_range_VAF)
	float_CBR_BE_range_stdev = statistics.stdev(list_CBR_BE_range_VAF)

	#exclude positions that deviate from mean by more than 2.5 standard deviations
	list_CBR_BE_range_corrected_VAF = [float_CBR_BE_corrected_VAF for float_CBR_BE_corrected_VAF in list_CBR_BE_range_VAF \
		if float_CBR_BE_range_mean - float_CBR_LOD_stdev_cutoff_BE*float_CBR_BE_range_stdev <= float_CBR_BE_corrected_VAF <= float_CBR_BE_range_mean + float_CBR_LOD_stdev_cutoff_BE*float_CBR_BE_range_stdev]

	float_CBR_BE_range_corrected_mean = statistics.mean(list_CBR_BE_range_corrected_VAF)
	float_CBR_BE_range_corrected_stdev = statistics.stdev(list_CBR_BE_range_corrected_VAF)

	return [[len(list_CBR_BE_range_corrected_VAF),len(list_CBR_BE_range_VAF)],[float_CBR_BE_range_corrected_mean,float_CBR_BE_range_corrected_stdev]]


def calc_lvafs(list_CLV_target_range_base_counts):
	dict_CLV_bases = {'A':0,'C':1,'G':2,'T':3,'N':4}

	list_CLV_target_range_LVAF = [(int_CLV_temp_range_position,division_zero_tolerant(max([list_CLV_temp_base_counts[1][int_CLV_i] for int_CLV_i in range(4) if int_CLV_i != dict_CLV_bases[list_CLV_temp_base_counts[0]]]), \
		sum(list_CLV_temp_base_counts[1][:4]))) for int_CLV_temp_range_position, list_CLV_temp_base_counts in list_CLV_target_range_base_counts]

	return list_CLV_target_range_LVAF


def calc_tvafs(list_CTV_target_range_base_counts,*list_CTV_target_base_changes):
	dict_CTV_bases = {'A':0,'C':1,'G':2,'T':3,'N':4}

	dict_CTV_base_changes = {}
	for list_CTV_temp_base_change in list_CTV_target_base_changes:
		if not list_CTV_temp_base_change[0] in dict_CTV_base_changes:
			dict_CTV_base_changes[list_CTV_temp_base_change[0]] = []
		dict_CTV_base_changes[list_CTV_temp_base_change[0]].append(dict_CTV_bases[list_CTV_temp_base_change[1]])

	list_CTV_target_range_TVAF = [(int_CTV_temp_range_position,division_zero_tolerant(max([list_CTV_temp_base_counts[1][int_CTV_i] for int_CTV_i in dict_CTV_base_changes[list_CTV_temp_base_counts[0]]]), \
		sum(list_CTV_temp_base_counts[1][:4]))) for int_CTV_temp_range_position, list_CTV_temp_base_counts in list_CTV_target_range_base_counts \
		if list_CTV_temp_base_counts[0] in dict_CTV_base_changes]

	return list_CTV_target_range_TVAF



def read_pileup_file(str_RPF_path_input_file,str_RPF_target_chr,int_RPF_target_range_start,int_RPF_target_range_len,int_RPF_target_position):
	with open(str_RPF_path_input_file,'r') as obj_RPF_input_file:
		it_RPF_input_csv = csv.reader(obj_RPF_input_file,delimiter='\t')

		dict_RPF_input_file = {}

		for list_RPF_temp_input_line in it_RPF_input_csv:
			if list_RPF_temp_input_line[0] == str_RPF_target_chr:
				#format indel counts
				list_RPF_temp_input_indels = [{str_RPF_temp_item_bchange:int(str_RPF_temp_item_count) for str_RPF_temp_item_bchange,str_RPF_temp_item_count \
					in [str_RPF_temp_indel_item.split(':') for str_RPF_temp_indel_item in list_RPF_temp_indel_type.split(' ') if str_RPF_temp_indel_item != '']} \
					for list_RPF_temp_indel_type in list_RPF_temp_input_line[8:10]]
				#create dictionary entry
				dict_RPF_input_file[int(list_RPF_temp_input_line[1])] = [list_RPF_temp_input_line[2],[int(str_RPF_temp_input_base) for str_RPF_temp_input_base \
					in list_RPF_temp_input_line[3:8]],list_RPF_temp_input_indels]

		#pull data for target region from dictionary, this checks if all positions have been read from file aswell
		list_RPF_target_range_base_counts = [(int_RPF_target_range_start+int_RPF_i,dict_RPF_input_file[int_RPF_target_range_start+int_RPF_i]) \
			for int_RPF_i in range(int_RPF_target_range_len) if int_RPF_target_range_start + int_RPF_i != int_RPF_target_position]

		list_RPF_target_position_base_counts = dict_RPF_input_file[int_RPF_target_position]

		return (list_RPF_target_position_base_counts,list_RPF_target_range_base_counts)


def formal_testing(float_FT_target_VAF,list_FT_target_range_VAF,int_FT_target_position,float_FT_LOD_stdev_test,float_FT_LOD_stdev_cutoff_BE,int_FT_NOP_target_range):
	list_FT_target_range_be_rate = calc_be_rate(list_FT_target_range_VAF,int_FT_target_position,float_FT_LOD_stdev_cutoff_BE)

	float_FT_target_LOD = list_FT_target_range_be_rate[1][0] + float_FT_LOD_stdev_test*list_FT_target_range_be_rate[1][1]
	bool_FT_target_T1 = True if float_FT_target_VAF > float_FT_target_LOD else False

	float_FT_NOP_LOD = float_FT_target_VAF - float_FT_LOD_stdev_test*list_FT_target_range_be_rate[1][1]
	list_FT_NOP_positions = [int_FT_temp_position for int_FT_temp_position,float_FT_temp_VAF in list_FT_target_range_VAF \
		if float_FT_temp_VAF >= float_FT_NOP_LOD and int_FT_temp_position != int_FT_target_position]
		
	#only account for NOPs in vicinity of 25 bases around target
	list_FT_NOP_positions = [int_FT_temp_position for int_FT_temp_position in list_FT_NOP_positions if abs(int_FT_temp_position - int_FT_target_position) <= int_FT_NOP_target_range]

	int_FT_NOP = len(list_FT_NOP_positions)
	bool_FT_target_T2 = True if int_FT_NOP == 0 else False

	bool_FT_target_final_T = all((bool_FT_target_T1,bool_FT_target_T2))

	return {'be':list_FT_target_range_be_rate,'lod':[float_FT_target_LOD,float_FT_NOP_LOD],'nop':int_FT_NOP,'t':[bool_FT_target_T1,bool_FT_target_T2,bool_FT_target_final_T]}

def calc_snp_score(list_CSS_target_position_base_counts,list_CSS_target_range_base_counts,int_CSS_target_position,list_CSS_target_base_change,float_CSS_LOD_stdev_test,float_CSS_LOD_stdev_cutoff_BE,int_CSS_NOP_target_range):
	if not len(list_CSS_target_base_change[0]) == len(list_CSS_target_base_change[1]) == 1:
		raise ValueError('SNP base change not recognized')

	dict_CSS_bases = {'A':0,'C':1,'G':2,'T':3,'N':4}

	int_CSS_target_Nsupp = list_CSS_target_position_base_counts[1][dict_CSS_bases[list_CSS_target_base_change[1]]]
	int_CSS_target_Rdepth = sum(list_CSS_target_position_base_counts[1][:4])

	float_CSS_target_VAF = division_zero_tolerant(int_CSS_target_Nsupp,int_CSS_target_Rdepth)

	list_CSS_target_range_VAF = [calc_lvafs(list_CSS_target_range_base_counts),calc_tvafs(list_CSS_target_range_base_counts,list_CSS_target_base_change)]

	list_CSS_target_formal_test = [formal_testing(float_CSS_target_VAF,list_CSS_temp_target_range_VAF,int_CSS_target_position,float_CSS_LOD_stdev_test,float_CSS_LOD_stdev_cutoff_BE,int_CSS_NOP_target_range) for list_CSS_temp_target_range_VAF in list_CSS_target_range_VAF]

	return [int_CSS_target_Nsupp,int_CSS_target_Rdepth,float_CSS_target_VAF,list_CSS_target_formal_test]

def calc_indel_score(list_CIS_target_position_base_counts,list_CIS_target_range_base_counts,int_CIS_target_position,list_CIS_target_base_change,float_CIS_LOD_stdev_test,float_CIS_LOD_stdev_cutoff_BE,int_CIS_NOP_target_range):
	if len(list_CIS_target_base_change[0]) == len(list_CIS_target_base_change[1]):
		raise ValueError('Indel not recognized')

	#get indel sequence from base change
	if len(list_CIS_target_base_change[0]) < len(list_CIS_target_base_change[1]):
		int_CIS_indel_type = 0
		str_CIS_indel_seq = list_CIS_target_base_change[1][1:]
	else:
		int_CIS_indel_type = 1
		str_CIS_indel_seq = list_CIS_target_base_change[0][1:]

	#get Nsupp from pileup
	if str_CIS_indel_seq in list_CIS_target_position_base_counts[2][int_CIS_indel_type]:
		int_CIS_target_Nsupp = list_CIS_target_position_base_counts[2][int_CIS_indel_type][str_CIS_indel_seq]
	else:
		int_CIS_target_Nsupp = 0

	int_CIS_target_Rdepth = sum(list_CIS_target_position_base_counts[1][:4])

	float_CIS_target_VAF = division_zero_tolerant(int_CIS_target_Nsupp,int_CIS_target_Rdepth)

	list_CIS_target_range_LVAF = calc_lvafs(list_CIS_target_range_base_counts)

	dict_CIS_target_formal_test = formal_testing(float_CIS_target_VAF,list_CIS_target_range_LVAF,int_CIS_target_position,float_CIS_LOD_stdev_test,float_CIS_LOD_stdev_cutoff_BE,int_CIS_NOP_target_range)

	return [int_CIS_target_Nsupp,int_CIS_target_Rdepth,float_CIS_target_VAF,[dict_CIS_target_formal_test]]


def calc_mrd_score(str_CMS_path_input_file,str_CMS_target_chr,int_CMS_target_range_start,int_CMS_target_range_len,int_CMS_target_position,list_CMS_target_base_change,float_CMS_LOD_stdev_test,float_CMS_LOD_stdev_cutoff_BE,int_CMS_NOP_target_range):
	#read pileup file
	list_CMS_target_position_base_counts, list_CMS_target_range_base_counts = read_pileup_file(str_CMS_path_input_file,str_CMS_target_chr, \
		int_CMS_target_range_start,int_CMS_target_range_len,int_CMS_target_position)

	#check target base change type
	if len(list_CMS_target_base_change[0]) == len(list_CMS_target_base_change[1]) == 1:
		list_CMS_output = calc_snp_score(list_CMS_target_position_base_counts,list_CMS_target_range_base_counts,int_CMS_target_position,list_CMS_target_base_change,float_CMS_LOD_stdev_test,float_CMS_LOD_stdev_cutoff_BE,int_CMS_NOP_target_range)
		
	else:
		list_CMS_output = calc_indel_score(list_CMS_target_position_base_counts,list_CMS_target_range_base_counts,int_CMS_target_position,list_CMS_target_base_change,float_CMS_LOD_stdev_test,float_CMS_LOD_stdev_cutoff_BE,int_CMS_NOP_target_range)
	
	return list_CMS_output

def main():

	CMS_argparse = argparse.ArgumentParser(description='Calculate MRD scores')
	CMS_argparse.add_argument('--pileup_file', type=str, required=True, help='Input pileup file')
	CMS_argparse.add_argument('--target_chr', type=str, required=True, help='target chromosome')
	CMS_argparse.add_argument('--target_position', type=int, required=True, help='target position')
	CMS_argparse.add_argument('--target_base_change', type=str, required=True, help='target mutation')
	CMS_argparse.add_argument('--gReg_start', type=int, required=True, help='start of genomic region to be considered for calculations')
	CMS_argparse.add_argument('--gReg_len', type=int, required=True, help='length of genomic region to be considered for calculations')

	CMR_args = CMS_argparse.parse_args()

	float_CS_LOD_stdev_test = 3
	float_CS_LOD_stdev_cutoff_BE = 2.5
	int_CS_NOP_target_range = 25

	list_CS_target_bchange = CMR_args.target_base_change.split('>')
	if not len(list_CS_target_bchange) == 2:
		raise ValueError('Target base change not recognized')

	list_CS_MRDscore_result = calc_mrd_score(CMR_args.pileup_file,CMR_args.target_chr,CMR_args.gReg_start,CMR_args.gReg_len,CMR_args.target_position,list_CS_target_bchange,float_CS_LOD_stdev_test,float_CS_LOD_stdev_cutoff_BE,int_CS_NOP_target_range)

	print('Target region:', f'{CMR_args.target_chr}:{CMR_args.gReg_start}-{CMR_args.gReg_start + CMR_args.gReg_len -1}')
	print('Target mutation:', f'{CMR_args.target_position}:{CMR_args.target_base_change}')
	print(list_CS_MRDscore_result)

if __name__ == '__main__':
	main()