#Plots target-VAF and background error for given pileup file, if target-mutation is SNP
#NOP criterion is shown aswell

import csv
from NGSMRDFuncs import division_zero_tolerant
import matplotlib.pyplot as plt
import statistics
import argparse

def read_pileup_file(str_RPF_path_input_file,str_RPF_target_chr):
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


	return dict_RPF_input_file

def calc_LVAF(list_CLV_base_count):
	dict_CLV_bases = {'A':0,'C':1,'G':2,'T':3}

	list_CLV_base_changes = [list_CLV_base_count[1][int_CLV_i] for int_CLV_i in range(4) if int_CLV_i != dict_CLV_bases[list_CLV_base_count[0]]]

	float_CLV_LVAF = division_zero_tolerant(max(list_CLV_base_changes), sum(list_CLV_base_count[1][:4]))

	return float_CLV_LVAF*100

def get_TVAF_base_changes(list_GTBC_base_changes):
	dict_GTBC_bases = {'A':0,'C':1,'G':2,'T':3}

	dict_GTBC_base_changes = {}

	for list_GTBC_temp_base_change in list_GTBC_base_changes:
		if not list_GTBC_temp_base_change[0] in dict_GTBC_base_changes:
			dict_GTBC_base_changes[list_GTBC_temp_base_change[0]] = []

		dict_GTBC_base_changes[list_GTBC_temp_base_change[0]].append(dict_GTBC_bases[list_GTBC_temp_base_change[1]])

	return dict_GTBC_base_changes

def calc_TVAF(list_CTV_base_count,dict_CTV_base_changes):
	if not list_CTV_base_count[0] in dict_CTV_base_changes:
		raise ValueError('No TVAF base at position')

	float_CTV_TVAF = division_zero_tolerant(max([list_CTV_base_count[1][int_CTV_temp_base] for int_CTV_temp_base in dict_CTV_base_changes[list_CTV_base_count[0]]]), \
		sum(list_CTV_base_count[1][:4]))

	return float_CTV_TVAF*100


def target_VAF_SNP(list_TVS_base_count,list_TVS_target_base_change):
	dict_TVS_bases = {'A':0,'C':1,'G':2,'T':3}

	if not list_TVS_base_count[0] == list_TVS_target_base_change[0]:
		raise ValueError('Unexpected reference base at target position')

	int_TVS_target_nsupp = list_TVS_base_count[1][dict_TVS_bases[list_TVS_target_base_change[1]]]

	float_TVS_target_VAF = division_zero_tolerant(int_TVS_target_nsupp,sum(list_TVS_base_count[1][:4]))

	return float_TVS_target_VAF*100

def target_VAF_indel(list_TVI_base_count,list_TVI_target_base_change):
	#insertion
	if len(list_TVI_target_base_change[1]) > len(list_TVI_target_base_change[0]):
		int_TVI_indel_type = 0
		str_TVI_indel_seq = list_TVI_target_base_change[1][1:]
	elif len(list_TVI_target_base_change[1]) < len(list_TVI_target_base_change[0]):
		int_TVI_indel_type = 1
		str_TVI_indel_seq = list_TVI_target_base_change[0][1:]
	else:
		raise ValueError('Target indel not recognized')

	if str_TVI_indel_seq in list_TVI_base_count[2][int_TVI_indel_type]:
		int_TVI_target_nsupp = list_TVI_base_count[2][int_TVI_indel_type][str_TVI_indel_seq]
	else:
		int_TVI_target_nsupp = 0

	float_TVI_target_VAF = division_zero_tolerant(int_TVI_target_nsupp,sum(list_TVI_base_count[1][:4]))

	return float_TVI_target_VAF*100


def get_LVAF_SNP(dict_GLVS_input_file,int_GLVS_target_range_start,int_GLVS_target_range_len,int_GLVS_target_position,list_GLVS_target_base_change):
	
	float_GLVS_target_VAF = target_VAF_SNP(dict_GLVS_input_file[int_GLVS_target_position],list_GLVS_target_base_change)

	list_GLVS_VAF = [(int_GLVS_i,calc_LVAF(dict_GLVS_input_file[int_GLVS_i])) for int_GLVS_i in range(int_GLVS_target_range_start,int_GLVS_target_position)]
	list_GLVS_VAF += [(int_GLVS_target_position,float_GLVS_target_VAF)]
	list_GLVS_VAF += [(int_GLVS_i,calc_LVAF(dict_GLVS_input_file[int_GLVS_i])) for int_GLVS_i in range(int_GLVS_target_position+1,int_GLVS_target_range_start+int_GLVS_target_range_len)]

	return float_GLVS_target_VAF,list_GLVS_VAF




def get_TVAF_SNP(dict_GTVS_input_file,int_GTVS_target_range_start,int_GTVS_target_range_len,int_GTVS_target_position,list_GTVS_target_base_change):
	dict_GTVS_base_changes = get_TVAF_base_changes([list_GTVS_target_base_change])

	float_GTVS_target_VAF = target_VAF_SNP(dict_GTVS_input_file[int_GTVS_target_position],list_GTVS_target_base_change)

	list_GTVS_VAF = [(int_GTVS_i,calc_TVAF(dict_GTVS_input_file[int_GTVS_i],dict_GTVS_base_changes)) for int_GTVS_i in range(int_GTVS_target_range_start,int_GTVS_target_position) if dict_GTVS_input_file[int_GTVS_i][0] in dict_GTVS_base_changes]
	list_GTVS_VAF += [(int_GTVS_target_position,float_GTVS_target_VAF)]
	list_GTVS_VAF += [(int_GTVS_i,calc_TVAF(dict_GTVS_input_file[int_GTVS_i],dict_GTVS_base_changes)) for int_GTVS_i in range(int_GTVS_target_position+1,int_GTVS_target_range_start+int_GTVS_target_range_len) if dict_GTVS_input_file[int_GTVS_i][0] in dict_GTVS_base_changes]

	return float_GTVS_target_VAF,list_GTVS_VAF


def get_VAF_indel(dict_GVI_input_file,int_GVI_target_range_start,int_GVI_target_range_len,int_GVI_target_position,list_GVI_target_base_change):

	float_GVI_target_VAF = target_VAF_indel(dict_GVI_input_file[int_GVI_target_position],list_GVI_target_base_change)

	list_GVI_VAF = [(int_GVI_i,calc_LVAF(dict_GVI_input_file[int_GVI_i])) for int_GVI_i in range(int_GVI_target_range_start,int_GVI_target_position)]
	list_GVI_VAF += [(int_GVI_target_position,float_GVI_target_VAF)]
	list_GVI_VAF += [(int_GVI_i,calc_LVAF(dict_GVI_input_file[int_GVI_i])) for int_GVI_i in range(int_GVI_target_position+1,int_GVI_target_range_start+int_GVI_target_range_len)]

	return float_GVI_target_VAF,list_GVI_VAF

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


def plot_pileup_SNP_LVAF(str_PPSL_path_input_file,str_PPSL_path_output_file,str_PPSL_target_range_chr,int_PPSL_target_range_start,int_PPSL_target_range_len,str_PPSL_target_chr, \
		int_PPSL_target_position,list_PPSL_target_base_change,float_PPSL_LOD_stdev_test,float_PPSL_LOD_stdev_cutoff_BE):
	dict_PPSL_input_file = read_pileup_file(str_PPSL_path_input_file,str_PPSL_target_range_chr)

	float_PPSL_target_VAF, list_PPSL_VAF = get_LVAF_SNP(dict_PPSL_input_file,int_PPSL_target_range_start,int_PPSL_target_range_len, \
		int_PPSL_target_position,list_PPSL_target_base_change)

	list_PPSL_BE_rate = calc_be_rate(list_PPSL_VAF,int_PPSL_target_position,float_PPSL_LOD_stdev_cutoff_BE)
	float_PPSL_be_av, float_PPSL_be_stdev = list_PPSL_BE_rate[1]

	float_PPSL_test_LOD = float_PPSL_be_av + float_PPSL_LOD_stdev_test*float_PPSL_be_stdev
	float_PPSL_NOP_LOD = float_PPSL_target_VAF - float_PPSL_LOD_stdev_test*float_PPSL_be_stdev

	list_PPSL_x, list_PPSL_y = list(zip(*list_PPSL_VAF))

	plt.figure()
	plt.plot(list_PPSL_x,list_PPSL_y,linewidth=3,color='blue',zorder=3)

	plt.yscale('symlog',linthreshy=0.03,linscaley=0.5)
	plt.hlines(float_PPSL_test_LOD,int_PPSL_target_range_start,int_PPSL_target_range_start+int_PPSL_target_range_len,zorder=2,linewidth=1,color='k')
	if float_PPSL_target_VAF > float_PPSL_test_LOD:
		plt.hlines(float_PPSL_NOP_LOD,int_PPSL_target_range_start,int_PPSL_target_range_start+int_PPSL_target_range_len,zorder=2,linewidth=1,color='k',linestyle='--')

	plt.axvline(int_PPSL_target_position,0.03,0.97,linewidth=1,color='k',zorder=1)
	plt.ylabel('VAF%')
	plt.xlabel('Chr Position')

	plt.plot([list_PPSL_x[0],list_PPSL_x[0]],[0,0.3],'o',color='none')

	plt.savefig(str_PPSL_path_output_file, format='svg')
	plt.close()


def plot_pileup_SNP_TVAF(str_PPST_path_input_file,str_PPST_path_output_file,str_PPST_target_range_chr,int_PPST_target_range_start,int_PPST_target_range_len,str_PPST_target_chr, \
	int_PPST_target_position,list_PPST_target_base_change,float_PPST_LOD_stdev_test,float_PPST_LOD_stdev_cutoff_BE):
	dict_PPST_input_file = read_pileup_file(str_PPST_path_input_file,str_PPST_target_range_chr)

	float_PPST_target_VAF, list_PPST_VAF = get_TVAF_SNP(dict_PPST_input_file,int_PPST_target_range_start,int_PPST_target_range_len, \
		int_PPST_target_position,list_PPST_target_base_change)

	list_PPST_BE_rate = calc_be_rate(list_PPST_VAF,int_PPST_target_position,float_PPST_LOD_stdev_cutoff_BE)

	float_PPST_be_av, float_PPST_be_stdev = list_PPST_BE_rate[1]

	float_PPST_test_LOD = float_PPST_be_av + float_PPST_LOD_stdev_test*float_PPST_be_stdev
	float_PPST_NOP_LOD = float_PPST_target_VAF - float_PPST_LOD_stdev_test*float_PPST_be_stdev

	list_PPST_x, list_PPST_y = list(zip(*list_PPST_VAF))

	plt.figure()
	plt.plot(list_PPST_x,list_PPST_y,'o',linewidth=3,color='blue',zorder=3)

	plt.yscale('symlog',linthreshy=0.03,linscaley=0.5)
	plt.hlines(float_PPST_test_LOD,int_PPST_target_range_start,int_PPST_target_range_start+int_PPST_target_range_len,zorder=2,linewidth=1,color='k')
	if float_PPST_target_VAF > float_PPST_test_LOD:
		plt.hlines(float_PPST_NOP_LOD,int_PPST_target_range_start,int_PPST_target_range_start+int_PPST_target_range_len,zorder=2,linewidth=1,color='k',linestyle='--')

	plt.axvline(int_PPST_target_position,0.03,0.97,linewidth=1,color='k',zorder=1)
	plt.ylabel('VAF%')
	plt.xlabel('Chr Position')

	plt.plot([list_PPST_x[0],list_PPST_x[0]],[0,0.3],'o',color='none')

	plt.savefig(str_PPST_path_output_file, format='svg')
	plt.close()

def plot_pileup_indel(str_PPI_path_input_file,str_PPI_path_output_file,str_PPI_target_range_chr,int_PPI_target_range_start,int_PPI_target_range_len,str_PPI_target_chr, \
	int_PPI_target_position,list_PPI_target_base_change,float_PPI_LOD_stdev_test,float_PPI_LOD_stdev_cutoff_BE):
	dict_PPI_input_file = read_pileup_file(str_PPI_path_input_file,str_PPI_target_range_chr)

	float_PPI_target_VAF, list_PPI_VAF = get_VAF_indel(dict_PPI_input_file,int_PPI_target_range_start,int_PPI_target_range_len, \
		int_PPI_target_position,list_PPI_target_base_change)

	list_PPI_BE_rate = calc_be_rate(list_PPI_VAF,int_PPI_target_position,float_PPI_LOD_stdev_cutoff_BE)

	float_PPI_be_av, float_PPI_be_stdev = list_PPI_BE_rate[1]

	float_PPI_test_LOD = float_PPI_be_av + float_PPI_LOD_stdev_test*float_PPI_be_stdev
	float_PPI_NOP_LOD = float_PPI_target_VAF - float_PPI_LOD_stdev_test*float_PPI_be_stdev

	list_PPI_x, list_PPI_y = list(zip(*list_PPI_VAF))

	plt.figure()
	plt.plot(list_PPI_x,list_PPI_y,linewidth=3,color='blue',zorder=3)

	plt.yscale('symlog',linthreshy=0.03,linscaley=0.5)
	plt.hlines(float_PPI_test_LOD,int_PPI_target_range_start,int_PPI_target_range_start+int_PPI_target_range_len,zorder=2,linewidth=1,color='k')
	plt.axvline(int_PPI_target_position,0.03,0.97,linewidth=1,color='k',zorder=1)
	plt.ylabel('VAF%')
	plt.xlabel('Chr Position')

	plt.plot([list_PPI_x[0],list_PPI_x[0]],[0,0.3],'o',color='none')

	plt.savefig(str_PPI_path_output_file, format='svg')
	plt.close()

def main():

	PLOT_argparse = argparse.ArgumentParser(description='Plot VAF Distribution')
	PLOT_argparse.add_argument('--pileup_file', type=str, required=True, help='Input pileup file')
	PLOT_argparse.add_argument('--output_prefix', type=str, required=True, help='Output SVG-file')
	PLOT_argparse.add_argument('--target_chr', type=str, required=True, help='target chromosome')
	PLOT_argparse.add_argument('--target_position', type=int, required=True, help='target position')
	PLOT_argparse.add_argument('--target_base_change', type=str, required=True, help='target mutation')
	PLOT_argparse.add_argument('--gReg_start', type=int, required=True, help='start of genomic region to be considered for calculations')
	PLOT_argparse.add_argument('--gReg_len', type=int, required=True, help='length of genomic region to be considered for calculations')

	PLOT_args = PLOT_argparse.parse_args()

	list_PLOT_target_base_change = PLOT_args.target_base_change.split('>')

	if not len(list_PLOT_target_base_change) == 2:
		raise ValueError('Target base change not recognized')

	if len(list_PLOT_target_base_change[0]) == 0 or len(list_PLOT_target_base_change[1]) == 0:
		raise ValueError('Target base change not recognized')

	int_PLOT_indel_length = abs(len(list_PLOT_target_base_change[1])-len(list_PLOT_target_base_change[0]))

	float_PLOT_LOD_stdev_cutoff_BE = 2.5
	float_PLOT_LOD_stdev_test = 3

	if int_PLOT_indel_length == 0:
		plot_pileup_SNP_LVAF(PLOT_args.pileup_file,PLOT_args.output_prefix+'_LVAF.svg',PLOT_args.target_chr,PLOT_args.gReg_start,PLOT_args.gReg_len,PLOT_args.target_chr, \
			PLOT_args.target_position,list_PLOT_target_base_change,float_PLOT_LOD_stdev_test,float_PLOT_LOD_stdev_cutoff_BE)

		plot_pileup_SNP_TVAF(PLOT_args.pileup_file,PLOT_args.output_prefix+'_TVAF.svg',PLOT_args.target_chr,PLOT_args.gReg_start,PLOT_args.gReg_len,PLOT_args.target_chr, \
			PLOT_args.target_position,list_PLOT_target_base_change,float_PLOT_LOD_stdev_test,float_PLOT_LOD_stdev_cutoff_BE)

	else:
		plot_pileup_indel(PLOT_args.pileup_file,PLOT_args.output_prefix+'_indel.svg',PLOT_args.target_chr,PLOT_args.gReg_start,PLOT_args.gReg_len,PLOT_args.target_chr, \
			PLOT_args.target_position,list_PLOT_target_base_change,float_PLOT_LOD_stdev_test,float_PLOT_LOD_stdev_cutoff_BE)


if __name__ == '__main__':
	main()