#Pileup read-families, takes sorted SAM/BAM-file as input, with RFs as QNAME-headers
#run process_file() function with following arguments:
#str_PF_path_input_file: path to input aligned bam file, bam file have RF-tags at the beginning of read name and must be sorted by read names
#str_PF_path_output_file: path to output file
#str_PF_path_reference_fasta: path to reference fasta file used for alignment
#str_PF_gReg_chr: target-region chromosome
#int_PF_gReg_start: target-region start position
#int_PF_gReg_len: length of target-region
#bool_PF_zero_based: target position base is 0 or 1
#int_PF_min_RF_size: min. number of reads to form a read family
#float_PF_min_consensus_fraction: min. fraction of reads to form consensus (must be >0.5)

import pysam
import itertools
from collections import defaultdict
from BaseCount import pileup_reads
import csv
import argparse


def process_file(str_PF_path_input_file,str_PF_path_output_file,str_PF_path_reference_fasta,str_PF_gReg_chr,int_PF_gReg_start,int_PF_gReg_len,bool_PF_zero_based,int_PF_min_RF_size,float_PF_min_consensus_fraction):
	#function calculates with 0-base coordinates
	if bool_PF_zero_based is False:
		int_PF_gReg_start -= 1

	str_PF_gReg_reference_fasta_output = pysam.faidx(str_PF_path_reference_fasta, \
		f'{str_PF_gReg_chr}:{int_PF_gReg_start+1}-{int_PF_gReg_start+1+int_PF_gReg_len-1}')
	str_PF_gReg_reference_seq = "".join(str_PF_gReg_reference_fasta_output.split('\n')[1:])

	print(f'Pileup file: {str_PF_path_input_file}')


	with pysam.AlignmentFile(str_PF_path_input_file) as obj_PF_input_alignment_file:
		it_PF_input_iterator = obj_PF_input_alignment_file.fetch(until_eof=True)
		list_PF_RF_base_count = pileup_RF(obj_PF_input_alignment_file,it_PF_input_iterator,str_PF_gReg_reference_seq,str_PF_gReg_chr,int_PF_gReg_start,int_PF_gReg_len,int_PF_min_RF_size,float_PF_min_consensus_fraction)


	if bool_PF_zero_based is False:
		int_PF_gReg_start += 1

	with open(str_PF_path_output_file, 'w') as obj_PF_output_file:
		obj_PF_output_csv = csv.writer(obj_PF_output_file,delimiter='\t')
		for int_PF_i in range(int_PF_gReg_len):
			obj_PF_output_csv.writerow([str_PF_gReg_chr,int_PF_gReg_start+int_PF_i,str_PF_gReg_reference_seq[int_PF_i],*list_PF_RF_base_count[int_PF_i][0], \
				" ".join([f'{str_PF_temp_indel}:{int_PF_temp_indel_count}' for str_PF_temp_indel, int_PF_temp_indel_count in list_PF_RF_base_count[int_PF_i][1]]), \
				" ".join([f'{str_PF_temp_indel}:{int_PF_temp_indel_count}' for str_PF_temp_indel, int_PF_temp_indel_count in list_PF_RF_base_count[int_PF_i][2]])])



def pileup_RF(obj_PRF_input_alignment_file,it_PRF_input_iterator,str_PRF_gReg_reference_seq,str_PRF_gReg_chr,int_PRF_gReg_start,int_PRF_gReg_len,int_PRF_min_RF_size,float_PRF_min_consensus_fraction):

	it_PRF_RF_iterator = filter(lambda obj_PRF_temp_RF:len(obj_PRF_temp_RF[0])>0, itertools.groupby(it_PRF_input_iterator, \
			key=lambda obj_PRF_temp_read:obj_PRF_temp_read.query_name.split('/',1)[0]))


	list_PRF_RF_pileup = [[[0,0,0,0,0],[],[]] for _ in range(int_PRF_gReg_len)]
	list_PRF_RF_temp_indel_pileup = [[defaultdict(int),defaultdict(int)] for _ in range(int_PRF_gReg_len)]

	for str_PRF_RF_tag, it_PRF_temp_RF_members in it_PRF_RF_iterator:

		list_PRF_temp_RF_members = list(it_PRF_temp_RF_members)

		if len(list_PRF_temp_RF_members) >= int_PRF_min_RF_size:
			list_PRF_temp_RF_consensus = collapse_RF(obj_PRF_input_alignment_file,list_PRF_temp_RF_members,str_PRF_gReg_reference_seq, \
				str_PRF_gReg_chr,int_PRF_gReg_start,int_PRF_gReg_len,int_PRF_min_RF_size,float_PRF_min_consensus_fraction)

			#add RF to pileup count
			for int_PRF_gReg_position in range(int_PRF_gReg_len):
				if not list_PRF_temp_RF_consensus[int_PRF_gReg_position][0] is None:
					list_PRF_RF_pileup[int_PRF_gReg_position][0][list_PRF_temp_RF_consensus[int_PRF_gReg_position][0]] += 1

				#inserts
				for str_PRF_temp_insert in list_PRF_temp_RF_consensus[int_PRF_gReg_position][1]:
					list_PRF_RF_temp_indel_pileup[int_PRF_gReg_position][0][str_PRF_temp_insert] += 1
				#delets
				for str_PRF_temp_delet in list_PRF_temp_RF_consensus[int_PRF_gReg_position][2]:
					list_PRF_RF_temp_indel_pileup[int_PRF_gReg_position][1][str_PRF_temp_delet] += 1

	for int_PRF_gReg_position in range(int_PRF_gReg_len):
		list_PRF_RF_pileup[int_PRF_gReg_position][1] = sorted(list_PRF_RF_temp_indel_pileup[int_PRF_gReg_position][0].items(),key=lambda int_PRF_tuple_item:int_PRF_tuple_item[1])
		list_PRF_RF_pileup[int_PRF_gReg_position][2] = sorted(list_PRF_RF_temp_indel_pileup[int_PRF_gReg_position][1].items(),key=lambda int_PRF_tuple_item:int_PRF_tuple_item[1])

	return list_PRF_RF_pileup



def collapse_RF(obj_CRF_input_alignment_file,list_CRF_RF_members,str_CRF_gReg_reference_seq,str_CRF_gReg_chr,int_CRF_gReg_start,int_CRF_gReg_len,int_CRF_min_RF_size,float_CRF_min_consensus_fraction):
	#number of members in RF
	int_CRF_RF_size = len(list_CRF_RF_members)
	#pileup RF
	list_CRF_RF_base_counts = pileup_reads(obj_CRF_input_alignment_file,list_CRF_RF_members,str_CRF_gReg_chr,int_CRF_gReg_start,str_CRF_gReg_reference_seq)

	#collapse RF to consensus sequence
	list_CRF_RF_consensus = [[None,[],[]] for _ in range(int_CRF_gReg_len)]
	for int_CRF_gReg_position in range(int_CRF_gReg_len):
		int_CRF_temp_depth = sum(list_CRF_RF_base_counts[int_CRF_gReg_position][0])
		#check read depth at position
		if int_CRF_temp_depth >= int_CRF_min_RF_size:
			int_CRF_temp_consensus_base = list_CRF_RF_base_counts[int_CRF_gReg_position][0].index(max(list_CRF_RF_base_counts[int_CRF_gReg_position][0]))

			#check if consensus base passes filter
			if list_CRF_RF_base_counts[int_CRF_gReg_position][0][int_CRF_temp_consensus_base] / int_CRF_temp_depth >= float_CRF_min_consensus_fraction:
				list_CRF_RF_consensus[int_CRF_gReg_position][0] = int_CRF_temp_consensus_base
			#else add 'N'
			else:
				list_CRF_RF_consensus[int_CRF_gReg_position][0] = 4


		#check inserts
		for tuple_CRF_temp_insert in list_CRF_RF_base_counts[int_CRF_gReg_position][1]:
			if tuple_CRF_temp_insert[1] / int_CRF_RF_size >= float_CRF_min_consensus_fraction:
				list_CRF_RF_consensus[int_CRF_gReg_position][1].append(tuple_CRF_temp_insert[0])

		#check delets
		for tuple_CRF_temp_delet in list_CRF_RF_base_counts[int_CRF_gReg_position][2]:
			if tuple_CRF_temp_delet[1] / int_CRF_RF_size >= float_CRF_min_consensus_fraction:
				list_CRF_RF_consensus[int_CRF_gReg_position][2].append(tuple_CRF_temp_delet[0])

	return list_CRF_RF_consensus


def main():

	PRF_argparse = argparse.ArgumentParser(description='Pileup NGS-Reads')
	PRF_argparse.add_argument('--input_file', type=str, required=True, help='Path to aligned SAM/BAM-file')
	PRF_argparse.add_argument('--output_file', type=str, required=True, help='Path to output file')
	PRF_argparse.add_argument('--ref_fasta_file', type=str, required=True, help='Path to reference fasta file used for alignment')
	PRF_argparse.add_argument('--gReg_chr', type=str, required=True, help='Target chromosome')
	PRF_argparse.add_argument('--gReg_target_start', type=int, required=True, help='First base of target region')
	PRF_argparse.add_argument('--gReg_target_length', type=int, required=True, help='First base of target region')
	PRF_argparse.add_argument('--gReg_target_base', type=int, choices=[0,1], required=True, help='Target start is 0-based or 1-based')
	PRF_argparse.add_argument('--min_member_count', type=int, required=True, help='min number of reads to form read family')
	PRF_argparse.add_argument('--min_consensus_percent', type=int, required=True, help='min percentage of bases to accept consensus')

	PRF_args = PRF_argparse.parse_args()

	bool_PRF_target_zero_based = False if PRF_args.gReg_target_base == 1 else True

	float_PRF_min_consensus_fraction = PRF_args.min_consensus_percent / 100

	process_file(PRF_args.input_file,PRF_args.output_file,PRF_args.ref_fasta_file,PRF_args.gReg_chr,PRF_args.gReg_target_start,PRF_args.gReg_target_length, \
		bool_PRF_target_zero_based,PRF_args.min_member_count,float_PRF_min_consensus_fraction)

if __name__ == '__main__':
	main()