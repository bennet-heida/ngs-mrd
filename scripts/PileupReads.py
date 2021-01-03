#Pileup NGS-reads to get base counts at each position of the target region
#run process_file() function with following arguments
#str_PF_path_input_file: path to input aligned bam file, bam file have RF-tags at the beginning of read name and must be sorted by read names
#str_PF_path_output_file: path to output file
#str_PF_path_reference_fasta: path to reference fasta file used for alignment
#str_PF_gReg_chr: target-region chromosome
#int_PF_gReg_start: target-region start position
#int_PF_gReg_len: length of target-region
#bool_PF_zero_based: target position base is 0 or 1

import pysam
import itertools
from BaseCount import pileup_reads
import csv
import argparse


def process_file(str_PF_path_input_file,str_PF_path_output_file,str_PF_path_reference_fasta,str_PF_gReg_chr,int_PF_gReg_start,int_PF_gReg_len,bool_PF_zero_based):
	#function calculates with 0-base coordinates
	if bool_PF_zero_based is False:
		int_PF_gReg_start -= 1

	str_PF_gReg_reference_fasta_output = pysam.faidx(str_PF_path_reference_fasta, \
		f'{str_PF_gReg_chr}:{int_PF_gReg_start+1}-{int_PF_gReg_start+1+int_PF_gReg_len-1}')
	str_PF_gReg_reference_seq = "".join(str_PF_gReg_reference_fasta_output.split('\n')[1:])

	print(f'Pileup file: {str_PF_path_input_file}')

	with pysam.AlignmentFile(str_PF_path_input_file) as obj_PF_input_alignment_file:
		it_PF_input_iterator = obj_PF_input_alignment_file.fetch(until_eof=True)
		list_PF_base_counts = pileup_reads(obj_PF_input_alignment_file,it_PF_input_iterator,str_PF_gReg_chr,int_PF_gReg_start,str_PF_gReg_reference_seq)


	if bool_PF_zero_based is False:
		int_PF_gReg_start += 1

	with open(str_PF_path_output_file, 'w') as obj_PF_output_file:
		obj_PF_output_csv = csv.writer(obj_PF_output_file,delimiter='\t')

		for int_PF_i in range(int_PF_gReg_len):

			obj_PF_output_csv.writerow([str_PF_gReg_chr,int_PF_gReg_start+int_PF_i,str_PF_gReg_reference_seq[int_PF_i], \
				*list_PF_base_counts[int_PF_i][0]," ".join(f'{str_PF_temp_indel}:{int_PF_temp_indel_count}' \
					for str_PF_temp_indel, int_PF_temp_indel_count in list_PF_base_counts[int_PF_i][1])," ".join(f'{str_PF_temp_indel}:{int_PF_temp_indel_count}' \
					for str_PF_temp_indel, int_PF_temp_indel_count in list_PF_base_counts[int_PF_i][2])])




def main():

	PR_argparse = argparse.ArgumentParser(description='Pileup NGS-Reads')
	PR_argparse.add_argument('--input_file', type=str, required=True, help='Path to aligned SAM/BAM-file')
	PR_argparse.add_argument('--output_file', type=str, required=True, help='Path to output file')
	PR_argparse.add_argument('--ref_fasta_file', type=str, required=True, help='Path to reference fasta file used for alignment')
	PR_argparse.add_argument('--gReg_chr', type=str, required=True, help='Target chromosome')
	PR_argparse.add_argument('--gReg_target_start', type=int, required=True, help='First base of target region')
	PR_argparse.add_argument('--gReg_target_length', type=int, required=True, help='First base of target region')
	PR_argparse.add_argument('--gReg_target_base', type=int, choices=[0,1], required=True, help='Target start is 0-based or 1-based')

	PR_args = PR_argparse.parse_args()

	bool_PR_target_zero_based = False if PR_args.gReg_target_base == 1 else True

	process_file(PR_args.input_file,PR_args.output_file,PR_args.ref_fasta_file,PR_args.gReg_chr,PR_args.gReg_target_start,PR_args.gReg_target_length,bool_PR_target_zero_based)


if __name__ == '__main__':
	main()