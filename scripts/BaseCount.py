#Pileup NGS-reads to get base counts at each position of the target region
import pysam
from collections import defaultdict

def pileup_reads(obj_PR_input_alignment_file,it_PR_input_iterator,str_PR_gReg_chr,int_PR_gReg_start,str_PR_gReg_reference_seq):

	int_PR_gReg_len = len(str_PR_gReg_reference_seq)

	list_PR_base_counts = [[[0,0,0,0,0],[],[]] for _ in range(int_PR_gReg_len)]
	list_PR_indel_counts = [[defaultdict(int),defaultdict(int)] for _ in range(int_PR_gReg_len)]

	for obj_PR_read in it_PR_input_iterator:
		list_PR_temp_base_count = count_bases(obj_PR_input_alignment_file,obj_PR_read,str_PR_gReg_chr,int_PR_gReg_start,str_PR_gReg_reference_seq)


		if len(list_PR_temp_base_count) != int_PR_gReg_len:
			raise ValueError('Wrong length of base count list')

		for int_PR_gReg_pos in range(int_PR_gReg_len):
			#matching bases
			for int_PR_temp_base_match in list_PR_temp_base_count[int_PR_gReg_pos][0]:
				list_PR_base_counts[int_PR_gReg_pos][0][int_PR_temp_base_match] += 1

			#insertions on read
			for str_PR_temp_insert in list_PR_temp_base_count[int_PR_gReg_pos][1]:
				list_PR_indel_counts[int_PR_gReg_pos][0][str_PR_temp_insert] += 1

			#deletions on read
			for str_PR_temp_delet in list_PR_temp_base_count[int_PR_gReg_pos][2]:
				list_PR_indel_counts[int_PR_gReg_pos][1][str_PR_temp_delet] += 1

	for int_PR_i in range(int_PR_gReg_len):
		list_PR_base_counts[int_PR_i][1] = sorted(list_PR_indel_counts[int_PR_i][0].items(),key=lambda tuple_PR_item:tuple_PR_item[1],reverse=True)
		list_PR_base_counts[int_PR_i][2] = sorted(list_PR_indel_counts[int_PR_i][1].items(),key=lambda tuple_PR_item:tuple_PR_item[1],reverse=True)

	return list_PR_base_counts

def count_bases(obj_CB_input_alignment_file,obj_CB_read,str_CB_gReg_chr,int_CB_gReg_start,str_CB_gReg_reference_seq):

	int_CB_gReg_len = len(str_CB_gReg_reference_seq)

	dict_CB_base_list = {'A':0,'C':1,'G':2,'T':3,'N':4}
	list_CB_base_counts = [[[],[],[]] for _ in range(int_CB_gReg_len)]

	if not obj_CB_read.is_unmapped and not obj_CB_read.is_secondary and not obj_CB_read.is_supplementary:

		str_CB_read_chr = obj_CB_input_alignment_file.get_reference_name(obj_CB_read.reference_id)
		int_CB_read_start = obj_CB_read.reference_start - int_CB_gReg_start
		str_CB_read_sequence = obj_CB_read.query_sequence
		int_CB_read_len = len(str_CB_read_sequence)
		#read matches target region
		if str_CB_gReg_chr == str_CB_read_chr and int_CB_read_start <= int_CB_gReg_len and int_CB_read_start + int_CB_read_len >= 0:
			if obj_CB_read.reference_start >= obj_CB_read.reference_end:
				raise ValueError('Expected SAM/BAM-file with reads converted to forward direction')

			int_CB_read_pointer = 0
			int_CB_gReg_pointer = int_CB_read_start

			for tuple_CB_cigar_tuple in obj_CB_read.cigartuples:
				#match region
				if tuple_CB_cigar_tuple[0] == 0:
					for _ in range(tuple_CB_cigar_tuple[1]):
						#base matched to target region
						if 0 <= int_CB_gReg_pointer < int_CB_gReg_len:
							list_CB_base_counts[int_CB_gReg_pointer][0].append(dict_CB_base_list[str_CB_read_sequence[int_CB_read_pointer]])
						int_CB_read_pointer += 1
						int_CB_gReg_pointer += 1

				#insertion on read
				elif tuple_CB_cigar_tuple[0] == 1:
					#match to target region
					if 0 < int_CB_gReg_pointer < int_CB_gReg_len:
						list_CB_base_counts[int_CB_gReg_pointer-1][1].append(str_CB_read_sequence[int_CB_read_pointer:int_CB_read_pointer+tuple_CB_cigar_tuple[1]])
					int_CB_read_pointer += tuple_CB_cigar_tuple[1]

				#deletion on read
				elif tuple_CB_cigar_tuple[0] == 2:
					#match to target region
					if 0 < int_CB_gReg_pointer < int_CB_gReg_len:
						#check if deletion spans end of reference-sequence
						if int_CB_gReg_pointer+tuple_CB_cigar_tuple[1] < int_CB_gReg_len:
							str_CB_temp_del_sequence = str_CB_gReg_reference_seq[int_CB_gReg_pointer:int_CB_gReg_pointer+tuple_CB_cigar_tuple[1]]
						else:
							str_CB_temp_del_sequence = str_CB_gReg_reference_seq[int_CB_gReg_pointer:] + ('N' * (int_CB_gReg_pointer + tuple_CB_cigar_tuple[1] - int_CB_gReg_len))
						list_CB_base_counts[int_CB_gReg_pointer-1][2].append(str_CB_temp_del_sequence)
					int_CB_gReg_pointer += tuple_CB_cigar_tuple[1]

				#soft clip
				elif tuple_CB_cigar_tuple[0] == 4:
					int_CB_read_pointer += tuple_CB_cigar_tuple[1]
					#remove from final build
					#print('Soft clip!')
				#hard clip
				elif tuple_CB_cigar_tuple[0] == 5:
					pass
					#remove from final build
					#print('Hard clip!')

				else:
					raise ValueError(f'Unexpected CIGAR operation: {tuple_CB_cigar_tuple}')

	return list_CB_base_counts
