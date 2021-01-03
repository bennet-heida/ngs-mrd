#functions for ngs-mrd scripts
import itertools
import pysam
from collections import defaultdict
import gzip
import csv

def get_analysis_range(str_GAR_gPrimer,str_GAR_MID_F,str_GAR_MID_R,str_GAR_path_MID_file,str_GAR_path_gPrimer_file, \
	int_GAR_length_illumina_read,int_GAR_length_CS,int_GAR_length_RF_tag,int_GAR_end_padding):

	dict_GAR_gPrimer = {}

	with open(str_GAR_path_gPrimer_file,'r') as obj_GAR_gPrimer_file:
		it_GAR_gPrimer_csv = csv.DictReader(obj_GAR_gPrimer_file,delimiter='\t')

		for dict_GAR_temp_gPrimer in it_GAR_gPrimer_csv:

			dict_GAR_gPrimer[dict_GAR_temp_gPrimer['Primer']] = {'chr':dict_GAR_temp_gPrimer['Chr_1'],'Positions':[int(dict_GAR_temp_gPrimer['FirstBase_1']),int(dict_GAR_temp_gPrimer['FirstBase_2'])], \
				'Seq_Lenghts':[len(dict_GAR_temp_gPrimer['Seq_1']),len(dict_GAR_temp_gPrimer['Seq_2'])]}


	list_GAR_MIDs = [{},{}]

	dict_GAR_MID_FR = {'F':0,'R':1}

	with open(str_GAR_path_MID_file,'r') as obj_GAR_MID_file:
		it_GAR_MID_csv = csv.DictReader(obj_GAR_MID_file,delimiter='\t')

		for dict_GAR_temp_MID in it_GAR_MID_csv:
			
			list_GAR_MIDs[dict_GAR_MID_FR[dict_GAR_temp_MID['Stockname'][-1]]][dict_GAR_temp_MID['Primername']] = \
				int(dict_GAR_temp_MID['Spacer_Length']) + len(dict_GAR_temp_MID['MID_Sequence'])


	int_GAR_const_length = int_GAR_length_illumina_read - (int_GAR_length_CS + int_GAR_length_RF_tag + int_GAR_end_padding)
	list_GAR_temp_MIDs = [list_GAR_MIDs[0][str_GAR_MID_F],list_GAR_MIDs[1][str_GAR_MID_R]]


	list_GAR_read_lengths = sorted([(int_GAR_temp_position, int_GAR_const_length - (int_GAR_temp_gPrimer_length + int_GAR_temp_MID_length)) \
		for int_GAR_temp_position,int_GAR_temp_gPrimer_length, int_GAR_temp_MID_length in zip(dict_GAR_gPrimer[str_GAR_gPrimer]['Positions'], \
		dict_GAR_gPrimer[str_GAR_gPrimer]['Seq_Lenghts'],list_GAR_temp_MIDs)], \
		key=lambda tuple_GAR_temp_seq_lenghts:tuple_GAR_temp_seq_lenghts[0])


	int_GAR_start_position = max(list_GAR_read_lengths[0][0],list_GAR_read_lengths[1][0]-list_GAR_read_lengths[1][1])
	int_GAR_end_position = min(list_GAR_read_lengths[1][0],list_GAR_read_lengths[0][0]+list_GAR_read_lengths[0][1])


	return (dict_GAR_gPrimer[str_GAR_gPrimer]['chr'],(int_GAR_start_position,int_GAR_end_position))

def format_time(float_FT_input_seconds):
	int_FT_seconds = int(round(float_FT_input_seconds))

	int_FT_output_hours = int_FT_seconds // 3600
	int_FT_seconds %= 3600
	int_FT_output_minutes = int_FT_seconds // 60
	int_FT_output_seconds = int_FT_seconds % 60

	str_FT_output = (f'{int_FT_output_hours:d}h' if int_FT_output_hours > 0 else '') \
		+ (f'{int_FT_output_minutes:02d}m' if int_FT_output_hours > 0 else f'{int_FT_output_minutes:d}m') + f'{int_FT_output_seconds:02d}s'

	return str_FT_output

def division_zero_tolerant(float_DZT_dividend,float_DZT_divisor):
	if float_DZT_divisor == 0:
		return 0
	else:
		return float_DZT_dividend / float_DZT_divisor

def format_fastq(dict_FF_input):
	str_FF_output = f'@{dict_FF_input["name"]}\n{dict_FF_input["read"]}\n+\n{dict_FF_input["quals"]}\n'
	return str_FF_output

def fastq_zip_equal(*list_ZE_iterables_input):
    obj_ZE_sentinel = object()
    for list_ZE_output in itertools.zip_longest(*list_ZE_iterables_input, fillvalue=obj_ZE_sentinel):
        if any(obj_ZE_temp is obj_ZE_sentinel for obj_ZE_temp in list_ZE_output):
            raise ValueError('Input Files have different lengths')
        yield list_ZE_output

def fastq_iterator(it_FI_input):
	for list_FI_item in it_FI_input:
		#check if names are equal
		if list_FI_item[0][0].rsplit(' ',1)[0] == list_FI_item[1][0].rsplit(' ',1)[0]:
			list_FI_item_output = [{'name':list_FI_temp_read[0],'read':list_FI_temp_read[1],'quals':list_FI_temp_read[2]} \
				for list_FI_temp_read in list_FI_item]
			yield list_FI_item_output
		else:
			raise ValueError('Files out of sync!')

def iter_chunked(it_IC_input,int_IC_chunksize):
	while True:
		list_IC_output = list(itertools.islice(it_IC_input,int_IC_chunksize))
		if len(list_IC_output) < 1:
			break
		yield list_IC_output

def iter_double_chunked(it_IDC_input,int_IDC_inner_chunksize,int_IDC_outer_chunksize):
	return iter_chunked(iter_chunked(it_IDC_input,int_IDC_inner_chunksize),int_IDC_outer_chunksize)

def RevSeq(str_RS_sequence):
    dict_RS_comps = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    str_RS_output = ""
    for str_RS_chr in str_RS_sequence[::-1]:
    	str_RS_output += dict_RS_comps[str_RS_chr]
    return str_RS_output


def input_is_gzip(*str_IG_file_paths):
	#check input file type
	dict_IG_file_endings = {'gz':True, 'gzip':True, 'fq':False, 'fastq':False}
	try:
		list_IG_is_gzip = [dict_IG_file_endings[str_IG_path.rsplit('.',1)[1]] for str_IG_path in str_IG_file_paths]
	except (IndexError, KeyError):
		raise ValueError('Input file endings not recognized')

	if all(list_IG_is_gzip):
		return True
	elif not any(list_IG_is_gzip):
		return False
	else:
		raise ValueError('Input files have different file type')

def open_file_gzip_tolerant(str_OFG_path_to_file,str_OFG_input_mode,bool_OFG_gzip):
	str_OFG_mode = str_OFG_input_mode + ('t' if bool_OFG_gzip else '')
	if bool_OFG_gzip:
		obj_OFG_file_handle = gzip.open(str_OFG_path_to_file,str_OFG_mode)
	else:
		obj_OFG_file_handle = open(str_OFG_path_to_file,str_OFG_mode)
	return obj_OFG_file_handle


def get_RF_distr(str_GRD_input_bam_file,int_GRD_min_members,int_GRD_range_start,int_GRD_range_end,int_GRD_range_step):
	dict_GRD_distr = defaultdict(int)

	with  pysam.AlignmentFile(str_GRD_input_bam_file) as obj_GRD_input_file:

		it_GRD_input_iterator = obj_GRD_input_file.fetch(until_eof=True)

		it_GRD_RF_iterator = filter(lambda obj_GRD_temp_RF:len(obj_GRD_temp_RF[0])>0, itertools.groupby(it_GRD_input_iterator, \
			key=lambda obj_GRD_temp_read:obj_GRD_temp_read.query_name.split('/',1)[0]))

		int_GRD_total_mapped_reads = 0
		int_GRD_total_RF_reads = 0
		int_GRD_RF_above_member_cutoff = 0

		for str_GRD_RF_tag, it_GRD_temp_RF_members in it_GRD_RF_iterator:

			int_GRD_count_members = 0

			for obj_GRD_temp_RF_read in it_GRD_temp_RF_members:
				if not obj_GRD_temp_RF_read.is_unmapped and not obj_GRD_temp_RF_read.is_secondary and not obj_GRD_temp_RF_read.is_supplementary:
					int_GRD_count_members += 1

			int_GRD_total_mapped_reads += int_GRD_count_members
			if int_GRD_count_members >= int_GRD_min_members:
				int_GRD_total_RF_reads += int_GRD_count_members
				int_GRD_RF_above_member_cutoff += 1

			dict_GRD_distr[int_GRD_count_members] += 1

	#list_GRD_output = sorted(dict_GRD_distr.items(),key=lambda tuple_GRD_temp_item:tuple_GRD_temp_item[0])
	list_GRD_output = [(int_GRD_i,dict_GRD_distr[int_GRD_i]) for int_GRD_i in range(int_GRD_range_start,int_GRD_range_end,int_GRD_range_step)]

	return int_GRD_total_RF_reads,int_GRD_total_mapped_reads,int_GRD_RF_above_member_cutoff,list_GRD_output


def fastq_get_read_number(*str_FGN_files):
	##Get Number of Reads in files, do basic integrity check:
	bool_FGN_input_gzip = input_is_gzip(*str_FGN_files)
	list_FGN_read_counts = []

	for str_FGN_temp_file in str_FGN_files:
		with open_file_gzip_tolerant(str_FGN_temp_file,'r',bool_FGN_input_gzip) as obj_FGN_file:
			int_FGN_line_count = 0
			for str_FGN_line in obj_FGN_file: int_FGN_line_count += 1
			#check if file length is multiple of 4, for 4 lines in a fastq entry
			if int_FGN_line_count % 4 == 0:
				list_FGN_read_counts.append(int_FGN_line_count // 4)
			else:
				raise ValueError(f'File {str_FGN_temp_file} line count is not multiple of 4')
	#check if both files have same length
	int_FGN_output_read_count = list_FGN_read_counts[0]
	if list_FGN_read_counts.count(int_FGN_output_read_count) == len(str_FGN_files):
		return int_FGN_output_read_count
	else:
		raise ValueError('Input Files don\'t have the same length')

def Chr2Phred(str_c2p_input_chr):
	if len(str_c2p_input_chr) != 1:
		raise ValueError('expects string of length 1')
	return ord(str_c2p_input_chr[0])-33

def Phred2Chr(int_p2c_input_phred):
	return chr(int_p2c_input_phred+33)

def get_average_phred(list_GA_input):
	return int(round(sum(list_GA_input)/len(list_GA_input)))