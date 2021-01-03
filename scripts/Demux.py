#Script for demultiplexing reads from NGS-MRD run on MiSeq
#to run this script 'process_files()' function should be run. This function expects the input as follows:
#
#list_PF_path_input_files: [path_input_file1,path_input_file2] list of paths to input fastq files with forward and reverse reads,
#	reads must be in same order, file can be compressed with gzip (.fastq.gz) or plain fastq
#list_PF_output_samples: 2-dimensional list, each row represents sample id and path to output file without file ending
#i.e. [[1, '/home/user/MRD2019_015/Demux/MRD2019_015_01_IDH2_NGS_145_e20210256'], [2, '/home/user/MRD2019_015/Demux/MRD2019_015_02_RUNX1_NGS_175_N786'],...]
#str_PF_output_prefix: path and file prefix to use for default output (unmatched, empty, ambiguous)
#i.e. '/home/user/MRD2019_015/Demux/MRD2019_015'
#
#list_PF_MIDs: contains two lists, one for forward and one for reverse MIDs, each item is a dictionary containing the MID info for one MID-primer (ID, seq, start)
#i.e. [[{'ID': 'MID1', 'seq': 'ACGAGTGCGT', 'start': 1}, {'ID': 'MID3', 'seq': 'AGACGCACTC', 'start': 2},...],[{'ID': 'MID2', 'seq': 'ACGCTCGACA', 'start': 1},...]]
#
#dict_PF_MID_grid: 2-dimensional dictionary, representing the MID-assignment for the run, primary keys being the forward MIDs, secondary keys
#are the reverse MIDs
#i.e. {'MID1': {'MID10': 9, 'MID12': 2, 'MID14': 4, 'MID16': 19}, 'MID3': {'MID10': 15, 'MID12': 3, 'MID14': 5, 'MID16': 22},...}
#
#int_PF_MIDs_maxPenalties: max. number of penalties (mismatch, insertion, deletion) allowed in MID-Sequence to still be considered a match
#bool_PF_check_MID_pos: True if aligned position of MID should be considered for max. Penalties
#
#int_PF_read_chunk: number of reads to be read from file and processed at a time, defaults to 10.000
#int_PF_number_processes: number of processes that should be run in parallel, 0 defaults to the number of CPUs on the system
#

#process_files(list_PF_path_input_files,list_PF_output_samples,str_PF_output_prefix,list_PF_MIDs,dict_PF_MID_grid, \
#	int_PF_MIDs_maxPenalties,bool_PF_check_MID_pos,int_PF_read_chunk=10000,int_PF_number_processes=0):

import sys
import os
from time import monotonic
import contextlib
from LocalNeedleAlnClass import LocalNeedleAln
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import concurrent.futures
import itertools
from operator import itemgetter
from NGSMRDFuncs import fastq_iterator, fastq_zip_equal, fastq_get_read_number, input_is_gzip, open_file_gzip_tolerant, iter_double_chunked, format_fastq, format_time

def process_read_chunk(list_PRC_reads,list_PRC_MIDs,dict_PRC_MID_grid,int_PRC_MIDs_maxPenalties,bool_PRC_check_MID_pos):
	#calc max. array size for local needle
	int_PRC_max_MID_len = max([max([len(dict_PRC_temp_MID['seq']) for dict_PRC_temp_MID in list_PRC_temp_MIDs]) for list_PRC_temp_MIDs in list_PRC_MIDs])

	obj_PRC_LocalNeedle = LocalNeedleAln(int_PRC_max_MID_len + 2*int_PRC_MIDs_maxPenalties, int_PRC_max_MID_len)

	list_PRC_output_reads = list(itertools.starmap(process_read_pair,zip(*map(itertools.repeat,[obj_PRC_LocalNeedle,list_PRC_MIDs,dict_PRC_MID_grid, \
		int_PRC_MIDs_maxPenalties,bool_PRC_check_MID_pos]),list_PRC_reads)))

	return list_PRC_output_reads

def process_read_pair(obj_PRP_LocalNeedle,list_PRP_MIDs,dict_PRP_MID_grid,int_PRP_MIDs_maxPenalties,bool_PRP_check_MID_pos,list_PRP_read):
	list_PRP_MID_matches = list(itertools.starmap(find_MIDs,zip(*map(itertools.repeat,[obj_PRP_LocalNeedle,int_PRP_MIDs_maxPenalties,bool_PRP_check_MID_pos]), \
		list_PRP_read,list_PRP_MIDs)))

	list_PRP_count_mid_match = list(map(len,list_PRP_MID_matches))

	list_PRP_cut_points = [0,0]

	int_PRP_sample_ID = None

	if all(int_PRP_temp_count_mid_match == 1 for int_PRP_temp_count_mid_match in list_PRP_count_mid_match):
		int_PRP_sample_ID = dict_PRP_MID_grid.get(list_PRP_MID_matches[0][0][0],{}).get(list_PRP_MID_matches[1][0][0],None)
		int_PRP_demux_state = 1 if int_PRP_sample_ID is not None else 2

		list_PRP_cut_points = [list_PRP_MID_matches[0][0][1],list_PRP_MID_matches[1][0][1]]

	elif any(int_PRP_temp_count_mid_match == 0 for int_PRP_temp_count_mid_match in list_PRP_count_mid_match):
		int_PRP_demux_state = 0
	else:
		int_PRP_demux_state = 3

	str_PRP_output_MID_match = "-".join(map(itemgetter(0),list_PRP_MID_matches[0]))+":"+"-".join(map(itemgetter(0),list_PRP_MID_matches[1]))

	list_PRP_output_read = [{'name':str_PRP_output_MID_match + "/" + dict_PRP_temp_read['name'],'read':dict_PRP_temp_read['read'][int_PRP_temp_cut_point:], \
	'quals':dict_PRP_temp_read['quals'][int_PRP_temp_cut_point:]} for dict_PRP_temp_read, int_PRP_temp_cut_point in zip(list_PRP_read,list_PRP_cut_points)]

	return (int_PRP_demux_state,int_PRP_sample_ID,list_PRP_output_read)

def find_MIDs(obj_FM_LocalNeedle,int_FM_MIDs_maxPenalties,bool_FM_check_MID_pos,dict_FM_read,list_FM_MIDs):
	#try exact match
	list_FM_matches = [(dict_FM_temp_MID['ID'],dict_FM_temp_MID['start']+len(dict_FM_temp_MID['seq'])) for dict_FM_temp_MID in list_FM_MIDs if dict_FM_temp_MID['seq'] == \
	dict_FM_read['read'][dict_FM_temp_MID['start']:dict_FM_temp_MID['start']+len(dict_FM_temp_MID['seq'])]]

	#if no matches were found, try error tolerant match
	if len(list_FM_matches) < 1:
		list_FM_matches = [(str_FM_temp_MID_ID,int_FM_temp_MID_cut_point) for str_FM_temp_MID_ID,(bool_FM_temp_MID_match, int_FM_temp_MID_cut_point) \
		in zip([dict_FM_temp_MID['ID'] for dict_FM_temp_MID in list_FM_MIDs],itertools.starmap(mid_match_needle,zip(*map(itertools.repeat, \
			[obj_FM_LocalNeedle,int_FM_MIDs_maxPenalties,dict_FM_read,bool_FM_check_MID_pos]),list_FM_MIDs))) if bool_FM_temp_MID_match]


	return list_FM_matches

def mid_match_needle(obj_MMN_LocalNeedle,int_MMN_MIDs_maxPenalties,dict_MMN_read,bool_MMN_check_MID_pos,dict_MMN_MID):
	list_MMN_expected_match_pos = [dict_MMN_MID['start'],dict_MMN_MID['start'] + len(dict_MMN_MID['seq'])]
	int_MMN_read_start_pos = max(0,dict_MMN_MID['start'] - int_MMN_MIDs_maxPenalties)

	list_MMN_needle_match = obj_MMN_LocalNeedle.alnQuery(dict_MMN_read['read'][int_MMN_read_start_pos:list_MMN_expected_match_pos[1] + \
		int_MMN_MIDs_maxPenalties],dict_MMN_MID['seq'])

	int_MMN_final_pen_count = list_MMN_needle_match[2]
	#calculate distance to expected match position
	if bool_MMN_check_MID_pos:
		int_MMN_final_pen_count += min([abs(int_MMN_temp_expected - int_MMN_temp_found) for int_MMN_temp_expected, int_MMN_temp_found in \
			zip(list_MMN_expected_match_pos,list_MMN_needle_match[:2])])

	bool_MMN_match = True if int_MMN_final_pen_count <= int_MMN_MIDs_maxPenalties else False
	int_MMN_cut_point = int_MMN_read_start_pos + list_MMN_needle_match[1] if int_MMN_final_pen_count <= int_MMN_MIDs_maxPenalties else None

	return (bool_MMN_match,int_MMN_cut_point)

def process_files(list_PF_path_input_files,list_PF_output_samples,str_PF_output_prefix,list_PF_MIDs,dict_PF_MID_grid, \
	int_PF_MIDs_maxPenalties,bool_PF_check_MID_pos,int_PF_read_chunk=10000,int_PF_number_processes=0):
	#number processes defaults to cpu count
	int_PF_number_processes = int_PF_number_processes if int_PF_number_processes > 0 else os.cpu_count()

	float_PF_time_start = monotonic()
	print('Open input file: ' + list_PF_path_input_files[0])
	print('Open input file: ' + list_PF_path_input_files[1])
	int_PF_total_reads = fastq_get_read_number(*list_PF_path_input_files)
	bool_PF_input_gzip = input_is_gzip(*list_PF_path_input_files)

	#open files
	list_PF_output_file_endings = ['_F.fq','_R.fq']
	list_PF_output_default_names = [[0,'_Unmatched'],[2,'_Empty'],[3,'_Ambiguous']]
	with contextlib.ExitStack() as stack:
		#Input files
		list_PF_input_files = [stack.enter_context(open_file_gzip_tolerant(str_PF_temp_path_input_file, 'r', bool_PF_input_gzip)) \
			for str_PF_temp_path_input_file in list_PF_path_input_files]
		#create fastq iterator for input files
		it_PF_input_fastq = fastq_iterator(fastq_zip_equal(*[FastqGeneralIterator(obj_PF_temp_input_file) \
			for obj_PF_temp_input_file in list_PF_input_files]))

		#open output files
		dict_PF_output_sample_files = {int_PF_temp_sample_ID:[stack.enter_context(open(str_PF_temp_path_output_file + str_PF_temp_output_file_ending, 'w')) \
			for str_PF_temp_output_file_ending in list_PF_output_file_endings] for int_PF_temp_sample_ID, str_PF_temp_path_output_file in list_PF_output_samples}
		#open default output files
		dict_PF_output_default_files = {int_PF_temp_default_ID:[stack.enter_context(open(str_PF_output_prefix + str_PF_temp_output_default_name + str_PF_temp_output_file_ending, 'w')) \
			for str_PF_temp_output_file_ending in list_PF_output_file_endings] for int_PF_temp_default_ID,str_PF_temp_output_default_name in list_PF_output_default_names}

		int_PF_reads_processed = 0


		for list_PF_reads_for_processing in iter_double_chunked(it_PF_input_fastq,int_PF_read_chunk,int_PF_number_processes):
			with concurrent.futures.ProcessPoolExecutor() as executor:
				list_PF_reads_processed = executor.map(process_read_chunk,list_PF_reads_for_processing,*(itertools.repeat(obj_PF_temp_arg) \
					for obj_PF_temp_arg in (list_PF_MIDs,dict_PF_MID_grid,int_PF_MIDs_maxPenalties,bool_PF_check_MID_pos)))
				#write processed reads to files
				for tup_PF_temp_read_output in itertools.chain.from_iterable(list_PF_reads_processed):
					if tup_PF_temp_read_output[0] == 1:
						for int_PF_temp_read_FR, dict_PF_temp_read in enumerate(tup_PF_temp_read_output[2]):
							dict_PF_output_sample_files[tup_PF_temp_read_output[1]][int_PF_temp_read_FR].write(format_fastq(dict_PF_temp_read))
					else:
						for int_PF_temp_read_FR, dict_PF_temp_read in enumerate(tup_PF_temp_read_output[2]):
							dict_PF_output_default_files[tup_PF_temp_read_output[0]][int_PF_temp_read_FR].write(format_fastq(dict_PF_temp_read))
					int_PF_reads_processed += 1


			str_PF_percent_processed = f'{round(int_PF_reads_processed*100/int_PF_total_reads,2):.2f}%' \
				if int_PF_total_reads > 0 else '100.00%'

			sys.stdout.write(f'\rProcessing Reads: {int_PF_reads_processed} / {int_PF_total_reads} ({str_PF_percent_processed})')
			sys.stdout.flush()


	float_PF_time_total = monotonic() - float_PF_time_start
	print(f'\nEnd of files: processed {int_PF_reads_processed} reads in {format_time(float_PF_time_total)}')


def main():
	print('No main function, please run \'process_files()\' function with arguments.')

if __name__ == '__main__':
	main()