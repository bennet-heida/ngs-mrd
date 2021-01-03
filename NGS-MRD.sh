#!/bin/bash
set -e
MRD_FOLDER="/home/user/MRD/"

RUN_NAME="MRD2019_000"

DEMUX_PATH_FILE1=${MRD_FOLDER}${RUN_NAME}"/Fastq/MRD2019_000_R1_001.fastq.gz"
DEMUX_PATH_FILE2=${MRD_FOLDER}${RUN_NAME}"/Fastq/MRD2019_000_R2_001.fastq.gz"

RUN_INFO_FILE="/home/user/MRD_Pipe/runs/"${RUN_NAME}".tsv"
MID_SEQS_FILE="/home/user/MRD_Pipe/seqs/MID_seqs.tsv"
GPRIMER_SEQS_FILE="/home/user/MRD_Pipe/seqs/gPrimer_seqs.tsv"
REF_FASTA_FILE="/home/user/ref_fasta/grch37.fa"
OUTPUT_FILE=${MRD_FOLDER}${RUN_NAME}"/MRDScores/"${RUN_NAME}".xlsx"

PILEUP_RF_QCUT="qcut"

DEMUX_OUTPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Demux/"
DEMUX_MID_SEQS=${MID_SEQS_FILE}
DEMUX_MAXPENMIDS=1
DEMUX_CHECK_MID_POS=1
DEMUX_READCHUNKS=10000
DEMUX_NUMBERPROCESSES=0

CUT_PRIMER_SEQS=${GPRIMER_SEQS_FILE}
CUT_INPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Demux/"
CUT_OUTPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Cut/"
CUT_MAXPENCS=3
CUT_MAXPENGPRIMER=3
CUT_MAXPENGPRIMERREV=5
CUT_BUFFERCS=5
CUT_BUFFERGPRIMER=5
CUT_BUFFERGPRIMERREV=10
CUT_READCHUNKS=10000

QFILTER_INPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Cut/"
QFILTER_OUTPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Cut/"
QFILTER_PHRED_CUTOFF=32
QFILTER_PERCENT_LOW_Q=10
QFILTER_PERCENT_MAX_READ_FILTER=20
QFILTER_START_END_TRIM=50
QFILTER_MIN_END_TRIM_LENGTH=10
QFILTER_END_TRIM_PERCENT_LOW_Q=50
QFILTER_READCHUNKS=10000

MERGE_INPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Cut/"
MERGE_OUTPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/R1R2/"
MERGE_MIN_CONSENSUS_LENGTH=40
MERGE_MIN_SCORE_FRACTION=60
MERGE_READCHUNKS=5000

PILEUP_R_READ_FILES_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Cut/"

PILEUP_R1R2_READ_FILES_FOLDER=${MRD_FOLDER}${RUN_NAME}"/R1R2/"
PILEUP_R1R2_INPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Pileup/"
PILEUP_R1R2_OUTPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Pileup/"
PILEUP_R1R2_REF_FASTA_FILE=${REF_FASTA_FILE}
PILEUP_R1R2_GPRIMER_SEQS_FILE=${GPRIMER_SEQS_FILE}


PILEUP_RF_READ_FILES_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Cut/"
PILEUP_RF_INPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Pileup/"
PILEUP_RF_OUTPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Pileup/"
PILEUP_RF_REF_FASTA_FILE=${REF_FASTA_FILE}
PILEUP_RF_GPRIMER_SEQS_FILE=${GPRIMER_SEQS_FILE}
PILEUP_RF_MIN_RF_MEMBERS=3
PILEUP_RF_PERCENT_CONSENSUS=70

WRITE_XLSX_PILEUP_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Pileup/"
WRITE_XLSX_FASTQ_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Demux/"
WRITE_XLSX_ALN_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Pileup/"
WRITE_XLSX_OUTPUT_FILE=${OUTPUT_FILE}

PLOT_OUTPUT_FOLDER=${MRD_FOLDER}${RUN_NAME}"/Plots/"

SCRIPT_PATH=$(dirname $(realpath -s $0))

mkdir -p $DEMUX_OUTPUT_FOLDER

python3 $SCRIPT_PATH/scripts/DemuxRun.py --input_file1 $DEMUX_PATH_FILE1 --input_file2 $DEMUX_PATH_FILE2 --output_folder $DEMUX_OUTPUT_FOLDER --run_name $RUN_NAME \
--run_info_file $RUN_INFO_FILE --MID_seq_file $DEMUX_MID_SEQS --maxPen_MID $DEMUX_MAXPENMIDS \
--check_MID_pos $DEMUX_CHECK_MID_POS --read_chunks $DEMUX_READCHUNKS

mkdir -p $CUT_OUTPUT_FOLDER

python3 $SCRIPT_PATH/scripts/CutAdaptersRun.py --run_info_file $RUN_INFO_FILE --run_name $RUN_NAME --primer_info_file $CUT_PRIMER_SEQS --input_folder $CUT_INPUT_FOLDER \
--output_folder $CUT_OUTPUT_FOLDER --maxPen_CS $CUT_MAXPENCS --maxPen_gPrimer $CUT_MAXPENGPRIMER --maxPen_gPrimer_rev $CUT_MAXPENGPRIMERREV --buffer_CS $CUT_BUFFERCS \
--buffer_gPrimer $CUT_BUFFERGPRIMER --buffer_gPrimer_rev $CUT_BUFFERGPRIMERREV --read_chunks $CUT_READCHUNKS

mkdir -p $QFILTER_OUTPUT_FOLDER

python3 $SCRIPT_PATH/scripts/QFilterRun.py --run_name $RUN_NAME --run_info_file $RUN_INFO_FILE --input_folder $QFILTER_INPUT_FOLDER --output_folder $QFILTER_OUTPUT_FOLDER \
--phred_cutoff $QFILTER_PHRED_CUTOFF --percent_low_q $QFILTER_PERCENT_LOW_Q --percent_max_read_filter $QFILTER_PERCENT_MAX_READ_FILTER --start_end_trim $QFILTER_START_END_TRIM \
--min_end_trim_length $QFILTER_MIN_END_TRIM_LENGTH --end_trim_percent_low_q $QFILTER_END_TRIM_PERCENT_LOW_Q --read_chunks $QFILTER_READCHUNKS

mkdir -p $MERGE_OUTPUT_FOLDER

python3 $SCRIPT_PATH/scripts/MergeRun.py --run_info_file $RUN_INFO_FILE --run_name $RUN_NAME --MID_seq_file $DEMUX_MID_SEQS --primer_info_file $CUT_PRIMER_SEQS --input_folder $MERGE_INPUT_FOLDER \
--output_folder $MERGE_OUTPUT_FOLDER --min_consensus_length $MERGE_MIN_CONSENSUS_LENGTH --min_score_fraction $MERGE_MIN_SCORE_FRACTION --read_chunks $MERGE_READCHUNKS


mkdir -p $PILEUP_R1R2_INPUT_FOLDER
mkdir -p $PILEUP_R1R2_OUTPUT_FOLDER

mkdir -p $PILEUP_RF_INPUT_FOLDER
mkdir -p $PILEUP_RF_OUTPUT_FOLDER

###align fastq files with bwa

#R1 and R2 - Reads
for Pileup_R_FPath in ${PILEUP_R_READ_FILES_FOLDER}*"_cut_F.fq"; do
Pileup_R_FBase=${Pileup_R_FPath##*/}
Pileup_R_FBase=${Pileup_R_FBase%_cut_F.fq}

bwa aln $PILEUP_R1R2_REF_FASTA_FILE ${PILEUP_R_READ_FILES_FOLDER}${Pileup_R_FBase}"_cut_F.fq" \
| bwa samse $PILEUP_R1R2_REF_FASTA_FILE - ${PILEUP_R_READ_FILES_FOLDER}${Pileup_R_FBase}"_cut_F.fq" \
| samtools view -S -b - | samtools sort -n - -o ${PILEUP_R1R2_INPUT_FOLDER}${Pileup_R_FBase}"_R1.bam"

bwa aln $PILEUP_R1R2_REF_FASTA_FILE ${PILEUP_R_READ_FILES_FOLDER}${Pileup_R_FBase}"_cut_R.fq" \
| bwa samse $PILEUP_R1R2_REF_FASTA_FILE - ${PILEUP_R_READ_FILES_FOLDER}${Pileup_R_FBase}"_cut_R.fq" \
| samtools view -S -b - | samtools sort -n - -o ${PILEUP_R1R2_INPUT_FOLDER}${Pileup_R_FBase}"_R2.bam"

done

##R1R2-Reads
for Pileup_R1R2_FPath in ${PILEUP_R1R2_READ_FILES_FOLDER}*; do
Pileup_R1R2_FBase=${Pileup_R1R2_FPath##*/}
Pileup_R1R2_FBase=${Pileup_R1R2_FBase%_R1R2.fq}

bwa aln $PILEUP_R1R2_REF_FASTA_FILE ${PILEUP_R1R2_READ_FILES_FOLDER}${Pileup_R1R2_FBase}"_R1R2.fq" \
| bwa samse $PILEUP_R1R2_REF_FASTA_FILE - ${PILEUP_R1R2_READ_FILES_FOLDER}${Pileup_R1R2_FBase}"_R1R2.fq" \
| samtools view -S -b - | samtools sort -n - -o ${PILEUP_R1R2_INPUT_FOLDER}${Pileup_R1R2_FBase}"_R1R2.bam"

done

##RF-Reads
for Pileup_RF_FPath in ${PILEUP_RF_READ_FILES_FOLDER}*"_"${PILEUP_RF_QCUT}"_F.fq"; do
Pileup_RF_FBase=${Pileup_RF_FPath##*/}
Pileup_RF_FBase=${Pileup_RF_FBase%"_"${PILEUP_RF_QCUT}"_F.fq"}

samtools cat ${PILEUP_R1R2_INPUT_FOLDER}${Pileup_RF_FBase}"_R1.bam" ${PILEUP_R1R2_INPUT_FOLDER}${Pileup_RF_FBase}"_R2.bam" \
| samtools sort -n - -o ${PILEUP_RF_INPUT_FOLDER}${Pileup_RF_FBase}"_RF.bam"

done

######


python3 $SCRIPT_PATH/scripts/PileupReadsRun.py --run_info_file $RUN_INFO_FILE --run_name $RUN_NAME --ref_fasta_file $PILEUP_R1R2_REF_FASTA_FILE \
--primer_info_file $PILEUP_R1R2_GPRIMER_SEQS_FILE --input_folder $PILEUP_R1R2_INPUT_FOLDER --output_folder $PILEUP_R1R2_OUTPUT_FOLDER



python3 $SCRIPT_PATH/scripts/PileupRFRun.py --run_info_file $RUN_INFO_FILE --run_name $RUN_NAME --ref_fasta_file $PILEUP_RF_REF_FASTA_FILE \
--primer_info_file $PILEUP_RF_GPRIMER_SEQS_FILE --input_folder $PILEUP_RF_INPUT_FOLDER \
--output_folder $PILEUP_RF_OUTPUT_FOLDER --min_member_count $PILEUP_RF_MIN_RF_MEMBERS --min_consensus_percent $PILEUP_RF_PERCENT_CONSENSUS


mkdir -p ${WRITE_XLSX_OUTPUT_FILE%/*}

##write to output file
python3 $SCRIPT_PATH/scripts/WriteXLSX.py --run_name $RUN_NAME --run_info_file $RUN_INFO_FILE --MID_file $DEMUX_MID_SEQS --gPrimer_file $CUT_PRIMER_SEQS \
--pileup_folder $WRITE_XLSX_PILEUP_FOLDER --fastq_folder $WRITE_XLSX_FASTQ_FOLDER --aln_folder $WRITE_XLSX_ALN_FOLDER --output_file $WRITE_XLSX_OUTPUT_FILE

mkdir -p $PLOT_OUTPUT_FOLDER

##plot distribution
python3 $SCRIPT_PATH/scripts/PlotRun.py --run_name $RUN_NAME --run_info_file $RUN_INFO_FILE --MID_file $DEMUX_MID_SEQS --gPrimer_file $CUT_PRIMER_SEQS \
--pileup_folder $WRITE_XLSX_PILEUP_FOLDER --output_folder $PLOT_OUTPUT_FOLDER

##zip files
gzip ${DEMUX_OUTPUT_FOLDER}*
gzip ${CUT_OUTPUT_FOLDER}*
gzip ${MERGE_OUTPUT_FOLDER}*
