from NGSMRDFuncs import Chr2Phred, Phred2Chr,get_average_phred

def compare_phreds(str_CP_phred_gap,str_CP_phred_surrounding):
	list_CP_phred_gap = [Chr2Phred(str_CP_phred_gap[int_CP_n]) for int_CP_n in range(len(str_CP_phred_gap))]
	list_CP_phred_surrounding = [Chr2Phred(str_CP_phred_surrounding[int_CP_n]) for int_CP_n in range(len(str_CP_phred_surrounding))]

	int_CP_average_phred_surrounding = get_average_phred(list_CP_phred_surrounding)

	if get_average_phred(list_CP_phred_gap) >= int_CP_average_phred_surrounding:
		list_CP_new_phreds = [(int_CP_temp_phred - int_CP_average_phred_surrounding if int_CP_temp_phred - int_CP_average_phred_surrounding > 0 else 0) for int_CP_temp_phred in list_CP_phred_gap]
		str_CP_new_phreds = "".join([Phred2Chr(int_CP_temp_phred) for int_CP_temp_phred in list_CP_new_phreds])
		return [True,str_CP_new_phreds]
	else:
		return[False,'']

class LocalGotohMerge:
	def __init__(self,int_LG_init_read1_max_length,int_LG_init_read2_max_length,\
		int_LG_init_matchScore=1,int_LG_init_mismatchPen=-1,int_LG_init_gapOpenPen=-3,int_LG_init_gapExtendPen=-1):
		
		self.int_LG_read1_max_length = int_LG_init_read1_max_length
		self.int_LG_read2_max_length = int_LG_init_read2_max_length
		self.int_LG_matchScore = int_LG_init_matchScore
		self.int_LG_mismatchPen = int_LG_init_mismatchPen
		self.int_LG_gapOpenPen = int_LG_init_gapOpenPen
		self.int_LG_gapExtendPen = int_LG_init_gapExtendPen

		self.list_LG_arrayM = [[0 for _ in range(self.int_LG_read1_max_length+1)] for _ in range(self.int_LG_read2_max_length+1)]
		self.list_LG_arrayUp = [[0 for _ in range(self.int_LG_read1_max_length+1)] for _ in range(self.int_LG_read2_max_length+1)]
		self.list_LG_arrayLeft = [[0 for _ in range(self.int_LG_read1_max_length+1)] for _ in range(self.int_LG_read2_max_length+1)]


		for int_LG_i in range(self.int_LG_read2_max_length+1):
			self.list_LG_arrayUp[int_LG_i][0] = 0 #start gaps are free
			self.list_LG_arrayM[int_LG_i][0] = 0 #start gaps are free
			self.list_LG_arrayLeft[int_LG_i][0] = float('-inf')
		for int_LG_j in range(self.int_LG_read1_max_length+1):
			self.list_LG_arrayLeft[0][int_LG_j] = 0 #start gaps are free
			self.list_LG_arrayM[0][int_LG_j] = 0 #start gaps are free
			self.list_LG_arrayUp[0][int_LG_j] = float('-inf')



	def alnMerge(self,list_LG_input_read1,list_LG_input_read2):
		
		self.list_LG_read1 = list_LG_input_read1
		self.list_LG_read2 = list_LG_input_read2
		#add first character for convenience

		for int_LG_num in range(2):
			self.list_LG_read1[int_LG_num] = '-' + self.list_LG_read1[int_LG_num]
			self.list_LG_read2[int_LG_num] = '-' + self.list_LG_read2[int_LG_num]

		self.int_LG_length_read1 = len(self.list_LG_read1[0])
		self.int_LG_length_read2 = len(self.list_LG_read2[0])
		
		if self.int_LG_length_read1 > self.int_LG_read1_max_length+1 or self.int_LG_length_read2 > self.int_LG_read2_max_length+1:
			raise ValueError('One of the reads exceeds max length')


		self.int_LG_maxScore_i = 0
		self.int_LG_maxScore_j = 0
		self.int_LG_maxScore = float('-inf')
		#calculate matrix
		for int_LG_i in range(1,self.int_LG_length_read2):
			for int_LG_j in range(1,self.int_LG_length_read1):

				int_LG_diagPen = self.int_LG_matchScore if self.list_LG_read1[0][int_LG_j] == self.list_LG_read2[0][int_LG_i] else self.int_LG_mismatchPen

				self.list_LG_arrayUp[int_LG_i][int_LG_j] = max(self.list_LG_arrayM[int_LG_i-1][int_LG_j]+self.int_LG_gapOpenPen, \
					self.list_LG_arrayUp[int_LG_i-1][int_LG_j]+self.int_LG_gapExtendPen)
				self.list_LG_arrayLeft[int_LG_i][int_LG_j] = max(self.list_LG_arrayM[int_LG_i][int_LG_j-1]+self.int_LG_gapOpenPen, \
					self.list_LG_arrayLeft[int_LG_i][int_LG_j-1]+self.int_LG_gapExtendPen)
				self.list_LG_arrayM[int_LG_i][int_LG_j] = max(self.list_LG_arrayM[int_LG_i-1][int_LG_j-1]+int_LG_diagPen, \
					self.list_LG_arrayUp[int_LG_i][int_LG_j],self.list_LG_arrayLeft[int_LG_i][int_LG_j])
		

				if (int_LG_i == self.int_LG_length_read2 - 1 or int_LG_j == self.int_LG_length_read1 - 1) \
					and self.list_LG_arrayM[int_LG_i][int_LG_j] >= self.int_LG_maxScore:
					self.int_LG_maxScore_i = int_LG_i
					self.int_LG_maxScore_j = int_LG_j
					self.int_LG_maxScore = self.list_LG_arrayM[int_LG_i][int_LG_j]


		#backtrack
		int_LG_backtrack_i = self.int_LG_maxScore_i
		int_LG_backtrack_j = self.int_LG_maxScore_j

		list_LG_consensus = ['','']

		while int_LG_backtrack_i > 0 and int_LG_backtrack_j > 0:
			bool_match = True if self.list_LG_read1[0][int_LG_backtrack_j] == self.list_LG_read2[0][int_LG_backtrack_i] else False
			int_LG_diagPen = self.int_LG_matchScore if bool_match else self.int_LG_mismatchPen

			#backtrack as match
			if self.list_LG_arrayM[int_LG_backtrack_i][int_LG_backtrack_j] \
				== self.list_LG_arrayM[int_LG_backtrack_i - 1][int_LG_backtrack_j - 1] + int_LG_diagPen:
				if bool_match is True:
					list_LG_consensus[0] = self.list_LG_read1[0][int_LG_backtrack_j] + list_LG_consensus[0]
					list_LG_consensus[1] = Phred2Chr(max(Chr2Phred(self.list_LG_read1[1][int_LG_backtrack_j]),Chr2Phred(self.list_LG_read2[1][int_LG_backtrack_i]))) \
						+ list_LG_consensus[1]
				else:
					list_LG_consensus[0] = 'N' + list_LG_consensus[0]
					list_LG_consensus[1] = Phred2Chr(min(Chr2Phred(self.list_LG_read1[1][int_LG_backtrack_j]),Chr2Phred(self.list_LG_read2[1][int_LG_backtrack_i]))) \
						+ list_LG_consensus[1]

				int_LG_backtrack_i -= 1
				int_LG_backtrack_j -= 1

			#deletion on read1, insert on read2
			elif self.list_LG_arrayM[int_LG_backtrack_i][int_LG_backtrack_j] \
				== self.list_LG_arrayUp[int_LG_backtrack_i][int_LG_backtrack_j]:
				for int_LG_k in range(1,int_LG_backtrack_i+1):
					if self.list_LG_arrayM[int_LG_backtrack_i - int_LG_k][int_LG_backtrack_j] + self.int_LG_gapOpenPen \
						> self.list_LG_arrayUp[int_LG_backtrack_i - int_LG_k][int_LG_backtrack_j] + self.int_LG_gapExtendPen:
						int_LG_gap_steps = int_LG_k
						break
				
				#calculate gap phred quality scores
				list_LG_gap_phreds = compare_phreds(self.list_LG_read2[1][int_LG_backtrack_i-int_LG_gap_steps+1:int_LG_backtrack_i+1], \
					self.list_LG_read1[1][int_LG_backtrack_j:int_LG_backtrack_j+2])
				if list_LG_gap_phreds[0] is True:
					list_LG_consensus[0] = 'N'*int_LG_gap_steps + list_LG_consensus[0]
					list_LG_consensus[1] = list_LG_gap_phreds[1] + list_LG_consensus[1]

				int_LG_backtrack_i -= int_LG_gap_steps


			#deletion on read2, insert on read1
			elif self.list_LG_arrayM[int_LG_backtrack_i][int_LG_backtrack_j] \
				== self.list_LG_arrayLeft[int_LG_backtrack_i][int_LG_backtrack_j]:
				for int_LG_k in range(1,int_LG_backtrack_j+1):
					if self.list_LG_arrayM[int_LG_backtrack_i][int_LG_backtrack_j - int_LG_k] + self.int_LG_gapOpenPen \
						> self.list_LG_arrayUp[int_LG_backtrack_i][int_LG_backtrack_j - int_LG_k] + self.int_LG_gapExtendPen:
						int_LG_gap_steps = int_LG_k
						break
				
				#calculate gap phred quality scores
				list_LG_gap_phreds = compare_phreds(self.list_LG_read1[1][int_LG_backtrack_j-int_LG_gap_steps+1:int_LG_backtrack_j+1], \
					self.list_LG_read2[1][int_LG_backtrack_i:int_LG_backtrack_i+2])
				if list_LG_gap_phreds[0] is True:
					list_LG_consensus[0] = 'N'*int_LG_gap_steps + list_LG_consensus[0]
					list_LG_consensus[1] = list_LG_gap_phreds[1] + list_LG_consensus[1]

				int_LG_backtrack_j -= int_LG_gap_steps

		return [self.int_LG_maxScore,list_LG_consensus]
