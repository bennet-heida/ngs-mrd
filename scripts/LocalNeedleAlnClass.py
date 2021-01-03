class LocalNeedleAln:
	def __init__(self,int_LN_init_max_read_length,int_LN_init_max_query_length,int_LN_init_matchScore=0, \
		int_LN_init_mismatchPen=-1,int_LN_init_gapPen=-1):
		self.int_LN_max_read_length = int_LN_init_max_read_length
		self.int_LN_max_query_length = int_LN_init_max_query_length
		self.int_LN_matchScore = int_LN_init_matchScore
		self.int_LN_mismatchPen = int_LN_init_mismatchPen
		self.int_LN_gapPen = int_LN_init_gapPen

		self.list_LN_array = [[0 for _ in range(self.int_LN_max_read_length+1)] for _ in range(self.int_LN_max_query_length+1)]

		for int_LN_init_j in range(1,self.int_LN_max_read_length+1):
			self.list_LN_array[0][int_LN_init_j] = 0
		for int_LN_init_i in range(1,self.int_LN_max_query_length+1):
			self.list_LN_array[int_LN_init_i][0] = int_LN_init_i*self.int_LN_gapPen



	def alnQuery(self,str_LN_input_read,str_LN_input_query):

		str_LN_read = '-' + str_LN_input_read
		str_LN_query = '-' + str_LN_input_query

		int_LN_len_read = len(str_LN_read)
		int_LN_len_query = len(str_LN_query)

		if int_LN_len_read > self.int_LN_max_read_length+1 or int_LN_len_query > self.int_LN_max_query_length + 1:
			raise ValueError('Read or query too long')

		for int_LN_i in range(1,int_LN_len_query):	
			for int_LN_j in range(1,int_LN_len_read):

				diagPen = self.int_LN_matchScore if str_LN_read[int_LN_j] == str_LN_query[int_LN_i] else self.int_LN_mismatchPen
				self.list_LN_array[int_LN_i][int_LN_j] = max(self.list_LN_array[int_LN_i-1][int_LN_j-1]+diagPen, \
					self.list_LN_array[int_LN_i-1][int_LN_j]+self.int_LN_gapPen,self.list_LN_array[int_LN_i][int_LN_j-1]+self.int_LN_gapPen)


		#backtrack
		int_LN_backtrack_i = int_LN_len_query - 1
		int_LN_max_score = max(self.list_LN_array[int_LN_len_query-1][:int_LN_len_read])
		int_LN_backtrack_j = max(filter(lambda int_LN_lambda_read_pos:self.list_LN_array[int_LN_len_query-1][int_LN_lambda_read_pos]==int_LN_max_score,range(int_LN_len_read)))

		int_LN_align_start = -1
		int_LN_align_end = int_LN_backtrack_j
		int_LN_penalties = 0

		while int_LN_backtrack_i > 0 and int_LN_backtrack_j > 0:
			bool_match = True if str_LN_read[int_LN_backtrack_j] == str_LN_query[int_LN_backtrack_i] else False
			diagPen = self.int_LN_matchScore if bool_match is True else self.int_LN_mismatchPen

			#diagonal
			if self.list_LN_array[int_LN_backtrack_i-1][int_LN_backtrack_j-1]+diagPen == self.list_LN_array[int_LN_backtrack_i][int_LN_backtrack_j]:
				if bool_match is False:
					int_LN_penalties += 1

				int_LN_align_start = int_LN_backtrack_j

				int_LN_backtrack_i -= 1
				int_LN_backtrack_j -= 1

			#deletion on read
			elif self.list_LN_array[int_LN_backtrack_i-1][int_LN_backtrack_j] + self.int_LN_gapPen \
				== self.list_LN_array[int_LN_backtrack_i][int_LN_backtrack_j]:
				int_LN_penalties += 1

				int_LN_align_start = int_LN_backtrack_j + 1

				int_LN_backtrack_i -= 1

			#insertion on read
			elif self.list_LN_array[int_LN_backtrack_i][int_LN_backtrack_j-1] + self.int_LN_gapPen \
				== self.list_LN_array[int_LN_backtrack_i][int_LN_backtrack_j]:
				int_LN_penalties += 1

				int_LN_align_start = int_LN_backtrack_j

				int_LN_backtrack_j -= 1

		int_LN_penalties += int_LN_backtrack_i

		return [int_LN_align_start-1,int_LN_align_end,int_LN_penalties]
