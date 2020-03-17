import regex as re
from collections import Counter

file_name = input('Enter multi-fasta file')

def getSortedTuplesFromDict(seq_dict,inc=True):
	'''
	Returns [(seq_id,sequence_string)] in sorted fashion indicated by inc.
	'''
	return sorted(seq_dict.items(),key=lambda x: len(x[1]), reverse=not inc)

def getOpenReadingFrames(seq_dict,frame,inc=True):
	'''
	Returns {<seq_id>:(start,<longest orf string>)}
	'''
	seq_orf_dict = {}
	for seq in seq_dict:
		seq_orf_dict[seq] = getLongestORF(seq_dict[seq],frame)

	##return sorted(seq_orf_dict.items(),key=lambda x: len(x[1][1]),reverse=not inc)
	return seq_orf_dict

def getLongestORF(seq_str,frame):
	'''
	Returns (start,<length of longest orf string>)
	Returns (0,0) if no ORF found
	'''
	orf = 0
	start_codons = ['ATG']
	stop_codons = ['TAA','TAG','TGA']
	start = -1


	for i in range(frame-1,len(seq_str),3):
		if(seq_str[i:i+3] in start_codons):
			for j in range(i+3,len(seq_str),3):
				if(seq_str[j:j+3] in stop_codons):
					if orf < (j + 3 - i):
						orf = j + 3 - i
						start = i+1
					break

	return (start, orf)

def getRepeats(seq_dict,n):
	'''
	Returns {<sequence_string_length_n>:<occurence_number>}
	'''
	repeat_dict = {}
	for seq_id in seq_dict:
		sequence = seq_dict[seq_id]
		repeat_dict_local = {}
		for i in range(len(sequence)-n+1):
			subseq = sequence[i:i+n]
			if(subseq not in repeat_dict_local):
				repeat_dict_local[subseq] = len(re.findall(subseq, sequence, overlapped=True))
		repeat_dict = dict(Counter(repeat_dict) + Counter(repeat_dict_local))		
	return repeat_dict


seq_dict = {}
curr_id = ""

try:
	multi_fasta = open(file_name)
	for line in multi_fasta:
		if(line.startswith('>')):
			curr_id = line.split(' ')[0].replace('>','').replace('\n','')
			seq_dict[curr_id] = ""
		else:
			seq_dict[curr_id] += line.replace('\n','')
	
	'''
	print("The file has {} records".format(len(seq_dict)))
	
	seq_dict_sorted = getSortedTuplesFromDict(seq_dict)
	print("The sequence with maximum length {} has id {}".\
		format(len(seq_dict_sorted[-1][1]),seq_dict_sorted[-1][0]))
	
	print("The sequence with minimum length {} has id {}".\
		format(len(seq_dict_sorted[0][1]),seq_dict_sorted[0][0]))

	if(len(seq_dict_sorted[-1][1]) == len(seq_dict_sorted[-2][1])):
		print('There\'s more than one longest sequence')
	if(len(seq_dict_sorted[0][1]) == len(seq_dict_sorted[1][1])):
		print('There\'s more than one shortest sequence')

	seq_id = 'gi|142022655|gb|EQ086233.1|16'
	orfs_frame1 = getOpenReadingFrames(seq_dict,1)
	orfs_frame1_val = orfs_frame1[seq_id][1]
	orfs_frame2 = getOpenReadingFrames(seq_dict,2)
	orfs_frame2_val = orfs_frame2[seq_id][1]
	orfs_frame3 = getOpenReadingFrames(seq_dict,3)
	orfs_frame3_val = orfs_frame3[seq_id][1]

	print(max(orfs_frame1_val, \
			orfs_frame2_val, \
			orfs_frame3_val))
	'''
	'''
	dict_repeats = getRepeats(seq_dict,6)
	print(sorted(dict_repeats.items(),key=lambda x:x[1])[-1][1])
	
	dict_repeats = getRepeats(seq_dict,12)
	dict_repeats_sorted = sorted(dict_repeats.items(),key=lambda x:x[1],reverse=True)
	max_occ = dict_repeats_sorted[0][1]
	for i in range(len(dict_repeats_sorted)):
		if(dict_repeats_sorted[i][1] != max_occ):
			break
	print('Max occurs {} times'.format(i))
	'''
	dict_repeats = getRepeats(seq_dict,7)
	print(sorted(dict_repeats.items(),key=lambda x:x[1])[-1][0])



except FileNotFoundError as e:
	print("File not found!")

