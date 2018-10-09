import sys # to read arguments from command line

def print_dp_matrix(pattern, text, dp_matrix):
  # print the dynamic programming matrix with row and cols titles
  header = list(" -"+text)
  rows = list("-"+pattern)
  print ("The dynamic programming table is:")
  #pretty print matrix from:
  #https://stackoverflow.com/questions/17870612/printing-a-two-dimensional-array-in-python
  print ('\n'.join([' '.join(['{:4}'.format(item) for item in header])]))
  for i in range(len(dp_matrix)):
    print rows[i],
    #print (rows[i], end = "") #python 3
    for val in dp_matrix[i]:
      print '{:4}'.format(val),
      #print ('{:4}'.format(val), end="") python 3
    print
  print

def backtrace_matrix(pattern, text, dp_matrix):
  # choose pathway in dynamic matrix and generate CIGAR (MDIMM...)
  i = len(pattern)
  j = len(text)
  CIGAR = []
  while i>0 and j>0:
    if dp_matrix[i][j] == dp_matrix[i-1][j] - 2: #Deletion
      i -= 1
      CIGAR.insert(0, "D")
    elif dp_matrix[i][j] == dp_matrix[i][j-1] - 4: #Insertion
      j -= 1
      CIGAR.insert(0,'I')
    else: #Substitution
      i -= 1
      j -= 1
      if pattern[i] == text[j]:
        CIGAR.insert(0, "M")
      else:
        CIGAR.insert(0, "X")
  if i > 0:
    for _ in range(i): CIGAR.insert(0, "D")
  if j > 0:
    for _ in range(j): CIGAR.insert(0, "I")
  return CIGAR

def score_match(pair, substitution_matrix):
  # read substitution matrix when edit_distance_dp function gets match-missmatch
  # this version only works with blosum45, blosum62 and blosum80
  from Bio.SubsMat import MatrixInfo # get matrix using biopython
  # choose matrix depending on passed argument
  if substitution_matrix == "blosum62":
    matrix = MatrixInfo.blosum62
  elif substitution_matrix == "blosum45":
    matrix = MatrixInfo.blosum45
  else:
    matrix = MatrixInfo.blosum80
  #read value from matrix
  if pair in matrix:
    return matrix[pair]
  else:
    return matrix[tuple(reversed(pair))]

def edit_distance_dp(pattern,text, substitution_matrix):
  # Init dynamic programming matrix
  dp_matrix = [[0 for i in range(len(text)+1)] for j in range(len(pattern)+1)]
  for i in range(len(pattern)+1):
    dp_matrix[i][0] = -2*i
  for j in range(len(text)+1):
    dp_matrix[0][j] = -4*j
  # Compute cells
  for i in range(1,len(pattern)+1):
    for j in range(1, len(text)+1):
      dp_matrix[i][j] = max(
        dp_matrix[i-1][j-1] + (score_match((pattern[i-1], text[j-1]), substitution_matrix)),
        dp_matrix[i][j-1] - 4, #insertion
        dp_matrix[i-1][j] - 2) #deletion
  # Print matrix, comment for big alignments
  #print_dp_matrix(pattern, text, dp_matrix)
  # Generate score and alignment
  score = dp_matrix[i][j]
  CIGAR = backtrace_matrix(pattern, text, dp_matrix)
  return (score, CIGAR)

def pretty_alignment(pattern,text,substitution_matrix):
  # add gaps to sequences to generate alignment and create line with | marking identities
  (score, CIGAR) = edit_distance_dp(pattern,text,substitution_matrix)
  line1 = ""
  line2 = ""
  line3 = ""
  pattern_index = 0
  text_index = 0
  for i in CIGAR:
    if i == "M":
      line1 += pattern[pattern_index]
      line2 += "|"
      line3 += text[text_index]
      pattern_index += 1
      text_index += 1
    if i == "I":
      line1 += "-"
      line2 += " "
      line3 += text[text_index]
      text_index += 1
    if i == "D":
      line1 += pattern[pattern_index]
      line2 += " "
      line3 += "-"
      pattern_index += 1
    if i == "X":
      line1 += pattern[pattern_index]
      line2 += " "
      line3 += text[text_index]
      pattern_index += 1
      text_index += 1
  return line1, line2, line3, score, CIGAR #CIGAR is returned only for testing

def read_fasta(file):
  # read fasta file and generate list of lists with name and seq
  list = []
  f = open(file, 'r')
  header = f.readline()
  while len(header) != 0:
    seq = []
    seq.append(header[1:-1])
    r = ""
    while True:
      s = f.readline()
      if len(s) == 0 or s[0] == ">":
        header = s
        break
      r = r + s [:-1]
    seq.append(r)
    list.append(seq)
  return list

def print_alignment(fasta_file, substitution_matrix, block_size):
  # get first 2 seq names and seqs from input file
  # generate a "paper-like" alignment with the read sequences
  pattern,text = read_fasta(fasta_file)[0][1], read_fasta(fasta_file)[1][1]
  pattern_name, text_name = read_fasta(fasta_file)[0][0], read_fasta(fasta_file)[1][0]
  line1, line2, line3, score, CIGAR = pretty_alignment(pattern,text,substitution_matrix) #alignment only for testing
  if block_size == "inf" or block_size == "all" or block_size == 0:
    return line1+"\n"+line2+"\n"+line3+"\n"
  else:
    aa1_len = 0
    aa2_len = 0
    result = ""
    for index in range(0, len(line1), block_size):
      aa1 = line1[index : index + block_size]
      aa2 = line3[index : index + block_size]
      aa1_len += len(aa1.replace("-",""))
      aa2_len += len(aa2.replace("-",""))
      result = result + "\n" +\
      pattern_name + "\t" + aa1 + "  " + str(aa1_len) + "\n" +\
      "\t" + line2[index : index + block_size] + "\n" +\
      text_name + "\t" + aa2 + "  " + str(aa2_len) + "\n"
      if (index + block_size) > len(line1):
        break
    return score, result, CIGAR #CIGAR is returned only for testing

def nw_protein(fasta_file, substitution_matrix = "blosum62", block_size = 70):
  # print alignment of the first 2 sequences in input file.fasta and the
  # score of the alignment
  # substitution_matrix and block_size are optional, default are blosum62 and 70
  (score, alignment, CIGAR) = print_alignment(fasta_file, substitution_matrix, block_size) #CIGAR only for testing. Comment next line unless testing
  #print CIGAR
  return "\nAlignment score is: " + str(score) + "\n" + alignment





##################################
# reading command line arguments #
##################################

#init dictionary of arguments
my_dict = {'--input': 0, '--subs_mat': 0, '--block_size': 0}

#add arguments to dictionary
if len(sys.argv) < 3: #exit script if not enough passed arguments
  sys.exit("syntax is caseofuse7-4.py --input filename [--subs_mat matrixname --block_size number]")
else:
  for i in range(1,len(sys.argv)):
    if sys.argv[i] in my_dict.keys():
      my_dict[sys.argv[i]] = sys.argv[i+1]

#print my_dict # this line is for testing

#pass arguments to function
if my_dict['--input'] == 0: #exit script if no filename provided
  sys.exit("syntax is caseofuse7-4.py --input filename [--subs_mat matrixname --block_size number]")
elif my_dict['--subs_mat'] != 0 and my_dict['--block_size'] != 0:
  print nw_protein(my_dict['--input'], my_dict['--subs_mat'], int(my_dict['--block_size']))
elif my_dict['--subs_mat'] != 0 and my_dict['--block_size'] == 0:
  print nw_protein(my_dict['--input'], substitution_matrix = my_dict['--subs_mat'])
elif my_dict['--subs_mat'] == 0 and my_dict['--block_size'] != 0:
  print nw_protein(my_dict['--input'], block_size = int(my_dict['--block_size']))
elif my_dict['--subs_mat'] == 0 and my_dict['--block_size'] == 0:
  print nw_protein(my_dict['--input'])
else: #exit script if the previous conditions do not hold
  sys.exit("syntax is caseofuse7-4.py --input filename [--subs_mat matrixname --block_size number]")

# py -2.7 align.py --input GHRs.fasta --subs_mat blosum62 --block_size 90
