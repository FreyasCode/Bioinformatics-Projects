import sys
from Bio import SeqIO

# same as part 1, align two sequences
def align(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    # Initialize the DP scoring matrix
    len_seq1 = len(seq1) - 1
    len_seq2 = len(seq2)

    rows, cols = len(seq1), len(seq2) + 1
    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    subtract = 0
    for i in range(rows):
        subtract += 1
        score_matrix[i][0] -= subtract
    for j in range(1, cols): # 1 for match and -1 for mismatch
        score_matrix[0][j] = 1 if seq1[0] == seq2[j - 1] else -1

    # Fill the scoring matrix
    max_score = 0
    max_position = (0, 0)

    for i in range(1, rows): 
        for j in range(1, cols):
            if seq1[i] == seq2[j - 1]: # 1 for match and -1 for mismatch
                diagonal_score = score_matrix[i - 1][j - 1] + match_score
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch_score

            top_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty

            score_matrix[i][j] = max(diagonal_score, top_score, left_score)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_position = (i, j)

    # Traceback to find the alignment
    seq2 = "-" + seq2
    max_position = score_matrix[len_seq1].index(max(score_matrix[len_seq1]))

    align1 = ""
    align2 = ""
    i, j = len_seq1, max_position
    

    while i >= 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i - 1][j - 1] if i > 0 and j > 0 else float('-inf')
        up_score = score_matrix[i][j - 1] if j > 0 else float('-inf')
        left_score = score_matrix[i - 1][j] if i > 0 else float('-inf')

        if current_score == left_score + gap_penalty: # if the current score is equal to the left score plus the gap penalty
            align1 = seq1[i] + align1
            align2 = "-" + align2
            i -= 1
        elif current_score == up_score + gap_penalty: # if the current score is equal to the up score plus the gap penalty
            align1 = "-" + align1
            align2 = seq2[j] + align2
            j -= 1
        elif current_score == diagonal_score + (match_score if seq1[i] == seq2[j] else mismatch_score): # if the current score is equal to the diagonal score plus the match score if the current character in seq1 is equal to the current character in seq2, else the mismatch score
            align1 = seq1[i] + align1
            align2 = seq2[j] + align2
            i -= 1
            j -= 1
        else:
            align1 = seq1[i] + align1
            align2 = seq2[j] + align2
            break

    return str(score_matrix[len_seq1][max_position]), align1, align2

# function to parse input 
def parse(file_name):
    content = []
    with open(file_name) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            content.append({'name': record.description, 'sequence': str(record.seq)})
    return content

# Running the program
query = parse(sys.argv[2])
sequences = parse(sys.argv[4])

max_score = float('-inf')
max_alignment = []

# Looping through sequences in the library and comparing with the query seq to find the best alignment
for seq in sequences:
    score, align1, align2 = align(query[0]['sequence'], seq['sequence'])
    score = int(score)

    if score > max_score:
        max_score = score
        max_alignment = [">" + seq['name'], align1, align2]

# write output to file 
with open(sys.argv[6], 'w') as file:
    file.write(max_alignment[0] + "\n")
    file.write(str(max_score) + "\n")
    file.write(max_alignment[1] + "\n")
    file.write(max_alignment[2] + "\n")