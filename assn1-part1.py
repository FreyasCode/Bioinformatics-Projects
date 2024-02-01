def align(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    # Initialize the scoring matrix
    len_seq1 = len(seq1) - 1
    len_seq2 = len(seq2)

    rows, cols = len(seq1), len(seq2) + 1
    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    subtract = 0
    for i in range(rows):
        subtract += 1
        score_matrix[i][0] -= subtract
        
    
    for j in range(1, cols):
        score_matrix[0][j] = 1 if seq1[0] == seq2[j - 1] else -1

    # Fill the scoring matrix
    max_score = 0
    max_position = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            # print(str(i) + "," + str(j) + ": comparing " + seq1[i] + " and " + seq2[j - 1])
            if seq1[i] == seq2[j - 1]:
                diagonal_score = score_matrix[i - 1][j - 1] + match_score
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch_score

            top_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty

            # print([score_matrix[i - 1][j - 1], score_matrix[i - 1][j], score_matrix[i][j - 1]])
            # print([diagonal_score, top_score, left_score])

            score_matrix[i][j] = max(diagonal_score, top_score, left_score)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_position = (i, j)
    # print(score_matrix[len_seq1][len_seq2])

    # Traceback to find the alignment
    seq2 = "-" + seq2
    max_position = score_matrix[len_seq1].index(max(score_matrix[len_seq1]))

    align1 = ""
    align2 = ""
    i, j = len_seq1, max_position
    print(score_matrix)
    

    while i >= 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i - 1][j - 1] if i > 0 and j > 0 else float('-inf')
        up_score = score_matrix[i][j - 1] if j > 0 else float('-inf')
        left_score = score_matrix[i - 1][j] if i > 0 else float('-inf')
        print([i,j] + [seq1[i], seq2[j]])

        if current_score == diagonal_score + (match_score if seq1[i] == seq2[j] else mismatch_score):
            print("diagonal")
            align1 = seq1[i] + align1
            align2 = seq2[j] + align2
            i -= 1
            j -= 1
        elif current_score == left_score + gap_penalty:
            print("left")
            align1 = seq1[i] + align1
            align2 = "-" + align2
            i -= 1
        elif current_score == up_score + gap_penalty:
            print("up")
            align1 = "-" + align1
            align2 = seq2[j] + align2
            j -= 1
        else:
            align1 = seq1[i] + align1
            align2 = seq2[j] + align2
            break

    return align1, align2, score_matrix[len_seq1][max_position]

# Example usage:
# sequence1 = "AAAA"
# sequence2 = "GGGGGGGGGGGG"
# sequence1 = "ACT"
# sequence2 = "CT"
# sequence1 = "AACCCCTAG"
# sequence2 = "TTAATCCCCAGGGTCGTTT"
# sequence1 = "ACCGTT"
# sequence2 = "ACGTAACCTTT"
# sequence1 = "ACCGTTACCGTT"
# sequence2 = "ACGTAACCTTT"
# sequence1 = "ACCGTTACCGTTACCGTT"
# sequence2 = "ACGTAACCTTT"
# sequence1 = "AAAAAAAACT"
# sequence2 = "CT"
# sequence1 = "ATGGGGGGG"
# sequence2 = "AT"
sequence1 = "ATGGGGGGG"
sequence2 = "TG"
alignment_result = align(sequence1, sequence2)

print(alignment_result[2])
print(alignment_result[0])
print(alignment_result[1])
