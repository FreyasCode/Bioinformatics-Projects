def smith_waterman(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    # Initialize the scoring matrix
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    rows, cols = len(seq1) + 1, len(seq2) + 1
    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    subtract = 0
    for i in range(rows):
        score_matrix[i][0] -= subtract
        subtract += 1

    # Fill the scoring matrix
    max_score = 0
    max_position = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = score_matrix[i - 1][j - 1] + match_score
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch_score

            top_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty

            score_matrix[i][j] = max(0, diagonal_score, top_score, left_score)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_position = (i, j)
    # print(score_matrix[len_seq1][len_seq2])

    # Traceback to find the alignment
    max_position = score_matrix[len_seq1].index(max(score_matrix[len_seq1]))

    align1 = ""
    align2 = ""
    i, j = len_seq1, max_position
    # print(score_matrix)

    while i > 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i - 1][j - 1] if i > 0 and j > 0 else float('-inf')
        up_score = score_matrix[i][j - 1] if j > 0 else float('-inf')
        left_score = score_matrix[i - 1][j] if i > 0 else float('-inf')

        if current_score == diagonal_score + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif current_score == left_score + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        elif current_score == up_score + gap_penalty:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2, score_matrix[len_seq1][max_position]

# Example usage:
sequence1 = "ACCGTT"
sequence2 = "ACGTAACCTTT"
alignment_result = smith_waterman(sequence1, sequence2)

print(alignment_result[2])
print(alignment_result[0])
print(alignment_result[1])
