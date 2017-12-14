import numpy as np

# defacult score matrix [match, mismatch, insertion, deletion]
SCORE_MATRIX = [5, -2, -6, -6]


def linear_align(seq_1, seq_2, matrix=SCORE_MATRIX):
    m = len(seq_1)
    n = len(seq_2)
    # empty array is  not really empty!
    scores = np.empty([2, n + 1])

    # initilization
    scores[0][0] = 0
    for i in range(1, n + 1):
        scores[0][i] = scores[0][i - 1] + matrix[2]

    # start to fill the table, we always fill the second row based on the first row result
    for i in range(1, m + 1):
        scores[1][0] = scores[0][0] + matrix[3]
        for j in range(1, n + 1):
            if seq_1[i - 1] == seq_2[j - 1]:
                # match
                scores[1][j] = max(scores[0][j - 1] + matrix[0],  # match +5
                                   scores[1][j - 1] + matrix[2],  # insertion -6
                                   scores[0][j] + matrix[3],  # deletion -6
                                   )
            else:
                scores[1][j] = max(scores[0][j - 1] + matrix[1],  # mismatch -2
                                   scores[1][j - 1] + matrix[2],  # insertion -6
                                   scores[0][j] + matrix[3],  # deletion -6
                                   )

        # finish one row, let first row equal to second row
        scores[0] = scores[1]

    return scores[1][n]


def linear_reverse(seq_1, seq_2, matrix=SCORE_MATRIX):
    m = len(seq_1)
    n = len(seq_2)
    # empty array is  not really empty!
    scores = np.empty([2, n + 1])

    # initilization
    scores[1][n] = 0
    for i in range(n - 1, -1, -1):
        scores[1][i] = scores[1][i + 1] + matrix[2]

    # start to fill the table, we always fill the second row based on the first row result
    for i in range(m - 1, -1, -1):
        scores[0][n] = scores[1][n] + matrix[3]
        for j in range(n - 1, -1, -1):
            if seq_1[i] == seq_2[j]:
                # match
                scores[0][j] = max(scores[1][j + 1] + matrix[0],  # match +5
                                   scores[0][j + 1] + matrix[2],  # insertion -6
                                   scores[1][j] + matrix[3],  # deletion -6
                                   )
            else:
                scores[0][j] = max(scores[1][j + 1] + matrix[1],  # mismatch -2
                                   scores[0][j + 1] + matrix[2],  # insertion -6
                                   scores[1][j] + matrix[3],  # deletion -6
                                   )

        # finish one row, let first row equal to second row
        scores[1] = scores[0]

    return scores[0][0]
