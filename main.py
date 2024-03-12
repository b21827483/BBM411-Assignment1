#Buğrahan Halıcı
import numpy as np
import sys

def calculate_alignment_score(alignedA, alignedB):

    middle = ""
    alignment_score = 0
    for i in range(len(alignedA)):
        if (alignedA[i] == alignedB[i] and alignedA[i] != "-" and alignedB[i] != "-"):
            alignment_score += scoring_matrix[alignedA[i]][alignedB[i]]
            middle = middle + "|"

        elif(alignedA[i] != "-" and alignedB[i] != "-"):
            alignment_score += scoring_matrix[alignedA[i]][alignedB[i]]
            middle += " "

        elif(alignedA[i] == "-" or alignedB[i] == "-"):
            if i == 0:
                alignment_score += gap_penalty
                middle += " "

            elif alignedA[i] == "-" and alignedA[i-1] == "-":
                alignment_score += gap_extend
                middle += " "

            elif alignedB[i] == "-" and alignedB[i-1] == "-":
                alignment_score += gap_extend
                middle += " "

            else:
                alignment_score += gap_penalty
                middle += " "

    return middle, alignment_score

def calculate_identity_value(middle):
    identity_value = 0.0
    for i in range(len(middle)):
        if middle[i] == "|":
            identity_value += 1

    identity_value /= len(middle)

    return identity_value

def fill_matrix():
    # insert gap penalties to the main matrix
    matrix[0][0] = float("-inf")
    Ix_matrix[0][0] = gap_penalty
    Iy_matrix[0][0] = gap_penalty

    for i in range(1, len(proteinB) + 1):
        matrix[i][0] = float("-inf")
        Ix_matrix[i][0] = float("-inf")
        Iy_matrix[i][0] = gap_penalty + i * gap_extend

    for i in range(1, len(proteinA) + 1):
        matrix[0][i] = float("-inf")
        Ix_matrix[0][i] = gap_penalty + i * gap_extend
        Iy_matrix[0][i] = float("-inf")


def global_alignment():
    for i in range(1, len(proteinB) + 1):
        for j in range(1, len(proteinA) + 1):
            score = scoring_matrix[proteinB[i-1]][proteinA[j-1]]
            matrix[i][j] = score + max(matrix[i - 1][j - 1], Ix_matrix[i - 1][j - 1], Iy_matrix[i - 1][j - 1])
            Ix_matrix[i][j] = max(matrix[i][j-1] + gap_penalty, Ix_matrix[i][j-1] + gap_extend, Iy_matrix[i][j-1] + gap_penalty)
            Iy_matrix[i][j] = max(matrix[i-1][j] + gap_penalty, Ix_matrix[i-1][j] + gap_penalty , Iy_matrix[i-1][j] + gap_extend)

    return max(matrix[len(proteinB)][len(proteinA)], Ix_matrix[len(proteinB)][len(proteinA)], Iy_matrix[len(proteinB)][len(proteinA)])


def local_alignmet():
    highest = 0
    location = (1, 1)

    for i in range(1, len(proteinB) + 1):
        for j in range(1, len(proteinA) + 1):
            score = scoring_matrix[proteinB[i-1]][proteinA[j-1]]
            matrix[i][j] = max(0, matrix[i - 1][j - 1] + score, Ix_matrix[i - 1][j - 1] + score, Iy_matrix[i - 1][j - 1] + score)
            Ix_matrix[i][j] = max(0, matrix[i][j-1] + gap_penalty, Ix_matrix[i][j-1] + gap_extend, Iy_matrix[i][j-1] + gap_penalty)
            Iy_matrix[i][j] = max(0, matrix[i-1][j] + gap_penalty, Ix_matrix[i-1][j] + gap_penalty , Iy_matrix[i-1][j] + gap_extend)

            if (max(matrix[i][j], Ix_matrix[i][j], Iy_matrix[i][j]) >= highest):
                highest = max(matrix[i][j], Ix_matrix[i][j], Iy_matrix[i][j])
                location = (i, j)

    return highest, location

def read_sequences():
    sequences_txt = sys.argv[1]
    in_file = open(sequences_txt, "r")
    lines = in_file.readlines()
    proteinA = ""
    proteinB = ""
    count = 0

    for line in lines:
        line = ''.join(i for i in line if i !='\n') #removes escape character
        if count == 0:
            count += 1
            proteinA = line
        else:
            proteinB = line
            count = 0

    in_file.close()
    return proteinA, proteinB

def create_scoring_matrix():
    count = 0
    in_file2 = sys.argv[3]
    scoring_matrix_txt = open(in_file2, "r")
    lines = scoring_matrix_txt.readlines()

    scoring_matrix = {}
    proteins = []

    for line in lines:
        line = ''.join(i for i in line if i != '\n')  # removes escape character
        scores = line.split()

        if count == 0:
            proteins = scores
            scoring_matrix = dict.fromkeys(scores)
            for i in range(len(scores)):
                temp = dict.fromkeys(scores)
                scoring_matrix[scores[i]] = temp
            count += 1
            continue

        for i in range(len(proteins)):
            scoring_matrix[scores[0]][proteins[i]] = int(scores[i+1])

    return scoring_matrix


def traceback(isGlobal):
    alignedA = ""
    alignedB = ""

    A_len = len(proteinA)
    B_len = len(proteinB)

    if isGlobal:
        pointer = global_alignment()

        if pointer == matrix[B_len][A_len]:
            pointer = matrix
        elif pointer == Ix_matrix[B_len][A_len]:
            pointer = Ix_matrix
        else:
            pointer = Iy_matrix

    if not isGlobal:

        highest_point, pointer = local_alignmet()
        A_len = pointer[1]
        B_len = pointer[0]
        if highest_point == matrix[B_len][A_len]:
            pointer = matrix
        elif highest_point == Ix_matrix[B_len][A_len]:
            pointer = Ix_matrix
        else:
            pointer = Iy_matrix

    while A_len > 0 or B_len > 0 :
        if np.array_equal(pointer, matrix):
            alignedA = proteinA[A_len - 1] + alignedA
            alignedB = proteinB[B_len - 1] + alignedB

            if pointer[B_len][A_len] == scoring_matrix[proteinB[B_len - 1]][proteinA[A_len - 1]] + matrix[B_len - 1][A_len - 1]:
                A_len -= 1
                B_len -= 1

            elif pointer[B_len][A_len] == scoring_matrix[proteinB[B_len - 1]][proteinA[A_len - 1]] + Ix_matrix[B_len - 1][A_len - 1]:
                pointer = Ix_matrix
                A_len -= 1
                B_len -= 1
            elif pointer[B_len][A_len] == scoring_matrix[proteinB[B_len - 1]][proteinA[A_len - 1]] + Iy_matrix[B_len - 1][A_len - 1]:
                pointer = Iy_matrix
                A_len -= 1
                B_len -= 1
            else:
                alignedA = alignedA[1:]
                alignedB = alignedB[1:]
                break

        elif np.array_equal(pointer, Ix_matrix):
            alignedA = proteinA[A_len - 1] + alignedA
            alignedB = "-" + alignedB

            if A_len == 0 and B_len > 1:
                break

            elif pointer[B_len][A_len] == matrix[B_len][A_len - 1] + gap_penalty:
                pointer = matrix
                A_len -= 1

            elif pointer[B_len][A_len] == Ix_matrix[B_len][A_len-1] + gap_extend:
                A_len -= 1

            elif pointer[B_len][A_len] == Iy_matrix[B_len][A_len-1] + gap_penalty:
                pointer = Iy_matrix
                A_len -= 1
            else:

                alignedA = alignedA[1:]
                alignedB = alignedB[1:]
                break

        elif np.array_equal(pointer, Iy_matrix):
            alignedA = "-" + alignedA
            alignedB = proteinB[B_len - 1] + alignedB
            if A_len == 0 and B_len == 1:
                break

            elif pointer[B_len][A_len] == matrix[B_len - 1][A_len] + gap_penalty:
                pointer = matrix
                B_len -= 1

            elif pointer[B_len][A_len] == Ix_matrix[B_len-1][A_len] + gap_penalty:
                B_len -= 1

            elif pointer[B_len][A_len] == Iy_matrix[B_len-1][A_len] + gap_extend:
                B_len -= 1
            else:

                alignedA = alignedA[1:]
                alignedB = alignedB[1:]
                break

    return alignedA, alignedB

proteinA, proteinB = read_sequences()
scoring_matrix = create_scoring_matrix()

#initialize matrices
matrix = np.zeros((len(proteinB) + 1, len(proteinA) + 1))
Ix_matrix = np.zeros((len(proteinB) + 1, len(proteinA) + 1))
Iy_matrix = np.zeros((len(proteinB) + 1, len(proteinA) + 1))

algorithm = sys.argv[2]
gap_penalty = int(sys.argv[4])
gap_extend = int(sys.argv[5])

def globalAl():
    # fill the matrix using Needleman_Wunsch algorithm
    fill_matrix()
    alignedA, alignedB = traceback(True)
    middle, alignment_Score = calculate_alignment_score(alignedA, alignedB)

    print()
    print("GLOBAL ALIGNMENT")
    print("Alignment Score:", alignment_Score)
    print("Identity Value:" , round(calculate_identity_value(middle) * 100, 1) , "%")

    print(alignedA)
    print(middle)
    print(alignedB)

    print()
    print("**********************************")
    print()

def localAl():
    #reassign matrix values for local alignment
    matrix = np.zeros((len(proteinB) + 1, len(proteinA) + 1))
    Ix_matrix = np.zeros((len(proteinB) + 1, len(proteinA) + 1))
    Iy_matrix = np.zeros((len(proteinB) + 1, len(proteinA) + 1))

    alignedA, alignedB = traceback(False)
    middle, alignment_Score = calculate_alignment_score(alignedA, alignedB)

    print("LOCAL ALIGNMENT")
    print("Alignment Score:", alignment_Score)
    print("Identity Value:" , round(calculate_identity_value(middle) * 100, 1) , "%")

    print(alignedA)
    print(middle)
    print(alignedB)

if algorithm == "global":
    globalAl()
if algorithm == "local":
    localAl()

