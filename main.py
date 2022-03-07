import numpy as np
#print("alaa")

def valid_seq(seq):
    seq = seq.upper()
    valid = 'ACTG'
    return(all(i in valid for i in seq))

def maximum(a, b, c):
    if (a >= b) and (a >= c):
        largest = a
    elif (b >= a) and (b >= c):
        largest = b
    else:
        largest = c
    return largest

def initializing_m(seq1,seq2):
    m=len(seq1)
    n=len(seq2)
    global matrix
    matrix = [[0 for x in range(m+1)] for y in range(n+1)] #setting zeros to our matrix
    print(matrix)
    #return matrix

def filling_m(seq1,seq2,gap):
    print(matrix)
    max= []
    m=len(seq1)
    n=len(seq2)
    for i in range(1,n+1):
        for j in range(1,m+1):
            if seq1[j-1] == seq2[i-1]: #specify the match and mismatch of the diagonal
                score= 5
            else:
                score=-3
            #assigning the intersaction cell to the max num the came from the diagonal or top or left
            #print(i,j)
            cell= maximum(matrix[i-1][j]+gap, matrix[i][j-1]+gap,matrix[i-1][j-1]+score)
            #print(cell)
            if cell<1:
               cell=0
            matrix[i][j]=cell
            #print(matrix[i][j])
            if len(max) == 0 or max[2]==cell: #putting the optimal cells
                max.append(i)
                max.append(j)
                max.append(cell)
            else:
                if max[2]<cell:
                    max*=0
                    max.append(i)
                    max.append(j)
                    max.append(cell)
    #print(matrix)
    #print(max)
    return (max)

def tracking_m (seq1,seq2,max,gap):
    global TBseq1
    global TBseq2
    TB1=[]
    TB2=[]
    for a in range(0, len(max), 3):
       # print(a)
        TBseq1 = ' '
        TBseq2 = ' '
        i = max[a]
        j = max[a + 1]
        while matrix[i][j] != 0:
            left = matrix[i][j - 1] + gap
            up = matrix[i - 1][j] + gap
            if seq1[j - 1] == seq2[i - 1]:
                score = 5
            else:
                score = -3
            diagonal = matrix[i - 1][j - 1] + score
            if matrix[i][j] == diagonal:
                TBseq1 += seq1[j - 1]
                TBseq2 += seq2[i - 1]
                i -= 1
                j -= 1
            elif matrix[i][j] == up:
                TBseq1 += '-'
                TBseq2 += seq2[i - 1]
                i -= 1

            elif matrix[i][j] == left:
                TBseq2 += '-'
                TBseq1 += seq1[j - 1]
                j -= 1
        TBseq1 = TBseq1[::-1]
        TBseq2 = TBseq2[::-1]
        TB1.append(TBseq1)
        TB2.append(TBseq2)


    print(TB1)
    print(TB2)

def local_alignment():
    seq1 = input("Enter the first DNA sequence: ")
    seq2 = input("Enter the second DNA sequence: ")
   # gap = input("Enter gap panilty: ")
    #match= input("Enter match value: ")
    #mismatch = input("Enter mismatch value: ")

    if (valid_seq(seq1)) and (valid_seq(seq2)):
        initializing_m(seq1, seq2)
        gap=-4
        max = filling_m(seq1, seq2,gap)
        tracking_m(seq1, seq2, max, gap)
    else:
        print("invalid DNA sequence")
local_alignment()









