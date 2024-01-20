def init_scores(m, n, scoring=(1, -1, -1)):
    '''Initialise the score matrix'''
    '''m = rows, n = columns, scoring = scoring scheme used -> tuple(match_score, mismatch_score, gap_penalty)'''
    mat = []

    # we make every element 0
    mat = [[0] * n for _ in range(m)]

    # the first row and first column correspond to indels from the previous scores
    
    # filling the first column
    for i in range(1, m):
        mat[i][0] = i * scoring[2]

    # filling the first row
    for j in range(1, n): 
        mat[0][j] = j * scoring[2]
    
    
    return mat

def match_score(a, b, scoring=(1, -1, -1)):
    '''Returns match/mismatch score of two characters; a, b are the characters checked, scoring is a tuple with the match and mismatch score'''
    if a == b:
        return scoring[0]
    else:
        return scoring[1]

def seq_align(s1, s2, scoring=(1, -1, -1)):
    '''s1, s2 - sequences to be aligned;
    scoring - scoring scheme of the form tuple(match_score, mismatch_score, gap_penalty)'''
    m = len(s1) + 1 # number of rows of the score matrix
    n = len(s2) + 1 # length of columns of the score matrix
    score_matrix = init_scores(m, n, scoring)

    # start filling the cells
    for i in range(1, m):
        for j in range(1, n):
            t1 = score_matrix[i-1][j-1] + match_score(s1[i-1], s2[j-1], scoring) # diagonal match/mismatch
            t2 = score_matrix[i-1][j] + scoring[2] # indel in s1
            t3 = score_matrix[i][j-1] + scoring[2] # indel in s2
            score_matrix[i][j] = max(t1, t2, t3)
            
    print("After filling the score matrix:\n")
    # for printing neatly
    for i in range(m):
        s = ""
        for j in range(n):
            val = score_matrix[i][j]
            if val > 0:
                s += "+"+str(val)+" | "
            elif val == 0:
                s += " "+str(val)+" | "
            else:
                s += str(val)+" | "
        print(s)

    # now we trace back to find the alignment
    a1 = "" # alignment of the first sequence
    a2 = "" # alignment of the second sequence

    # starting from the bottom left element
    i = m-1
    j = n-1

    while i > 0 and j > 0:
        current = score_matrix[i][j]
        diagonal = score_matrix[i-1][j-1]
        left = score_matrix[i][j-1]
        top = score_matrix[i-1][j]

        if current == (diagonal + match_score(s1[i-1], s2[j-1], scoring)):
            a1 += s1[i-1]
            a2 += s2[j-1]
            i -= 1
            j -= 1
        elif current == (left + scoring[2]):
            a1 += "-"
            a2 += s2[j-1]
            j -= 1
        elif current == (top + scoring[2]):
            a1 += s1[i-1]
            a2 += "-"
            i -= 1
    print("Alignment:")
    # since we started from the bottom and moved towards the top, the strings will be reversed
    print(a1[::-1])
    print(a2[::-1]) 

    return score_matrix[m-1][n-1] # return final score

if __name__ == "__main__":
    l = seq_align("GCATGCGA", "GATTACA", (1, -1, -1)) # using the same scoring scheme proposed by Needleman and Wunsch
    print("\nThe final alignment score is {}.".format(l))

