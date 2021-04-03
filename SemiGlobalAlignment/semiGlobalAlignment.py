# Since edit distance is the same as alignment distance
# it just makes sense to compute the global alignment
# distance

import sys

def getSequencesFromFastaFormat(data):
    _, E1,E2 = data.split(">")
    S1 = "".join(E1.strip().split("\n")[1:])
    S2 = "".join(E2.strip().split("\n")[1:])
    return (S1,S2)

def getOptimalHorizontalGaps( Mat, i, j, S1, S2, gap):
    if(j == len(S2) or j == 0):
        return max( [ Mat[i][j-k] for k in range(1,j+1) ] )
    return max( [ Mat[i][j-k] + gap(k) for k in range(1,j+1) ] )

def getOptimalVerticalGaps( Mat, i, j, S1, S2, gap):
    if(i == len(S1) or i == 0):
        return max( [ Mat[i-k][j] for k in range(1,i+1) ] )
    return max( [ Mat[i-k][j] + gap(k) for k in range(1,i+1) ] )
    
def computeSimilarityMatrix(S1, S2, Mat, gap = lambda g: g):
    m = len(S1)
    n = len(S2)
    
    for i in range(1, m+1):
        for j in range(1,n+1):
            matchCost = -1
            if(S1[i-1] == S2[j-1]):
                matchCost = 1

            verticalGapCost = gap(1)
            HorizontalGapCost = gap(1)
            if(i == m or i == 0):
                HorizontalGapCost = 0
            if(j == n or j == 0):
                verticalGapCost = 0

            Mat[i][j] = max([
                              Mat[i][j-1] + HorizontalGapCost,
                              Mat[i-1][j] + verticalGapCost,
                              Mat[i-1][j-1] + matchCost
                            ])

def getSimilarityScore(S1, S2, EditMatrix):
    m = len(S1)
    n = len(S2)
    return EditMatrix[m][n]
    
#Linear gap cost by default.
def initializeSimilarityMatrix(S1, S2, gap = lambda g: g ):
    m = len(S1)
    n = len(S2)
    
    SimilarityMatrix = []

    for i in range(m+1):
        SimilarityMatrix.append( [0] * (n+1) )
    
    return SimilarityMatrix

def getOptimalAlignment(S1_, S2_, EditMatrix, gap = lambda g: -g):
    Align_S1 = ""
    Align_S2 = ""
    i = len(S1_)
    j = len(S2_)
    S1 = "-" + S1_
    S2 = "-" + S2_
    while( i != 0 or j != 0 ):
        verticalGapCost = gap(1)
        HorizontalGapCost = gap(1)
        if(i == len(S1_) or i == 0):
            HorizontalGapCost = 0
        if(j == len(S2_) or j == 0):
            verticalGapCost = 0
        Elem = EditMatrix[i][j]
        if  ( i != 0 and Elem == (EditMatrix[i-1][j] + verticalGapCost) ):
            Align_S1 = S1[i] + Align_S1
            Align_S2 = "-" + Align_S2
            i = i - 1
        elif( j != 0 and Elem == (EditMatrix[i][j-1] + HorizontalGapCost) ):
            Align_S1 = "-" + Align_S1
            Align_S2 = S2[j] + Align_S2
            j = j - 1
        else:
            matchCost = -1
            if(S1[i] == S2[j]):
                matchCost = 1
            if(Elem == (EditMatrix[i-1][j-1] + matchCost) ):
                Align_S1 = S1[i] + Align_S1
                Align_S2 = S2[j] + Align_S2
                i = i - 1
                j = j - 1
            
    return Align_S1, Align_S2

if( __name__ == "__main__"):
    if(len(sys.argv) < 2 ):
        print("usage: python3 editDistanceAlignment.py <data_set_file_name>")
        sys.exit(-1)

    dataset_filename = sys.argv[1]
    data = ""

    with open(dataset_filename,"r") as f:
        data = f.read()

    S1, S2 = getSequencesFromFastaFormat(data)
    
    SimilarityMatrix = initializeSimilarityMatrix(S1,S2, gap = lambda g: -g)
    computeSimilarityMatrix(S1, S2, SimilarityMatrix, gap = lambda g: -g)

    SimilarityScore = getSimilarityScore(S1, S2, SimilarityMatrix)
    Align_S1, Align_S2 = getOptimalAlignment(S1, S2, SimilarityMatrix, gap = lambda g: -g)

    
    print(SimilarityScore)
    print(Align_S1)
    print(Align_S2)
