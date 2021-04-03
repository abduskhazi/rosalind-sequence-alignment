# Since edit distance is the same as alignment distance
# it just makes sense to compute the global alignment
# distance

import sys
import math

"""
PAM250 =
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  2 -2  0  0 -3  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3
C -2 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0
D  0 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4
E  0 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4
F -3 -4 -6 -5  9 -5 -2  1 -5  2  0 -3 -5 -5 -4 -3 -3 -1  0  7
G  1 -3  1  0 -5  5 -2 -3 -2 -4 -3  0  0 -1 -3  1  0 -1 -7 -5
H -1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0
I -1 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1
K -1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4
L -2 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1
M -1 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2
N  0 -4  2  1 -3  0  2 -2  1 -3 -2  2  0  1  0  1  0 -2 -4 -2
P  1 -3 -1 -1 -5  0  0 -2 -1 -3 -2  0  6  0  0  1  0 -1 -6 -5
Q  0 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4
R -2 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4
S  1  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3
T  1 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3
V  0 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2
W -6 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0
Y -3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10
"""
ScoringMatrix = ""

def getScoringMatrix(scoringData):
    ScoringMatrix = {}
    list = scoringData.strip().split("\n")
    AminoAcids = list[0].strip().split()
    for i in range(1,len(list)):
        elems = list[i].strip().split()
        aa = elems[0]
        ScoringMatrix[aa] = {}
        for j in range(1, len(elems)):
            ScoringMatrix[aa][AminoAcids[j-1]] = int(elems[j])
    return ScoringMatrix;
    
def getSequencesFromFastaFormat(data):
    _, E1, E2 = data.split(">")
    S1 = "".join(E1.strip().split("\n")[1:])
    S2 = "".join(E2.strip().split("\n")[1:])
    return (S1,S2)

def getMatchCost(ScoringMatrix, a,b):
    return ScoringMatrix[a][b]
    
def getOptimalHorizontalGaps( Mat, i, j, gap):
    return max( [ Mat[i][j-k] + gap(k) for k in range(1,j+1) ] )

def getOptimalVerticalGaps( Mat, i, j, gap):
    return max( [ Mat[i-k][j] + gap(k) for k in range(1,i+1) ] )

def computeSimilarityMatrix(S1, S2, Mat, ScoringMatrix, gap = lambda g: g):
    m = len(S1)
    n = len(S2)
        
    for i in range(1, m+1):
        for j in range(1,n+1):
            matchCost = getMatchCost(ScoringMatrix, S1[i-1],S2[j-1])
            Mat[i][j] = max([ 0,
                              Mat[i-1][j-1] + matchCost,
                              Mat[i][j-1] + gap(1),
                              Mat[i-1][j] + gap(1)
                            ])

def getSimilarity(S1, S2, SimilarityMatrix):
    m = len(S1)
    n = len(S2)
    return SimilarityMatrix[m][n]
    
#Linear gap cost by default.
def initializeSimilarityMatrix(S1, S2, gap = lambda g: g ):
    m = len(S1)
    n = len(S2)
    
    SimilarityMatrix = []

    for i in range(m+1):
        SimilarityMatrix.append([])
        for j in range(n+1):
            SimilarityMatrix[i].append(0);
    
    for i in range(1,m+1): #GAP cost for the first column
        SimilarityMatrix[i][0] = gap(i)
    
    for j in range(1,n+1): #GAP cost for the first row
        SimilarityMatrix[0][j] = gap(j)
    
    return SimilarityMatrix

def getAlignment(S1, S2, SimilarityMatrix, scoringMatrix, P, Q, gap):
    Alignment_S1 = ""
    Alignment_S2 = ""
    
    m = len(S1)
    n = len(S2)
    S1 = "-" + S1
    S2 = "-" + S2
    
    while( m != 0 or n != 0 ):
        if( m != 0 and P[m][n] == SimilarityMatrix[m][n] ):
            Alignment_S1 = S1[m] + Alignment_S1
            Alignment_S2 = "-" + Alignment_S2
            while(P[m][n] == (P[m - 1][n] + gap(2) - gap(1)) ):
                m = m - 1
                Alignment_S1 = S1[m] + Alignment_S1
                Alignment_S2 = "-" + Alignment_S2
            m = m - 1
        elif( n != 0 and Q[m][n] == SimilarityMatrix[m][n]):
            Alignment_S1 = "-" + Alignment_S1
            Alignment_S2 = S2[n] + Alignment_S2
            while(Q[m][n] == (Q[m][n - 1] + gap(2) - gap(1))):
                n = n - 1
                Alignment_S1 = "-" + Alignment_S1
                Alignment_S2 = S2[n] + Alignment_S2
            n = n - 1
        elif( m != 0 and n != 0):
            matchCost = scoringMatrix[S1[m]][S2[n]]
            if(SimilarityMatrix[m][n] == SimilarityMatrix[m-1][n-1] + matchCost):
                Alignment_S1 = S1[m] + Alignment_S1
                Alignment_S2 = S2[n] + Alignment_S2
                m = m - 1
                n = n - 1

    return Alignment_S1, Alignment_S2

def getBestAlignment(S1, S2, SimilarityMatrix, scoringMatrix, gap):
    Mat = SimilarityMatrix
    
    for i in range(len(Mat)):
        for j in range(len(Mat[0])):
            Mat[i][j] = (Mat[i][j], i, j)
    
    maxPerRow = [ max(x, key = lambda a: a[0]) for x in Mat]
    maxScore = max(maxPerRow, key = lambda a: a[0]) #Best Alignment score
    
    for i in range(len(Mat)):
        for j in range(len(Mat[0])):
            Mat[i][j] = Mat[i][j][0]
    
    i = maxScore[1]
    j = maxScore[2]
    Alignment_S1 = ""
    Alignment_S2 = ""
    S1 = "-" + S1
    S2 = "-" + S2
    while(Mat[i][j] != 0):
        if(Mat[i][j] == Mat[i-1][j-1] + scoringMatrix[ S1[i] ][ S2[j] ]):
            Alignment_S1 = S1[i] + Alignment_S1
            Alignment_S2 = S2[j] + Alignment_S2
            i = i - 1
            j = j - 1
        elif(Mat[i][j] == Mat[i-1][j] + gap(1)):
            Alignment_S1 = S1[i] + Alignment_S1
            Alignment_S2 = "-" + Alignment_S2
            i = i - 1
        elif(Mat[i][j] == Mat[i][j-1] + gap(1)):
            Alignment_S1 = "-" + Alignment_S1
            Alignment_S2 = S2[j] + Alignment_S2
            j = j - 1

    return maxScore[0], Alignment_S1, Alignment_S2

def gapPenalty(i):
    #Linear gap penalty
    return -5 * i
    
if( __name__ == "__main__"):
    if(len(sys.argv) < 3 ):
        print("usage: python3 LocalAlignment.py <scoring_matrix_file_name> <data_set_file_name> ")
        sys.exit(-1)

    scoringMatrix_filename = sys.argv[1]
    scoringData = ""
    with open(scoringMatrix_filename,"r") as f:
        scoringData = f.read()
        
    dataset_filename = sys.argv[2]
    data = ""
    with open(dataset_filename,"r") as f:
        data = f.read()
    
    scoringMatrix = getScoringMatrix(scoringData)
    S1, S2 = getSequencesFromFastaFormat(data)
    
    #Data Input Completes....
    
    SimilarityMatrix = initializeSimilarityMatrix(S1,S2, gapPenalty)
    computeSimilarityMatrix(S1, S2, SimilarityMatrix, scoringMatrix, gapPenalty)
    alignmentScore, Alignment_S1, Alignment_S2 = getBestAlignment(S1, S2, SimilarityMatrix, scoringMatrix, gapPenalty)
    
    print(alignmentScore)
    print(Alignment_S1.replace("-", ""))
    print(Alignment_S2.replace("-", ""))
