# Since edit distance is the same as alignment distance
# it just makes sense to compute the global alignment
# distance

import sys
import math

"""
BLOSUM62 =
  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
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
    
    P = [] #Storing the gap history of S2
    Q = [] #Storing the gap history of S1
    P.append( [-math.inf] * (n+1) )
    for i in range(0,m+1):
        if(i != 0):
            P.append( [0] * (n+1) )
        Q.append( [0] * (n+1) )
        Q[-1][0] = -math.inf

    maxScore = (Mat[0][0], 0, 0)
    
    for i in range(1, m+1):
        for j in range(1,n+1):
            P[i][j] = max([ Mat[i-1][j] + gap(1),
                            P[i-1][j] + gap(2) - gap(1)
                          ])
            Q[i][j] = max([ Mat[i][j-1] + gap(1),
                            Q[i][j-1] + gap(2) - gap(1)
                          ])
            matchCost = getMatchCost(ScoringMatrix, S1[i-1],S2[j-1])
            Mat[i][j] = max([ 0,
                              Mat[i-1][j-1] + matchCost,
                              P[i][j], # getOptimalHorizontalGaps(Mat, i, j, gap),
                              Q[i][j] # getOptimalVerticalGaps(Mat, i, j, gap),
                            ])
                            
            if(maxScore[0] < Mat[i][j]):
                maxScore = (Mat[i][j], i, j)
                
    return P, Q, maxScore

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
        SimilarityMatrix[i][0] = 0 # gap(i) now local alignment
    
    for j in range(1,n+1): #GAP cost for the first row
        SimilarityMatrix[0][j] = 0 # gap(i) now local alignment
    
    return SimilarityMatrix

def getAlignment(S1, S2, m, n, SimilarityMatrix, ScoringMatrix, P, Q, gap):
    Alignment_S1 = ""
    Alignment_S2 = ""

    S1 = "-" + S1
    S2 = "-" + S2
    
    while( SimilarityMatrix[m][n] != 0 ):
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
            matchCost = getMatchCost(ScoringMatrix, S1[m],S2[n])
            if(SimilarityMatrix[m][n] == SimilarityMatrix[m-1][n-1] + matchCost):
                Alignment_S1 = S1[m] + Alignment_S1
                Alignment_S2 = S2[n] + Alignment_S2
                m = m - 1
                n = n - 1

    return Alignment_S1, Alignment_S2

def getBestAlignmentData(SimilarityMatrix):
    Mat = SimilarityMatrix
    
    maxScore = (Mat[0][0], 0, 0)
    for i in range(len(Mat)):
        for j in range(len(Mat[0])):
            if(maxScore[0] < Mat[i][j]):
                maxScore = (Mat[i][j], i, j)
    
    maxPerRow = [ max(x, key = lambda a: a[0]) for x in Mat]
    maxScore = max(maxPerRow, key = lambda a: a[0]) #Best Alignment score
    
    for i in range(len(Mat)):
        for j in range(len(Mat[0])):
            Mat[i][j] = Mat[i][j][0]
    
    return maxScore

def gapPenalty(i):
    if(i == 1):
        return -11
    return -11 - (i-1)*1
    
if( __name__ == "__main__"):
    if(len(sys.argv) < 3 ):
        print("usage: python3 AffineLocalAlignment.py <scoring_matrix_file_name> <data_set_file_name> ")
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
    P, Q, maxScore = computeSimilarityMatrix(S1, S2, SimilarityMatrix, scoringMatrix, gapPenalty)
    alignmentScore, m, n = maxScore
    Alignment_S1, Alignment_S2 = getAlignment(S1, S2, m, n, SimilarityMatrix, scoringMatrix, P, Q, gapPenalty)
    
    print(alignmentScore)
    print(Alignment_S1.replace("-", ""))
    print(Alignment_S2.replace("-", ""))
