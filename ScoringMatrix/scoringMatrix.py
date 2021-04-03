# Since edit distance is the same as alignment distance
# it just makes sense to compute the global alignment
# distance

import sys

"""
S =
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
    _, E1,E2 = data.split(">")
    S1 = "".join(E1.strip().split("\n")[1:])
    S2 = "".join(E2.strip().split("\n")[1:])
    return (S1,S2)

def getMatchCost(ScoringMatrix, a,b):
    return ScoringMatrix[a][b]
    
def computeSimilarityMatrix(S1, S2, Mat, ScoringMatrix, gap = lambda g: g):
    m = len(S1)
    n = len(S2)
    for i in range(1, m+1):
        for j in range(1,n+1):
            matchCost = getMatchCost(ScoringMatrix, S1[i-1],S2[j-1])
            Mat[i][j] = max([ Mat[i][j-1] + gap(1) ,
                              Mat[i-1][j] + gap(1) ,
                              Mat[i-1][j-1] + matchCost]
                            )

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
        
if( __name__ == "__main__"):
    if(len(sys.argv) < 3 ):
        print("usage: python3 scoringMatrix.py <data_set_file_name> <scoring_matrix_file_name>")
        sys.exit(-1)

    dataset_filename = sys.argv[1]
    data = ""
    with open(dataset_filename,"r") as f:
        data = f.read()

    scoringMatrix_filename = sys.argv[2]
    scoringData = ""
    with open(scoringMatrix_filename,"r") as f:
        scoringData = f.read()
    
    scoringMatrix = getScoringMatrix(scoringData)
    S1, S2 = getSequencesFromFastaFormat(data)
    
    SimilarityMatrix = initializeSimilarityMatrix(S1,S2, gap = lambda g : -5 * g)
    computeSimilarityMatrix(S1, S2, SimilarityMatrix, scoringMatrix,
                            gap = lambda g : -5 * g)
    
    print(getSimilarity(S1, S2, SimilarityMatrix))
