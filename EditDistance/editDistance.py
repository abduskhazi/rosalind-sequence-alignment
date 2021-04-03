# Since edit distance is the same as alignment distance
# it just makes sense to compute the global alignment
# distance

import sys

def getSequencesFromFastaFormat(data):
    _, E1,E2 = data.split(">")
    S1 = "".join(E1.strip().split("\n")[1:])
    S2 = "".join(E2.strip().split("\n")[1:])
    return (S1,S2)

def computeEditMatrix(S1, S2, Mat, m, n):
    for i in range(1, m+1):
        for j in range(1,n+1):
            matchCost = 1
            if(S1[i-1] == S2[j-1]):
                matchCost = 0
            Mat[i][j] = min([ Mat[i][j-1] + 1 , Mat[i-1][j] + 1, Mat[i-1][j-1] + matchCost])

def getEditDistance(S1, S2):
    m = len(S1)
    n = len(S2)
    
    EditMatrix = []

    for i in range(m+1):
        EditMatrix.append([])
        for j in range(n+1):
            EditMatrix[i].append(0);
    
    for i in range(1,m+1):
        EditMatrix[i][0] = i
    
    for j in range(1,n+1):
        EditMatrix[0][j] = j

    computeEditMatrix(S1, S2, EditMatrix, m, n) 
    return EditMatrix[m][n]

if( __name__ == "__main__"):
    if(len(sys.argv) < 2 ):
        print("usage: python3 transcription.py <data_set_file_name>")
        sys.exit(-1)

    dataset_filename = sys.argv[1]
    data = ""

    with open(dataset_filename,"r") as f:
        data = f.read()

    S1, S2 = getSequencesFromFastaFormat(data)
    print(getEditDistance(S1, S2))
