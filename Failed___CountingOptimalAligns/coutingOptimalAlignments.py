# Since edit distance is the same as alignment distance
# it just makes sense to compute the global alignment
# distance

import sys

def getSequencesFromFastaFormat(data):
    _, E1,E2 = data.split(">")
    S1 = "".join(E1.strip().split("\n")[1:])
    S2 = "".join(E2.strip().split("\n")[1:])
    return (S1,S2)

def computeEditDistanceMatrix(S1, S2, Mat):
    m = len(S1)
    n = len(S2)
    for i in range(1, m+1):
        for j in range(1,n+1):
            matchCost = 1
            if(S1[i-1] == S2[j-1]):
                matchCost = 0
            Mat[i][j] = min([ Mat[i][j-1] + 1 , Mat[i-1][j] + 1, Mat[i-1][j-1] + matchCost])

def getEditDistance(S1, S2, EditMatrix):
    m = len(S1)
    n = len(S2)
    return EditMatrix[m][n]
    
#Linear gap cost by default.
def initializeEditMatrix(S1, S2, gap = lambda g: g ):
    m = len(S1)
    n = len(S2)
    
    EditMatrix = []

    for i in range(m+1):
        EditMatrix.append([])
        for j in range(n+1):
            EditMatrix[i].append(0);
    
    for i in range(1,m+1): #GAP cost for the first column
        EditMatrix[i][0] = gap(i)
    
    for j in range(1,n+1): #GAP cost for the first row
        EditMatrix[0][j] = gap(j)
    
    return EditMatrix

def getOptimalAlignment(S1_, S2_, EditMatrix):
    Align_S1 = ""
    Align_S2 = ""
    i = len(S1_)
    j = len(S2_)
    S1 = "-" + S1_
    S2 = "-" + S2_
    while( i != 0 or j != 0 ):
        Elem = EditMatrix[i][j]
        if  ( i != 0 and Elem == (EditMatrix[i-1][j] + 1) ):
            Align_S1 = S1[i] + Align_S1
            Align_S2 = "-" + Align_S2
            i = i - 1
        elif( j != 0 and Elem == (EditMatrix[i][j-1] + 1) ):
            Align_S1 = "-" + Align_S1
            Align_S2 = S2[j] + Align_S2
            j = j - 1
        else:
            matchCost = 1
            if(S1[i] == S2[j]):
                matchCost = 0
            if(Elem == (EditMatrix[i-1][j-1] + matchCost) ):
                Align_S1 = S1[i] + Align_S1
                Align_S2 = S2[j] + Align_S2
                i = i - 1
                j = j - 1
            
    return Align_S1, Align_S2

def countOptimalAlignments(S1, S2, m, n, EditMatrix, cacheList ):
    if( m == 0 or n == 0):
        return 1

    count = 0
    Elem = EditMatrix[m][n]
    if( Elem == (EditMatrix[m-1][n] + 1) ):
        if(cacheList[m-1][n] == -1):
            cacheList[m-1][n] = countOptimalAlignments(S1, S2, m-1, n, EditMatrix, cacheList)
        count += cacheList[m-1][n]
        count = count % 134217727
        
    if( Elem == (EditMatrix[m][n-1] + 1) ):
        if( cacheList[m][n-1] == -1):
            cacheList[m][n-1] = countOptimalAlignments(S1, S2, m, n-1, EditMatrix, cacheList)
        count += cacheList[m][n-1]
        count = count % 134217727
    
    matchCost = 1
    if(S1[m-1] == S2[n-1]):
        matchCost = 0
    if( Elem == (EditMatrix[m-1][n-1] + matchCost) ):
        if( cacheList[m-1][n-1] == -1):
            cacheList[m-1][n-1] = countOptimalAlignments(S1, S2, m-1, n-1, EditMatrix, cacheList)
        count += cacheList[m-1][n-1]
        count = count % 134217727
        
    return count
        
if( __name__ == "__main__"):
    if(len(sys.argv) < 2 ):
        print("usage: python3 coutingOptimalAlignments.py <data_set_file_name>")
        sys.exit(-1)

    dataset_filename = sys.argv[1]
    data = ""

    with open(dataset_filename,"r") as f:
        data = f.read()

    S1, S2 = getSequencesFromFastaFormat(data)
    
    EditMatrix = initializeEditMatrix(S1,S2)
    computeEditDistanceMatrix(S1, S2, EditMatrix)
    
    # print(len(S1), len(S2))
    m = len(S1)
    n = len(S2)
    count = countOptimalAlignments(S1, S2, len(S1), len(S2), EditMatrix,  [[-1] * n] * m )
    print(count)
