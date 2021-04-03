# Multiple Sequence Alignment using Needleman-Wunsh
# algorithm with linear gap penalty.
# Considering all possible gaps?

import sys

def getSequencesFromFastaFormat(data):
    E = data.split(">")
    sequencelist = []
    for i in range(1,len(E)):
        S = "".join(E[i].strip().split("\n")[1:])
        sequencelist.append(S)
    return sequencelist

def getNDimentionalMatrix(sequencelist, Mat = []):
    if(len(sequencelist) == 0):
        return 0
    n = len(sequencelist[0])
    for i in range(n+1):
        Mat.append([])
        Mat[i] = getNDimentionalMatrix(sequencelist[1:], Mat[i])
        
    return Mat
            

def initializeSimilarityMatrix(sequencelist, gap = lambda g: -g ):

    SimilarityMatrix = getNDimentionalMatrix(sequencelist) #Initializing n-Dimentional Matrix

    return SimilarityMatrix
    

def score(S1, S2, S3, S4):
    scorelist = [S1, S2, S3, S4]
    score = 0
    for i in range(len(scorelist)):
        for j in range(i+1,len(scorelist)):
            if( scorelist[i] != scorelist[j] ):
                score += -1
    return(score)

def computeSimilarityMatrix(sequencelist, Mat, gap = lambda g: -g):
    S1, S2, S3, S4 = sequencelist
    
    a = len(S1)
    b = len(S2)
    c = len(S3)
    d = len(S4)

    S1 = "-" + S1
    S2 = "-" + S2
    S3 = "-" + S3
    S4 = "-" + S4

    for i in range(0, a+1):
        for j in range(0,b+1):
            for k in range(0,c+1):
                for l in range(0,d+1):
                    scorelist = []
                    if(i!=0 and j!=0 and k!=0 and l!= 0):
                        scorelist.append (Mat[i-1][j-1][k-1][l-1] + score( S1[i], S2[j], S3[k], S4[l] ) )
                    if(i!=0 and j!=0 and k!=0 and True ):
                        scorelist.append( Mat[i-1][j-1][k-1][l  ] + score( S1[i], S2[j], S3[k],   "-" ) )
                    if(i!=0 and j!=0 and True and l != 0):
                        scorelist.append( Mat[i-1][j-1][k  ][l-1] + score( S1[i], S2[j],   "-", S4[l] ) )
                    if(i!=0 and j!=0 and True and True ):
                        scorelist.append( Mat[i-1][j-1][k  ][l  ] + score( S1[i], S2[j],   "-",   "-" ) )
                    if(i!=0 and True and k!=0 and l!=0 ):
                        scorelist.append( Mat[i-1][j  ][k-1][l-1] + score( S1[i],   "-", S3[k], S4[l] ) )
                    if(i!=0 and True and k!=0 and True ):
                        scorelist.append( Mat[i-1][j  ][k-1][l  ] + score( S1[i],   "-", S3[k],   "-" ) )
                    if(i!=0 and True and True and l!=0 ):
                        scorelist.append( Mat[i-1][j  ][k  ][l-1] + score( S1[i],   "-",   "-", S4[l] ) )
                    if(i!=0 and True and True and True ):
                        scorelist.append( Mat[i-1][j  ][k  ][l  ] + score( S1[i],   "-",   "-",   "-" ) )
                    if(True and j!=0 and k!=0 and l!=0 ):
                        scorelist.append( Mat[i  ][j-1][k-1][l-1] + score(   "-", S2[j], S3[k], S4[l] ) )
                    if(True and j!=0 and k!=0 and True ):
                        scorelist.append( Mat[i  ][j-1][k-1][l  ] + score(   "-", S2[j], S3[k],   "-" ) )
                    if(True and j!=0 and True and l!=0 ):
                        scorelist.append( Mat[i  ][j-1][k  ][l-1] + score(   "-", S2[j],   "-", S4[l] ) )
                    if(True and j!=0 and True and True):
                        scorelist.append( Mat[i  ][j-1][k  ][l  ] + score(   "-", S2[j],   "-",   "-" ) )
                    if(True and True and k!=0 and l!=0 ):
                        scorelist.append( Mat[i  ][j  ][k-1][l-1] + score(   "-",   "-", S3[k], S4[l] ) )
                    if(True and True and k!=0 and True ):
                        scorelist.append( Mat[i  ][j  ][k-1][l  ] + score(   "-",   "-", S3[k],   "-" ) )
                    if(True and True and True and l!=0 ):
                        scorelist.append( Mat[i  ][j  ][k  ][l-1] + score(   "-",   "-",   "-", S4[l] ) )
                    
                    if(i!=0 or j!=0 or k!=0 or l!= 0):
                        Mat[i][j][k][l] = max(scorelist)

def getSimilarity(sequencelist, Mat):
    S1, S2, S3, S4 = sequencelist
    a = len(S1)
    b = len(S2)
    c = len(S3)
    d = len(S4)
    return Mat[a][b][c][d]

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

if( __name__ == "__main__"):
    if(len(sys.argv) < 2 ):
        print("usage: python3 multipleAlignment.py <fasta_file_name>")
        sys.exit(-1)

    dataset_filename = sys.argv[1]
    data = ""

    with open(dataset_filename,"r") as f:
        data = f.read()

    sequencelist = getSequencesFromFastaFormat(data)
    print(sequencelist)
    
    SimilarityMatrix = initializeSimilarityMatrix(sequencelist) # This is a n-D Cuboid.
    computeSimilarityMatrix(sequencelist, SimilarityMatrix)
    
    Similarity = getSimilarity(sequencelist, SimilarityMatrix)
    # Align_S1, Align_S2 = getOptimalAlignment(S1, S2, EditMatrix)
    
    print(Similarity)
    """
    print()
    print(Align_S1)
    print()
    print(Align_S2)
    """
