"""
RNA Codon Table:
UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G
"""

import sys

def getmRNACodonTable( codontable_filename ):
    CodonTable = {}
    codonData = ""
    
    with open(codontable_filename,"r") as f:
        codonData = f.read()

    list = codonData.strip().replace("\n"," ").split(" ")
    list = [x for x in list if x != "" ]
    for i in range(0, len(list), 2):
        CodonTable[list[i]] = list[i+1]

    return CodonTable

def getNum_mRNA(proteinSequence, codonTable):
    numMap = {}
    for k in codonTable:
        prtResidue = codonTable[k]
        if(prtResidue not in numMap):
            numMap[prtResidue] = 0
        numMap[prtResidue] += 1
    numMap["M"] = numMap["Stop"]
    #print(numMap)
    res = 1
    for p in proteinSequence:
        res = (res * numMap[p]) % 1000000
    
    return res

if( __name__ == "__main__"):
    if(len(sys.argv) < 3 ):
        print("usage: python3 protein_to_mRNA.py <mRNA_CodonTable_filename> <protein_file_name>")
        sys.exit(-1)

    mRNA_CodonTable_filename = sys.argv[1]
    codonTable = getmRNACodonTable( mRNA_CodonTable_filename )

    protein_filename = sys.argv[2]
    proteinSequence = ""

    with open(protein_filename,"r") as f:
        proteinSequence = f.read().strip()
        
    num = getNum_mRNA(proteinSequence, codonTable)

    print(num)

