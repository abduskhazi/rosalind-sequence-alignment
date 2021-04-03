#!/usr/bin/env pypy3

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

def ribosome(mRNASequence, codonTable):
    proteinSequence = ""
    
    for i in range(0, len(mRNASequence), 3):
        codon = mRNASequence[i:i+3]
        if(codonTable[codon] != "Stop"):
            proteinSequence = proteinSequence + codonTable[codon]
    
    return proteinSequence
    
if( __name__ == "__main__"):
    if(len(sys.argv) < 2 ):
        print("usage: python3 rna_to_protein.py <mRNA_CodonTable_filename> <mRNA_file_name>")
        sys.exit(-1)

    mRNA_CodonTable_filename = sys.argv[1]
    codonTable = getmRNACodonTable( mRNA_CodonTable_filename )

    mRNA_filename = sys.argv[2]
    mRNASequence = ""

    with open(mRNA_filename,"r") as f:
        mRNASequence = f.read().strip()
        
    proteinSequence = ribosome(mRNASequence, codonTable)

    print(proteinSequence)
