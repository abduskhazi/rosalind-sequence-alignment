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

def getSequenceFromFastaFormat(data):
    _, E1 = data.split(">")
    S1 = "".join(E1.strip().split("\n")[1:])
    return S1

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

def getProteinSeq(dna):
    proteinSeqs = ["Stop"]
    
    for i in range(3):
        for j in range(i, len(dna), 3):
            if( j + 2 < len(dna) ):
                codon = dna[j:j+3]
                if(proteinSeqs[-1][-4:] == "Stop" and codonTable[codon] == "M"):
                    proteinSeqs.append("")
                if(proteinSeqs[-1][-4:] != "Stop"):
                    proteinSeqs[-1] += codonTable[codon]
        if(proteinSeqs[-1][-4:] != "Stop"):
            proteinSeqs[-1] = "Stop"
    
    proteinSeqs = [x[:-4] for x in proteinSeqs if x != "Stop"]
    
    extraSeq = []
    for seq in proteinSeqs:
        indices = [i for i in range(1,len(seq)) if seq[i] == 'M']
        for i in indices:
            extraSeq.append(seq[i:])

    return proteinSeqs + extraSeq
    
def translate(dna, codonTable):
    revDNA = ""
    for i in range(len(dna)):
        if(dna[i] == "A"):
           revDNA = "T" + revDNA
        elif(dna[i] == "T"):
           revDNA = "A" + revDNA
        elif(dna[i] == "G"):
           revDNA = "C" + revDNA
        elif(dna[i] == "C"):
           revDNA = "G" + revDNA
    
    return set(getProteinSeq(dna) + getProteinSeq(revDNA))
    
if( __name__ == "__main__"):
    if(len(sys.argv) < 3 ):
        print("usage: python3 DNA_to_protein.py <DNA_CodonTable_filename> <DNA_file_name>")
        sys.exit(-1)

    DNA_CodonTable_filename = sys.argv[1]
    codonTable = getmRNACodonTable( DNA_CodonTable_filename )

    dna_filename = sys.argv[2]
    dna_data_fasta = ""
    with open(dna_filename,"r") as f:
        dna_data_fasta = f.read()

    dna = getSequenceFromFastaFormat(dna_data_fasta)
        
    proteinSeqs = translate(dna, codonTable)

    for prt in proteinSeqs:
        print(prt)

