import sys

def complement(ch):
    if(ch == 'A'):
        return "T"
    if(ch == 'T'):
        return "A"
    if(ch == 'C'):
        return "G"
    if(ch == 'G'):
        return "C"
    return ""

if( __name__ == "__main__"):
    if(len(sys.argv) < 2 ):
        print("usage: python3 transcription.py <data_set_file_name>")
        sys.exit(-1)

    dataset_filename = sys.argv[1]
    dna_sequence = ""

    with open(dataset_filename,"r") as f:
        dna_sequence = f.read()

    complement_sequence = ""
    for ch in reversed(dna_sequence):
        complement_sequence += complement(ch)

    print(complement_sequence)
