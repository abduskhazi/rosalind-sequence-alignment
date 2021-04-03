import sys

def hamming_distance(dna1,dna2):
    distance = 0
    for i in range(len(dna1)):
        if(dna1[i] != dna2[i]):
            distance += 1
    return distance

if( __name__ == "__main__"):
    if(len(sys.argv) < 2 ):
        print("usage: python3 transcription.py <data_set_file_name>")
        sys.exit(-1)

    dataset_filename = sys.argv[1]
    dna_sequence = ""

    with open(dataset_filename,"r") as f:
        sequence = f.read()

    dna1, dna2 = sequence.strip().split('\n')
    print(hamming_distance(dna1,dna2))
