import sys

if(len(sys.argv) < 2 ):
    print("usage: python3 transcription.py <data_set_file_name>")
    sys.exit(-1)

dataset_filename = sys.argv[1]
dna_sequence = ""

with open(dataset_filename,"r") as f:
    dna_sequence = f.read()

print(dna_sequence.replace("T","U"))
