import sys

dataset_filename = sys.argv[1]
dna_sequence = ""

with open(dataset_filename,"r") as f:
    dna_sequence = f.read()

print("%d %d %d %d" % (dna_sequence.count('A'), \
                       dna_sequence.count('C'), \
                       dna_sequence.count('G'), \
                       dna_sequence.count('T')) \
      )
