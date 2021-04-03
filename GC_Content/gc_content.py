import sys

def getBestGCContent(dna_sequence):
    sequence_list = dna_sequence.split('>')[1:] 
    gc_list = []
    for seq in sequence_list:
        dna_id = seq.split('\n')[0]
        dna = "".join(seq.split('\n')[1:])
        gc = 100*(dna.count('G') + dna.count('C')) / len(dna)
        gc_list.append((round(gc,6), dna_id))
    gc_list.sort(reverse = True)
    return gc_list[0]
    # return (60.919540,"Rosalind_0808")

if( __name__ == "__main__"):
    if(len(sys.argv) < 2 ):
        print("usage: python3 transcription.py <data_set_file_name>")
        sys.exit(-1)

    dataset_filename = sys.argv[1]
    dna_sequence = ""

    with open(dataset_filename,"r") as f:
        dna_sequence = f.read()

    gc, dna_id = getBestGCContent(dna_sequence)
    print(dna_id)
    print(gc)
