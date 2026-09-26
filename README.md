# Rosalind Sequence Alignment

**Bioinformatics Algorithms — Implemented From Scratch**  
**Result**: 9 DP alignment algorithms + 8 sequence utilities, all verified against Rosalind hidden test cases — zero library calls.

---

## What's Here
| Alignment (9) | Utilities (8) |
|---------------|---------------|
| Global / Local / Semi-Global | GC Content |
| Edit Distance + Reconstruction | Point Mutations (Hamming) |
| Constant Gap Penalty | DNA → RNA (transcription) |
| Affine Gap / Affine Local | RNA → Protein (translation) |
| Multiple Alignment (progressive) | DNA Rev Complement |
| Scoring Matrices (BLOSUM/PAM) | DNA → Protein (6-frame) |
| | Protein → mRNA (back-translation) |
| | Nucleotide helpers |

---

## What Implementing Teaches You
- **Affine gaps = 3 DP matrices** (M, Ix, Iy). Gap open vs extend changes everything.
- **Local alignment = reset to 0**. Finding subsequences, not full alignment.
- **Space optimization** (Hirschberg): 2 rows for score, full matrix for traceback.
- **Multiple alignment is NP-hard**. Progressive heuristic (guide tree + pairwise) is the only practical way.

---

## Verification
Every algorithm passes Rosalind's hidden test cases. Sample inputs in `rosalind_*.txt`.

---

## Layout
```
├── *Alignment/           # 9 algos
├── GC_Content/  PointMutations/
├── DNA_to_RNA/  RNA_to_protein/
├── DNA_strand_Comp/  DNA_to_proteins/
├── Proteins_to_mRNA/  Nucleotids.py
└── rosalind_*.txt
```

---

## Run Any Algorithm
```bash
cd GlobalAlignment && python main.py < ../rosalind_*.txt
```

---

**Author**: Abdus Salam Khazi  
**Contact**: abduskhazi@gmail.com