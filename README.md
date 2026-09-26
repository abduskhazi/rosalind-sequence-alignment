# Rosalind Sequence Alignment

> Implementations of **sequence-alignment algorithms from first principles**, developed through the Rosalind bioinformatics problem set.

This repository explores how dynamic programming and scoring models can be used to compare biological sequences and find optimal alignments.

```text
Sequences
    ↓
Scoring Model
    ↓
Dynamic Programming
    ↓
Optimal Alignment
```

## Implemented

The repository contains implementations covering:

* **Global, local and semi-global alignment**
* **Edit distance and alignment reconstruction**
* **Constant and affine gap penalties**
* **Multiple sequence alignment**
* **Scoring matrices**
* Related DNA/RNA and protein sequence transformations

The more advanced alignment variants explore how changing the **scoring model** changes the underlying dynamic-programming formulation.

## Key Idea

A sequence-alignment problem can be viewed as a path through an alignment matrix:

```text
        Sequence B
      ──────────────→
     ┌───────────────┐
Seq A│   DP Matrix   │
     │               │
     │       ↘       │
     │         ↘     │
     └───────────────┘
```

Matches, mismatches and gaps correspond to different transitions, while the dynamic-programming recurrence determines the best-scoring path.

The repository was built to understand these algorithms by implementing them directly rather than treating alignment as a black-box library operation.

## Repository Structure

Each major alignment technique is kept in its own directory, including:

```text
AffineGapPenalty/
AffineLocalAlignment/
ConstantGapPenalty/
EditDistance/
EditDistanceAlignment/
LocalAlignment/
MultipleAlignment/
SemiGlobalAlignment/
ScoringMatrix/
```

The project was developed as a collection of solutions to problems from **Rosalind**.

**Abdus Salam Khazi**
