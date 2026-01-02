# WPGMA Phylogenetic Tree (From Scratch, R)

This repository contains an implementation of the
**WPGMA (Weighted Pair Group Method with Arithmetic Mean)**
algorithm using **Jukes–Cantor distances** to construct
a phylogenetic tree.

---

## 🧬 Features
- Pairwise Jukes–Cantor distance calculation
- WPGMA hierarchical clustering
- Manual distance matrix updates
- Tree construction and visualization
- Designed for instructional use

---

## 🚀 Usage

```r
source("wpgma_phylogeny_from_scratch.R")

## ⏱️ Runtime & Space Complexity

Let **n** be the number of sequences and **L** be the sequence length.

### Distance Calculation
- Computing Jukes–Cantor distance between two sequences takes **O(L)** time.
- Computing all pairwise distances takes **O(n² × L)** time.

### WPGMA Clustering
- Each iteration finds the minimum distance in the matrix: **O(n²)**.
- The distance matrix is updated **n − 1** times.
- Overall clustering time complexity: **O(n³)**.

### Space Complexity
- Distance matrix storage requires **O(n²)** space.
- Tree and bookkeeping structures require **O(n)** space.
- Total space complexity: **O(n²)**.

### Practical Notes
- WPGMA assigns equal weight to merged clusters regardless of cluster size.
- This makes WPGMA sensitive to unequal sampling and less biologically realistic
  than UPGMA in many scenarios.
- Like UPGMA, WPGMA assumes a molecular clock (ultrametric distances).
