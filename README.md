# 🧬 CS 482 – Computational Biology Projects

This repository contains four major assignments completed for **CS 482: Computational Biology** at the University of Waterloo. The projects explore core methods in bioinformatics, including **sequence alignment**, **peptide classification**, **variant analysis**, and **mass spectrometry–based peptide identification**, using Python and widely used computational biology tools.

---

## 📁 Assignment 1: Pairwise Sequence Alignment

**Goal:**  
Align a query DNA sequence against a set of target sequences using a custom dynamic programming algorithm.

### Features:
- Global alignment with customizable scoring:
  - Match = +1  
  - Mismatch = -1  
  - Gap = -1
- Reconstructs aligned sequences via traceback
- Parses input FASTA files using `Biopython`
- Outputs best match and alignment score to a text file

## 📁 Assignment 2: Peptide Naturalness Classification

**Goal:**  
Classify peptides as biologically "natural" or randomly shuffled using statistical scoring and optional machine learning models.

### 🔧 Features:
- **Score1:** Computes a log-likelihood score based on **k-mer frequency** (typically k=3) comparing natural vs. random peptide datasets.
- **Score2:** Implements a user-defined method, optionally using **machine learning** (e.g. logistic regression with `scikit-learn`) to improve classification.
- Uses **pseudocount smoothing** to avoid divide-by-zero and ensure robustness when k-mers are rare or unseen.
## 📁 Assignment 3: SARS-CoV-2 Variant RBD Extraction & Alignment

**Goal:**  
Extract the **receptor-binding domain (RBD)** of the spike protein from various **SARS-CoV-2 variants** and perform a **multiple sequence alignment (MSA)** to identify and visualize mutations relative to the original Wuhan-Hu-1 strain.

---

### 🔬 Background

The spike (S) protein of SARS-CoV-2 contains a receptor-binding domain (RBD) responsible for attaching to human ACE2 receptors. Monitoring mutations in the RBD is critical for understanding variant behavior and vaccine resistance.

---

### 🔧 Features

- Retrieves complete spike protein sequences for the following variants:
  - **Alpha (B.1.1.7)**
  - **Beta (B.1.351)**
  - **Gamma (P.1)**
  - **Delta (B.1.617.2)**
  - **Omicron (B.1.1.529)**
- Extracts the **RBD subsequence** from each variant using annotated residue ranges (e.g., Wuhan-Hu-1: positions 330–583)
- Aligns sequences using **Clustal Omega** to compare mutations
- Highlights differences from the reference strain (Wuhan-Hu-1)

---

### 📄 Input/Output

#### ✅ Input:
- Variant names (used to query NCBI Protein database)
- Full spike protein sequences in FASTA format (downloaded manually or programmatically)

#### 📤 Output:
- RBD sequences in FASTA format
- Multiple sequence alignment (MSA) highlighting mutated residues

---
## 📁 Assignment 4: Peptide-Spectrum Matching with Cosine Similarity

**Goal:**  
Match experimental **MS/MS spectra** to candidate peptides (target and decoy) by comparing theoretical and observed spectra using **cosine similarity**.

---

### 🔬 Background

Mass spectrometry-based proteomics relies on matching observed spectra to theoretical fragment spectra generated from candidate peptides. This assignment implements a scoring-based pipeline for identifying the best-matching peptide for each input spectrum.

---

### 🔧 Features

- **Reads Input Data:**
  - MS/MS spectra from `.mgf` files
  - Candidate peptide lists (target and decoy)
- **Generates Theoretical Spectra:**
  - Calculates **b- and y-ion fragment masses** for each peptide
  - Assumes charge +1 for simplification
- **Discretizes Spectra:**
  - Both experimental and theoretical spectra are binned using a fixed **bin size = 0.5**
  - Intensity-normalized vectors are created for cosine comparison
- **Precursor Mass Filtering:**
  - Skips comparisons if the **precursor mass difference > 0.1**, improving efficiency
- **Scoring:**
  - Uses **cosine similarity** between binned spectra vectors to determine best match
- **Postprocessing:**
  - Reports number of identifications at **1% false discovery rate (FDR)**

---
