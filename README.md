# 🧬 CS 482 – Computational Biology Projects

This repository contains four major assignments completed for **CS 482: Computational Biology** at the University of Waterloo. The projects explore core topics in bioinformatics, including **sequence alignment**, **peptide classification**, **variant analysis**, and **mass spectrometry–based peptide identification**, using Python and common computational biology tools.

---

## 📁 Project 1: Pairwise Sequence Alignment

**Goal:**  
Align a query DNA sequence against a set of target sequences using a custom **dynamic programming** algorithm.

### 🔧 Features
- Global alignment with customizable scoring:
  - Match = +1  
  - Mismatch = -1  
  - Gap = -1
- Reconstructs aligned sequences via traceback
- Parses FASTA input using `Biopython`
- Outputs best-scoring aligned sequence from a sequence database

---

## 📁 Project 2: Peptide Naturalness Classification

**Goal:**  
Classify peptides as biologically “natural” or randomly shuffled using **k-mer frequency statistics** and optional **machine learning** techniques.

### 🔧 Features
- **Score1:** Computes log-likelihood based on **k-mer frequency** comparisons between natural and random datasets
- **Score2:** Supports a custom classifier (e.g. logistic regression via `scikit-learn`)
- Applies **pseudocount smoothing** to avoid divide-by-zero issues in rare k-mers
- Provides interpretability through scoring breakdown

---

## 📁 Project 3: SARS-CoV-2 Variant RBD Extraction & Alignment

**Goal:**  
Extract the **receptor-binding domain (RBD)** from the spike protein of various **SARS-CoV-2 variants** and perform a **multiple sequence alignment (MSA)** to highlight mutations relative to the original Wuhan-Hu-1 strain.

### 🔧 Features
- Retrieves full spike sequences for:
  - **Alpha**, **Beta**, **Gamma**, **Delta**, and **Omicron**
- Extracts RBD regions using annotated residue positions (e.g. Wuhan-Hu-1: 330–583)
- Performs **MSA** using Clustal Omega
- Highlights amino acid substitutions across variants

### 📄 Input / Output
- **Input:** Variant protein sequences from the NCBI Protein database (in FASTA format)
- **Output:** Aligned RBD sequences showing mutations from the reference

---

## 📁 Project 4: Peptide-Spectrum Matching with Cosine Similarity

**Goal:**  
Match **MS/MS spectra** to candidate peptides using discretized theoretical spectra and **cosine similarity** scoring.

### 🔧 Features
- Reads spectra from `.mgf` files and peptides from target/decoy lists
- Generates **theoretical b- and y-ion spectra**
- Discretizes spectra with fixed **bin size = 0.5**
- Computes **cosine similarity** between observed and theoretical spectra
- Applies **precursor mass filtering** (tolerance ≤ 0.1) for speedup
- Reports number of identifications at **1% FDR**

---

## 🛠 Installation

Install all required packages with:
```bash
pip install numpy scipy biopython scikit-learn
