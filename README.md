
# Project Overview

This project evaluates structural similarity-based search methods for identifying Diels-Alder (DA) enzyme homologs, with a particular focus on spirotetronate cyclases. The study contrasts traditional sequence similarity searches with advanced structural methods such as Foldseek and AlphaFold clusters to assess their effectiveness in retrieving relevant enzyme homologs.

# Key Objectives

1. **Evaluate Search Methods**: Compare traditional sequence similarity searches against structural similarity methods.
2. **Performance Analysis**:
   - **Method A**: Applies Foldseek in conjunction with TM-align and 3Di alignment methods.
   - **Method B**: Utilizes AlphaFold clusters for a structured clustering approach.
3. **Identify Relevant Homologs**: Assess the ability of each method to successfully identify relevant DA homologs.
4. **Limitations Analysis**: Explore the limitations of relying solely on structural similarity and highlight the benefits of incorporating sequence similarity.
5. **Parameter Testing**: Investigate challenges in predicting closed conformations for DA enzymes using AlphaFold.

# Data Sources

- **Foldseek Results**: m8 format files generated from Foldseek searches against the AlphaFold proteome database (AFDB50) and SwissProt databases.
- **UniProt Database**: Used for retrieving protein sequences based on UniProt IDs.

# Methodology

The project employs a series of Python scripts to process and analyze data obtained from structural searches. These scripts include:

- Extracting UniProt IDs from Foldseek results.
- Finding unique IDs from the results.
- Retrieving protein sequences in FASTA format from UniProt IDs.

## Visualization and Analysis

Visualizations of the structural data and clustering results were generated using PyMOL and Cytoscape, allowing for effective representation of the enzyme homolog structures.

# Results

The study concluded that combining structural and sequence-based methods enhances the discovery of DA enzyme homologs. Method B (AlphaFold clusters) demonstrated higher reliability in retrieving accurate hits compared to Method A (Foldseek).

**Key Findings**:
- Foldseek's performance is influenced by the choice of alignment method.
- Incorporating sequence similarity improves search comprehensiveness, identifying additional relevant hits.
- Predicting closed conformations for DA enzymes using AlphaFold remains challenging due to biases in training data.

# Repository Contents

- **scripts/**: Python scripts for data processing and analysis.
- **figures/**: Visualizations and structural representations of enzyme homologs.

# Future Directions

Future research should further explore the methods employed, investigate additional parameters, and conduct experimental validations to confirm potential DA enzyme candidates.

---
