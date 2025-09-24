# Hydrophobicity Profile of Proteins
 This project computes the **hydrophobicity profile** of a protein sequence using multiple well-known scales: Kyte-Doolittle, Hopp-Woods, Cornette, Eisenberg, Rose, Janin, and Engelman GES. The tool uses a **sliding window** approach to calculate the average hydrophobicity score and visualizes the hydrophobic and hydrophilic regions. Users can choose the hydrophobicity method and window size, then save and analyze results.

 ## **Features**
- Fetch protein sequences directly from the **UniProt database** using UniProt IDs.
- Compute hydrophobicity scores using selected scales.
- Generate a plot highlighting **hydrophobic (positive score)** regions.
- Save the positions of hydrophobic amino acids to a text file.
- Save the hydrophobicity profile plot as an image (.png) and optionally as a PDF.

## **How It Works**
1. **Fetch Protein Sequence**: Input a UniProt ID (e.g., `P11229`) and fetch the amino acid sequence from UniProt.
2. **Select Hydrophobicity Method**: Choose one of the popular hydrophobicity scales (Kyte-Doolittle, Hopp-Woods, Cornette, Eisenberg, Rose, Janin, Engelman GES).
3. **Sliding Window Calculation**: Specify a window size (recommended 7 or 11). The script computes the average hydrophobicity score for each window along the sequence.
4. **Hydrophobicity Plot**: The computed scores are plotted along the sequence, showing hydrophobic and hydrophilic regions. A horizontal line indicates the threshold.
5. **Save Results**: 
- **Positions file**: Text file listing indices of hydrophobic amino acids (score > threshold).
- **Plot file**: PNG (or PDF) file with the hydrophobicity profile.

## **Example Usage**
- **Input**: 
    - UniProt ID: `P11229`
    - Hydrophobicity scale: Kyte-Doolittle
    - Sliding window size: 9
- **Output**: 
    - `results/P11229_hydro_positions.txt`: positions of hydrophobic amino acids.
    - `results/P11229_hydro_plot.png`: hydrophobicity profile plot.
    - Optional: PDF version of the plot.

## **Requirements:**
- Python 3.x
- `pandas`
- `numpy`
- `matplotlib`
- `biopython`
- `requests`

## **Installation:**
Clone the repository and install the necessary dependencies:

```bash
git clone https://github.com/RemyKlein/Hydrophobicity-Profile-of-Proteins.git
cd Hydrophobicity-Profile-of-Proteins
pip install -r requirements.txt
```

## **Running the Script**
```bash
python Hydrophobicity-profile-of-proteins.py P11229 --method Kyte-Doolittle --window 9 --threshold 0 --outdir results
```
- `--method`: Hydrophobicity scale (default: Kyte-Doolittle).
- `--window`: Sliding window size (default: 11).
- `--threshold`: Score threshold for hydrophobic positions (default: 0).
- `--outdir`: Directory to save outputs (default: results).