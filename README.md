# Hydrophobicity Profile of Proteins
 This project computes the hydrophobicity profile of a given protein sequence using multiple well-known scales: Kyte-Doolittle, Hopp-Woods, Cornette, Eisenberg, Rose, Janin, and Engelman GES. The tool uses a sliding window approach to calculate the average hydrophobicity score and visualizes the hydrophobic and hydrophilic regions of the protein sequence. The user can choose the method and window size, then analyze and save the results.

 ## **Features:**
- Fetch protein sequences directly from the UniProt database using the protein's UniProt ID.
- Compute hydrophobicity scores for the sequence based on selected methods.
- Display the hydrophobicity profile as a plot, highlighting hydrophobic (positive score) regions.
- Save the positions of hydrophobic amino acids in a text file for further analysis.

## **How It Works:**
1. **Fetch Protein Sequence**: The user inputs the UniProt ID of the protein of interest. The script then fetches the amino acid sequence from the UniProt database.
2. **Select Hydrophobicity Method**: The user selects one of the popular hydrophobicity scales, including Kyte-Doolittle, Hopp-Woods, Cornette, Eisenberg, Rose, Janin, or Engelman GES.
3. **Sliding Window Calculation**: The user specifies a window size (recommended size: 7 or 11). The script calculates the average hydrophobicity score within the selected window, then slides across the entire sequence.
4. **Hydrophobicity Plot**: The hydrophobicity scores are plotted along the sequence. The plot visually represents the hydrophobic and hydrophilic regions.
5. **Save Results**: The positions of hydrophobic amino acids (with a score greater than 0) are saved in a text file for further review.

## **Example:**
- **Input**: 
    - UniProt ID: `P11229`
    - Hydrophobicity scale: Kyte-Doolittle
    - Sliding window size: 9
- **Output**: 
    - A graph showing the hydrophobicity profile of the protein sequence.
    - A text file (`hydrophobic_positions.txt`) containing the positions of hydrophobic amino acids.
    - A PDF file (`hydrophobic_profile.pdf`) containing the graph of hydrophobic profile of the protein.

## **Requirements:**
- Python 3.x
- `pandas`
- `numpy`
- `matplotlib`
- `biopython`

## **Installation:**
Clone the repository and install the necessary dependencies:

```bash
git clone https://github.com/RemyKlein/hydrophobicity-profile-proteins.git
cd hydrophobicity-profile-proteins
pip install -r requirements.txt