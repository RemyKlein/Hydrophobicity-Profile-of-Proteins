import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from io import StringIO
import argparse
import requests
import os

# DataFrame with hydrophobicity scale
hydrophobicity_scale = pd.DataFrame({
    "Kyte-Doolittle": {"A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3},
    "Hopp-Woods": {"A": -0.5, "C": -1.0, "D": 3.0, "E": 3.0, "F": -2.5, "G": 0.0, "H": -0.5, "I": -1.8, "K": 3.0, "L": -1.8, "M": -1.3, "N": 0.2, "P": 0.0, "Q": 0.2, "R": 3.0, "S": 0.3, "T": -0.4, "V": -1.5, "W": -3.4, "Y": -2.3},
    "Cornette": {"A": 0.2, "C": 4.1, "D": -3.1, "E": -1.8, "F": 4.4, "G": 0.0, "H": 0.5, "I": 4.8, "K": -3.1, "L": 5.7, "M": 4.2, "N": -0.5, "P": -2.2, "Q": -2.8, "R": 1.4, "S": -0.5, "T": -1.9, "V": 4.7, "W": 1.0, "Y": 3.2},
    "Eisenberg": {"A": 0.62, "C": 0.29, "D": -0.9, "E": -0.74, "F": 1.19, "G": 0.48, "H": -0.4, "I": 1.38, "K": -1.5, "L": 1.06, "M": 0.64, "N": -0.78, "P": 0.12, "Q": -0.85, "R": -2.53, "S": -0.18, "T": -0.05, "V": 1.08, "W": 0.81, "Y": 0.26},
    "Rose": {"A": 0.74, "C": 0.91, "D": 0.62, "E": 0.62, "F": 0.88, "G": 0.72, "H": 0.78, "I": 0.88, "K": 0.52, "L": 0.85, "M": 0.85, "N": 0.63, "P": 0.64, "Q": 0.62, "R": 0.64, "S": 0.66, "T": 0.70, "V": 0.86, "W": 0.85, "Y": 0.76},
    "Janin": {"A": 0.3, "C": 0.9, "D": -0.6, "E": -0.7, "F": 0.5, "G": 0.3, "H": -0.1, "I": 0.7, "K": -1.8, "L": 0.5, "M": 0.4, "N": -0.5, "P": -0.3, "Q": -0.7, "R": -1.4, "S": -0.1, "T": -0.2, "V": 0.6, "W": 0.3, "Y": -0.4},
    "Engelman GES": {"A": 1.6, "C": 2.0, "D": -9.2, "E": -8.2, "F": 3.7, "G": 1.0, "H": -3.0, "I": 3.1, "K": -8.8, "L": 2.8, "M": 3.4, "N": -4.8, "P": -0.2, "Q": -4.1, "R": -12.3, "S": 0.6, "T": 1.2, "V": 2.6, "W": 1.9, "Y": -0.7}
})

hydrophobicity_scale.index.name = "Amino acids"

methods = {
    1: "Kyte-Doolittle",
    2: "Hopp-Woods",
    3: "Cornette",
    4: "Eisenberg",
    5: "Rose",
    6: "Janin",
    7: "Engelman GES" 
}

def fetch_uniprot_sequence(uniprot_id):
    """
    Fetch the protein sequence from UniProt using its REST API.

    Parameters:
        uniprot_id (str): UniProt identifier of the protein.

    Returns:
        str: Amino acid sequence of the protein.

    Raises:
        ValueError: If the UniProt ID cannot be fetched.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError(f"Could not fetch UniProt ID {uniprot_id}")
    fast_io = StringIO(response.text) # Creation temporary folder to be read by SeqIO.
    record = SeqIO.read(fast_io, "fasta") # Reading of FASTA protein file.
    return str(record.seq)

def compute_hydrophobicity(sequence, method, window_size):
    """
    Compute hydrophobicity scores for a protein sequence using a given scale and sliding window.

    Parameters:
        sequence (str): Amino acid sequence of the protein.
        method (pd.Series): Hydrophobicity scale for each amino acid.
        window_size (int): Size of the sliding window used to calculate the average score.

    Returns:
        tuple:
            - scores (list of float): Hydrophobicity score for each window along the sequence.
            - scores_hydro (list of int): Positions (indices) where the hydrophobicity score > 0.
    """
    scores = []
    scores_hydro = []
    for i in range(len(sequence) - window_size + 1):
        subseq = sequence[i:i+window_size]
        score = np.mean([method.loc[aa] for aa in subseq])
        scores.append(score)
        if score > 0:
            scores_hydro.append(i)
    return scores, scores_hydro

def plot_hydrophobicity(sequence, method, window_size, threshold, output_file):
    """
    Generate and save a hydrophobicity profile plot for a protein sequence.

    Parameters:
        sequence (str): Protein sequence.
        method (pd.Series): Hydrophobicity scale for amino acids.
        window_size (int): Sliding window size for averaging scores.
        threshold (float): Hydrophobicity threshold to draw on the plot.
        output_file (str): Path to save the plot image.

    Returns:
        None. The plot is saved to the specified file.
    """
    scores, _ = compute_hydrophobicity(sequence, method, window_size) 
    positions = range(len(scores))
    plt.plot(positions, scores, color="blue", label="Hydrophobicity")
    plt.axhline(y=threshold, color="k", linestyle='--', label="Hydrophilic/Hydrophobic threshold")
    plt.scatter(positions, scores, color="blue", s=20, alpha=0.4)
    plt.title("Hydrophobicity profile", fontweight="bold")
    plt.xlabel("Position in sequence")
    plt.ylabel("Hydrophobicity score")
    plt.tight_layout()
    plt.legend()
    plt.savefig(output_file)
    plt.show()

def save_positions_to_file(sequence, method, window_size, filename="hydrophobic_positions.txt"):
    """
    Save positions of hydrophobic amino acids in a protein sequence to a text file.

    Parameters:
        sequence (str): Protein sequence.
        method (pd.Series): Hydrophobicity scale for amino acids.
        window_size (int): Sliding window size used to calculate scores.
        filename (str): Output file path. Default is "hydrophobic_positions.txt".

    Returns:
        None. The file is written with positions where the hydrophobicity score > 0.

    Notes:
        - Each position corresponds to the index of the first amino acid in the window 
          whose average hydrophobicity score is positive.
    """
    _, positions = compute_hydrophobicity(sequence, method, window_size)
    with open(filename, "w") as file:
        file.write("Positions of hydrophobic amino acids:\n")
        for pos in positions:
            file.write(f"{pos}\n")

def main():
    parser = argparse.ArgumentParser(description="Compute protein hydrophobicity profiles from UniProt sequences.")
    parser.add_argument("ids", nargs="+", help="UniProt protein IDs (e.g. P69905 P68871)")
    parser.add_argument("--method", choices=hydrophobicity_scale.columns, default="Kyte-Doolittle", help="Hydrophobicity scale")
    parser.add_argument("--window", type=int, default=11, help="Sliding window size")
    parser.add_argument("--threshold", type=float, default=0.0, help="Hydrophobicity threshold")
    parser.add_argument("--outdir", default="results", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for uid in args.ids:
        print(f"Processing {uid}...")
        seq = fetch_uniprot_sequence(uid)
        scores, hydro_pos = compute_hydrophobicity(seq, hydrophobicity_scale[args.method], args.window)

        # Utiliser la fonction pour sauvegarder les positions
        pos_file = os.path.join(args.outdir, f"{uid}_hydro_positions.txt")
        save_positions_to_file(seq, hydrophobicity_scale[args.method], args.window, filename=pos_file)
        
        plot_file = os.path.join(args.outdir, f"{uid}_hydro_plot.png")
        plot_hydrophobicity(seq, hydrophobicity_scale[args.method], args.window, threshold=args.threshold, output_file=plot_file)

if __name__ == "__main__":
    main()
