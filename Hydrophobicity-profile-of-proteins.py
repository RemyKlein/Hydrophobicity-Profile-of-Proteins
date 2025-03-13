import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, Entrez

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
# Name of the index (not mandatory)
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

def load_protein_sequence(id_uniprot):
    """
    This function takes the id of the protein from UniProt database and return the amino acids sequence.

    Parameters:
    id_uniprot (str): ID of the protein from UniProt

    Returns:
    str: The amino acids sequence
    """
    # Email to access database UniProt
    Entrez.email = input("Enter your email to access the UniProt database: ")
    with Entrez.efetch(db="protein", id=id_uniprot, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
        return str(record.seq)

def compute_hydrophobicity(sequence, method, window_size):
    """
    This function computes the hydrophobicity score for a protein sequence using a selected method and window size. 
    It returns the hydrophobicity scores for each amino acid in the sequence, and the positions where amino acids 
    are considered hydrophobic based on the chosen method.

    Parameters:
    sequence (str): The sequence of amino acids of the protein selected.
    method (DataFrame): The hydrophobicity profile method chosen (e.g., Kyte-Doolittle, Hopp-Woods, etc.).
    window_size (int): The size of the sliding window used to calculate the hydrophobicity score.

    Returns:
    scores (list): A list containing the hydrophobicity score for each window of amino acids in the sequence.
    scores_hydro (list): A list of positions where the hydrophobicity score exceeds the threshold (i.e., positive score).
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

def plot_hydrophobicity():
    """
    This function generates and displays a plot of the hydrophobicity profile of a protein sequence. 
    It uses the hydrophobicity scores calculated for the sequence, based on a selected method and 
    window size. The plot visualizes the hydrophobicity scores along the sequence, marking the 
    hydrophilic/hydrophobic threshold with a dashed line and highlighting the positions of amino 
    acids with varying hydrophobicity.

    Parameters:
    None (this function relies on external variables such as "sequence", "choice_hydrophobicity_scale", 
    and "window_size_choice" for its operation).

    Returns:
    None (The function directly shows the plot).

    Notes:
    - The plot will display the hydrophobicity scores as a line graph with scatter points.
    - The threshold for hydrophobicity (score > 0) is marked by a horizontal dashed line at y=0.
    """
    scores, _ = compute_hydrophobicity(sequence, choice_hydrophobicity_scale, window_size_choice) 
    positions = range(len(scores))
    plt.plot(positions, scores, color="blue", label="Hydrophobicity")
    plt.axhline(y=0, color="k", linestyle='--', label="Hydrophilic/Hydrophobic threshold")
    plt.scatter(positions, scores, color="blue", s=20, alpha=0.4)
    plt.title("Hydrophobicity profile", fontweight="bold")
    plt.xlabel("Position in sequence")
    plt.ylabel("Hydrophobicity score")
    plt.tight_layout()
    plt.legend()
    plt.savefig(f"hydrophobic_profile.pdf")
    plt.show()

def save_positions_to_file(filename="hydrophobic_positions.txt"):
    """
    This function saves the positions of hydrophobic amino acids in the protein sequence to a text file. 
    The positions are determined based on the hydrophobicity profile calculated using a selected method 
    and a specified window size. Only positions with a positive hydrophobicity score are saved.

    Parameters:
    filename (str): The name of the file where the positions will be saved. Default is "hydrophobic_positions.txt".

    Returns:
    None (The function writes the positions to the file and does not return any value).
    
    Notes:
    - Each position is written on a new line in the text file.
    - The positions correspond to the indices in the protein sequence where the hydrophobicity score exceeds the threshold (positive score).
    """
    _, positions = compute_hydrophobicity(sequence, choice_hydrophobicity_scale, window_size_choice)
    with open(filename, "w") as file:
        file.write("Positions of hydrophobic amino acids:\n")
        for pos in positions:
            file.write(f"{pos}\n")

# Loading of protein ID
id_uniprot = input("Enter the UniProt ID of your protein: ")

# Loading of sliding window size
print("Select a size for the size of the sliding window (recommended size: 7 or 11)")
while True:
    try:
        window_size_choice = int(input("Your window size choice: "))
        if window_size_choice <= 0:
            print("Error: The window size must be a positive integer greater than zero!")
        else:
            break
    except ValueError:
        print("Error: You must enter a valid number!")

# Loading the method for the hydrophobicity profile
print("Select a numeric value for the choice of method:\n1 - Kyte-Doolittle\n2 - Hopp-Woods\n3 - Cornette\n4 - Eisenberg\n5 - Rose\n6 - Janin\n7 - Engelman GES")
while True:
    try:
        choice_user = int(input("Your choice (1 - 7): "))
        if 1 <= choice_user <= 7:
            break
        else:
            print("Error: The number must be between 1 and 7.")
    except ValueError:
        print("Error: You must enter a valid number!")
choice_method = methods[choice_user]
choice_hydrophobicity_scale = hydrophobicity_scale[choice_method]

sequence = load_protein_sequence(id_uniprot)
plot_hydrophobicity()
save_positions_to_file()