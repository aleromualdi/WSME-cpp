from typing import List, Tuple
import matplotlib.pyplot as plt
import numpy as np
import partition_function as pf  # Import the C++ extension
from Bio import PDB

from src.utils import compute_contact_map, extract_C_alpha

kB = 1.987e-3  # kcal/mol·K

def compute_free_energy(
    Z_Q: np.ndarray, temperature: float = 310
) -> Tuple[List[float], List[float]]:
    """Computes the free-energy landscape F(Q).

    Parameters:
        Z_Q (np.ndarray): Restricted partition function Z(Q) stored as a NumPy array.
        temperature (float): Temperature in Kelvin (default = 310K).

    Returns:
        Q_values (list): Native contact fractions.
        F_values (list): Free energy per Q.
    """
    Z_total = np.sum(Z_Q)
    F_values = []
    Q_values = []

    for q in range(len(Z_Q)):
        if Z_Q[q] > 0:
            F_q = -kB * temperature * np.log(Z_Q[q] / Z_total)
            F_values.append(F_q)
            Q_values.append(q / (len(Z_Q) - 1))

    # Shift the free energy scale so the unfolded state (Q=0) is set to 0
    max_FQ = max(F_values)
    F_values = [F - max_FQ for F in F_values]

    return Q_values, F_values


# Load the protein structure
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("Protein", "7urj.pdb")

# Extract C-alpha coordinates and residues
coords, residues = extract_C_alpha(structure)

contact_map = compute_contact_map(coords, cutoff=6.0)

# Parameters for partition function calculation
entropy_penalty = -4. / 1000  # kcal/mol per residue
contact_energy = -1.0  # kcal/mol per contact
temperature = 310.0  # Kelvin

Z_Q = pf.compute_partition_function_Q_DSA(contact_map, entropy_penalty, contact_energy, temperature)

# Compute free energy based on partition function Q
Q, F = compute_free_energy(Z_Q, temperature=temperature)

plt.figure(figsize=(7, 5))
plt.plot(Q, F)
plt.xlabel("Fraction of Native Contacts (Q)")
plt.ylabel("Free Energy F(Q) (kcal/mol)")
plt.title("WSME Free Energy Landscape (DSA)")
plt.grid()
plt.show()