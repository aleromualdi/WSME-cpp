from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB

kB = 1.987e-3  # kcal/mol·K



def extract_C_alpha(structure: PDB.Structure, model_index=0, chain_id: str = None):
    """Extracts alpha-carbon (Cα) coordinates from a specific model."""
    coords = []
    residues = []

    model = structure[model_index]

    for chain in model:
        if chain_id and chain.id != chain_id:
            continue

        for residue in chain:
            if "CA" in residue:
                atom = residue["CA"]
                coords.append(atom.coord)
                residues.append(residue)

    return np.array(coords), residues


def compute_contact_map(coords: np.ndarray, cutoff: float):
    """Computes the contact map for a given set of atomic coordinates."""
    N = len(coords)
    contact_map = np.zeros((N, N))

    for i in range(N):
        for j in range(i + 2, N):
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance < cutoff:
                contact_map[i, j] = 1
                contact_map[j, i] = 1

    return contact_map


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
