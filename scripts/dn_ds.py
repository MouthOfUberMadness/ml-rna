#!/usr/bin/env python3

import numpy as np
import pyvolve
from compute_dnds_from_mutsel import *

constant_map_codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC",
                       "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
constant_map_codons_dict = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T", "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S", "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I", "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P", "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", "CTA": "L", "CTC": "L",
                            "CTG": "L", "CTT": "L", "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D", "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A", "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G", "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V", "TAC": "Y", "TAT": "Y", "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", "TGC": "C", "TGG": "W", "TGT": "C", "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}


def kimura(rate, alpha_1):
    kimura_star = -(0.25 + 0.25*alpha_1 + 0.25*alpha_1)
    mat = rate * np.array([[kimura_star, 0.25*alpha_1, 0.25, 0.25*alpha_1],
                           [0.25*alpha_1, kimura_star, 0.25*alpha_1, 0.25],
                           [0.25, 0.25*alpha_1, kimura_star, 0.25*alpha_1],
                           [0.25*alpha_1, 0.25, 0.25*alpha_1, kimura_star]])
    return mat


def codonsMatrix(nucleotides_matrix):
    return np.kron(nucleotides_matrix, nucleotides_matrix)


def aminoAcidToCodons(aminoAcid):
    global constant_map_codons_dict
    codons = []
    for codon, aa in constant_map_codons_dict.items():
        if (aa == aminoAcid):
            codons.append(codon)
    return codons


def codonsToIndices(codons):
    global constant_map_codons
    indices = []
    for codon in codons:
        indices.append(constant_map_codons.index(codon))
    return indices


def codonFitness(aminoAcidsWeights, syn_fitness_weight, scale):
    # codon_fitness = np.random.normal(size=61)
    codon_fitness = 0.1*np.ones(61)

    for aa, weight in aminoAcidsWeights.items():
        # print(aa)
        my_codons = aminoAcidToCodons(aa)
        # print(my_codons)
        my_codons_indices = codonsToIndices(my_codons)
        # print(my_codons_indices)

        cursor = 0
        for idx in my_codons_indices:
            # print(idx)
            if (cursor >= 0):
                codon_fitness[idx] = scale*weight
            else:
                codon_fitness[idx] = syn_fitness_weight*scale*weight
            cursor = cursor + 1

    return codon_fitness


# aminoAcids = {"A": 1.0, "R": 1.0, "N": 1.0, "D": 1.0, "C": 1.0, "Q": 1.0, "E": 1.0, "G": 1.0, "H": 1.0,
#               "I": 1.0, "L": 1.0, "K": 1.0, "M": 1.0, "F": 1.0, "P": 1.0, "S": 1.0, "T": 1.0, "W": 1.0, "Y": 1.0, "V": 1.0}
aminoAcids = {"A": 0.0, "R": 10.0, "N": 0.0, "D": 0.0, "C": 0.0, "Q": 0.0, "E": 0.0, "G": 0.0, "H": 0.0,
              "I": 0.0, "L": 0.0, "K": 0.0, "M": 0.0, "F": 0.0, "P": 0.0, "S": 0.0, "T": 0.0, "W": 0.0, "Y": 0.0, "V": 0.0}

# codon_fitness = codonFitness(
#     aminoAcidsWeights=aminoAcids, syn_fitness_weight=(1 - 0.18), scale=25.0)
# codon_fitness = codonFitness(
#     aminoAcidsWeights=aminoAcids, syn_fitness_weight=(1 - 0.18), scale=35.0)
codon_fitness = codonFitness(
    aminoAcidsWeights=aminoAcids, syn_fitness_weight=(1 - 0.18), scale=1.0)


print(codon_fitness)

mutsel_codon_model_fits = pyvolve.Model("MutSel", {"fitness": codon_fitness})
mutsel_codon_frequencies = np.array(
    mutsel_codon_model_fits.params["state_freqs"])

print(mutsel_codon_frequencies)

c = dNdS_from_MutSel(mutsel_codon_frequencies)
dnds = c.compute_dnds()  # to return dN/dS
dn = c.compute_dn()   # to return dN
ds = c.compute_ds()   # to return dS
print("dnds = ", dnds)
print("dn = ", dn)
print("ds = ", ds)


my_partition = pyvolve.Partition(
    models=mutsel_codon_model_fits, root_sequence="GATAGAACT")
#my_tree = pyvolve.read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762_m1_):0.921_m2_):0.207);")
my_tree = pyvolve.read_tree(tree="(r:0.207);")
# pyvolve.print_tree(my_tree)

my_evolver = pyvolve.Evolver(tree=my_tree, partitions=my_partition)
my_evolver()
