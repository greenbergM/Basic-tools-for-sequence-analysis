from typing import List, Union

# 3-letter with corresponding 1-letter residues names
RESIDUES_NAMES = {'ALA': 'A',
                  'ARG': 'R',
                  'ASN': 'N',
                  'ASP': 'D',
                  'CYS': 'C',
                  'GLN': 'Q',
                  'GLU': 'E',
                  'GLY': 'G',
                  'HIS': 'H',
                  'ILE': 'I',
                  'LEU': 'L',
                  'LYS': 'K',
                  'MET': 'M',
                  'PHE': 'F',
                  'PRO': 'P',
                  'SER': 'S',
                  'THR': 'T',
                  'TRP': 'W',
                  'TYR': 'Y',
                  'VAL': 'V'
                  }

# 1-letter with corresponding 3-letter residues names
RESIDUES_NAMES_REVERSED = {'A': 'ALA',
                           'R': 'ARG',
                           'N': 'ASN',
                           'D': 'ASP',
                           'C': 'CYS',
                           'Q': 'GLN',
                           'E': 'GLU',
                           'G': 'GLY',
                           'H': 'HIS',
                           'I': 'ILE',
                           'L': 'LEU',
                           'K': 'LYS',
                           'M': 'MET',
                           'F': 'PHE',
                           'P': 'PRO',
                           'S': 'SER',
                           'T': 'THR',
                           'W': 'TRP',
                           'Y': 'TYR',
                           'V': 'VAL'
                           }

# first value is hydrophobicity index, second is pKa (pKa1, pKa2, pKa3 respectively), third is molecular mass in Da
RESIDUES_CHARACTERISTICS = {'A': [1.8, [2.34, 9.69, 0], 89],
                            'R': [-4.5, [2.17, 9.04, 12.48], 174],
                            'N': [-3.5, [2.02, 8.80, 0], 132],
                            'D': [-3.5, [1.88, 9.60, 3.65], 133],
                            'C': [2.5, [1.96, 10.28, 8.18], 121],
                            'Q': [-3.5, [2.17, 9.13, 0], 146],
                            'E': [-3.5, [2.19, 9.67, 4.25], 147],
                            'G': [-0.4, [2.34, 9.60, 0], 75],
                            'H': [-3.2, [1.82, 9.17, 6.00], 155],
                            'I': [4.5, [2.36, 9.60, 0], 131],
                            'L': [3.8, [2.36, 9.60, 0], 131],
                            'K': [-3.9, [2.18, 8.95, 10.53], 146],
                            'M': [1.9, [2.28, 9.21, 0], 149],
                            'F': [2.8, [1.83, 9.13, 0], 165],
                            'P': [-1.6, [1.99, 10.60, 0], 115],
                            'S': [-0.8, [2.21, 9.15, 0], 105],
                            'T': [-0.7, [2.09, 9.10, 0], 119],
                            'W': [-0.9, [2.83, 9.39, 0], 204],
                            'Y': [-1.3, [2.20, 9.11, 0], 181],
                            'V': [4.2, [2.32, 9.62, 0], 117]}

# amino acid with corresponding degenerate codon/codons
AMINO_ACID_TO_MRNA = {'A': 'GCN',
                      'R': '(CGN/AGR)',
                      'N': 'AAY',
                      'D': 'GAY',
                      'C': 'UGY',
                      'Q': 'CAR',
                      'E': 'GAR',
                      'G': 'GGN',
                      'H': 'CAY',
                      'I': 'AUH',
                      'L': '(CUN/UUR)',
                      'K': 'AAR',
                      'M': 'AUG',
                      'F': 'UUY',
                      'P': 'CCN',
                      'S': '(UCN/AGY)',
                      'T': 'ACN',
                      'W': 'UGG',
                      'Y': 'UAY',
                      'V': 'GUN'}


def make_three_letter(seq: str):
    """
    Get amino acid list from protein in 3-letter encoding
    :param seq: protein seq in 3-letter encoding (str)
    :return: list of amino acids (list)
    """
    seq = seq.replace(" ", "")
    three_letter_seq_list = [seq[counter:counter + 3].upper() for counter in range(0, len(seq), 3)]
    return three_letter_seq_list


def is_protein(seq: str, curr_encoding: int) -> bool:
    """
    Check if sequence is protein or not by identify invalid seq elements, which are not presented in dicts above.
    :param seq: protein seq in 3-letter or 1-letter encoding (str)
    :param curr_encoding: protein encoding that is used in input (int)
    :return: if seq is correct protein seq or not (bool)
    """
    if curr_encoding == 1:
        if set(seq.upper()).difference((RESIDUES_NAMES.values())):
            return False
        else:
            return True
    elif curr_encoding == 3:
        amino_acids = make_three_letter(seq)
        if set(amino_acids).difference((RESIDUES_NAMES.keys())):
            return False
        else:
            return True
    else:
        raise ValueError(f'{curr_encoding}-letter amino acid encoding is not available! Use 1 or 3.')


def change_encoding(seqs: tuple[str], curr_encoding: int) -> list[str]:
    """
    Get protein sequence in 1-letter encoding
    :param seqs: protein seq in 3-letter or 1-letter encoding (str)
    :param curr_encoding: protein encoding that is used in input (int)
    :return: protein sequence in 1-letter encoding (str)
    """
    renamed_seqs = []
    for seq in seqs:
        renamed_seq = str()
        if curr_encoding == 3:
            amino_acids = make_three_letter(seq)
            for amino_acid in amino_acids:
                renamed_seq += RESIDUES_NAMES[amino_acid]
            renamed_seqs.append(renamed_seq)
        else:
            renamed_seqs.append(seq.upper())
    return renamed_seqs


def get_seq_characteristic(seq: str) -> dict:
    """
    Count entry of each residue type in your seq. Get description of amino acid composition in dict format.
    :param seq: protein seq in 1-letter encoding (str)
    :return: each residue type in seq in 3-letter code and its amount in current seq (dict)
    """
    residue_count = {}
    for residue in set(seq):
        residue_entry = seq.count(residue)
        residue_count[RESIDUES_NAMES_REVERSED[residue]] = residue_entry
    return residue_count


def find_site(seq: str, site: str) -> str:
    """
    Find if seq contains certain site and get positions of this site
    :param seq: protein seq in 1-letter encoding (str)
    :param site: specify site of interest (str)
    :return: positions of residues for each certain site in seq (str)
    """
    site = site.upper()
    if set(site).difference(RESIDUES_NAMES.values()):
        return f'Site {site} is not a protein!'
    if site in seq:
        site_full_position = []
        site_count = seq.count(site)
        site_start_position = [(coordinate + 1) for coordinate in range(len(seq)) if seq.startswith(site, coordinate)]
        site_end_position = [(coordinate + len(site) - 1) for coordinate in site_start_position]
        for counter in range(len(site_start_position)):
            site_full_position.append(f'{site_start_position[counter]}:{site_end_position[counter]}')
        return f'Site entry in sequence = {site_count}. Site residues can be found at positions: {site_full_position}'
    else:
        return f'{site} site is not in sequence!'


def calculate_protein_mass(seq: str) -> float:
    """
    Get mass of residues in your seq in Da
    :param seq: protein seq in 1-letter encoding (str)
    :return: mass in Da (float)
    """
    total_mass = 0
    for res in seq:
        total_mass += RESIDUES_CHARACTERISTICS[res][2]
    return total_mass


def calculate_average_hydrophobicity(seq: str) -> float:
    """
    Get hydrophobicity index for protein seq as sum of index for each residue in your seq divided by its length
    :param seq: protein seq in 1-letter encoding (str)
    :return: average hydrophobicity (float)
    """
    sum_hydrophobicity_ind = 0
    for res in seq:
        sum_hydrophobicity_ind += RESIDUES_CHARACTERISTICS[res][0]
    return round(sum_hydrophobicity_ind / len(seq), 2)


def get_mrna(seq: str) -> str:
    """
    Get encoding mRNA nucleotides for your seq
    :param seq: protein seq in 1-letter encoding (str)
    :return: potential encoding mRNA sequence with multiple choice for some positions (str)
    """
    mrna_seq = str()
    for res in seq:
        mrna_seq += AMINO_ACID_TO_MRNA[res]
    return mrna_seq


def calculate_isoelectric_point(seq: str) -> float:
    """
    Find isoelectrinc point as sum of known pI for residues in your seq
    :param seq: protein seq in 1-letter encoding (str)
    :return: isoelectric point (float)
    """
    sum_pka = 0
    pka_amount = 0
    for ind, res in enumerate(seq, 1):
        if ind == 1:
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][1]
            pka_amount += 1
        elif RESIDUES_CHARACTERISTICS[res][1][2] != 0:
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][2]
            pka_amount += 1
        elif ind == len(seq):
            sum_pka += RESIDUES_CHARACTERISTICS[res][1][0]
            pka_amount += 1
    pi = sum_pka / pka_amount
    return round(pi, 2)
