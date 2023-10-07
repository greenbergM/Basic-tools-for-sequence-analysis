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
    Count entry of each residue type in your seq. Get description of amino acid composition.
    :param seq: protein seq in 1-letter encoding (str)
    :return: each residue type in seq in 3-letter code and its amount in current seq (dict)
    """
    seq = seq.upper()
    res_count = {}
    for res in seq:
        res_count[[tl_code for tl_code in RESIDUES_NAMES if RESIDUES_NAMES[tl_code] == res][0]] = 0
    for res in seq:
        res_count[[tl_code for tl_code in RESIDUES_NAMES if RESIDUES_NAMES[tl_code] == res][0]] += 1
    return res_count


def find_res(seq: str, res_of_interest: str) -> str:
    """
    Find all positions of certain residue in your seq
    :param seq: protein seq in 1-letter encoding (str)
    :param res_of_interest: specify the residue of interest (str)
    :return: positions of specified residue in your seq (str)
    """
    res_of_interest = res_of_interest.upper()
    seq = seq.upper()
    if len(res_of_interest) == 3:
        res_of_interest = RESIDUES_NAMES[res_of_interest]
    res_of_interest_position = []
    for ind, res in enumerate(seq, 1):
        if res == res_of_interest:
            res_of_interest_position.append(ind)
    return f'{res_of_interest} positions: {res_of_interest_position}'


def find_site(seq: str, site: str) -> str:
    """
    Find if seq contains certain site and get positions of its site
    :param seq: protein seq in 1-letter encoding (str)
    :param site: specify site of interest (str)
    :return: positions of residues for each certain site in seq (str)
    """
    site = change_residues_encoding(site).upper()
    seq = seq.upper()
    if not is_protein(site):
        return f'Site {site} is not a protein!'
    if site in seq:
        site_full_position = []
        site_count = seq.count(site)
        site_start_position = [(coordinate + 1) for coordinate in range(len(seq)) if seq.startswith(site, coordinate)]
        site_end_position = [(coordinate + len(site)) for coordinate in site_start_position]
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
    for res in seq.upper():
        total_mass += RESIDUES_CHARACTERISTICS[res][2]
    return total_mass


def calculate_average_hydrophobicity(seq: str) -> float:
    """
    Get hydrophobicity index for protein seq as sum of index for each residue in your seq divided by its length
    :param seq: protein seq in 1-letter encoding (str)
    :return: average hydrophobicity (float)
    """
    sum_hydrophobicity_ind = 0
    for res in seq.upper():
        sum_hydrophobicity_ind += RESIDUES_CHARACTERISTICS[res][0]
    return sum_hydrophobicity_ind / len(seq)


def get_mrna(seq: str) -> str:
    """
    Get encoding mRNA nucleotides for your seq
    :param seq: protein seq in 1-letter encoding (str)
    :return: potential encoding mRNA sequence with multiple choice for some positions (str)
    """
    mrna_seq = str()
    for res in seq.upper():
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
    for ind, res in enumerate(seq.upper(), 1):
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
    return pi


def analyze_secondary_structure(seq: str) -> list[str]:
    """
    Calculate the percentage of amino acids found in the three main
    types of protein secondary structure: beta-turn, beta-sheet and alpha-helix
    :param seq: protein seq in 1-letter encoding (str)
    :return: percentage of amino acids belonging to three types of secondary structure (list[str])
    """
    b_turn_set = {'G', 'P', 'N', 'D'}
    b_sheet_set = {'F', 'Y', 'I', 'V', 'C', 'W'}
    alpha_helix_set = {'M', 'A', 'L', 'E', 'K'}
    result = []
    res_for_seq = []
    count = 0
    protein_length = len(seq)
    for aa in b_turn_set:
        count += seq.upper().count(aa)
    b_turn_exp = str(count / protein_length * 100)
    res_for_seq += ['b-turn amino acids in protein' + ' is ' + b_turn_exp + '%']
    count = 0
    for aa in b_sheet_set:
        count += seq.upper().count(aa)
    b_sheet_exp = str(count / protein_length * 100)
    res_for_seq += ['b-sheet amino acids in protein' + ' is ' + b_sheet_exp + '%']
    count = 0
    for aa in alpha_helix_set:
        count += seq.upper().count(aa)
    alpha_helix_exp = str(count / protein_length * 100)
    res_for_seq += ['alpha_helix amino acids in protein' + ' is ' + alpha_helix_exp + '%']
    result += res_for_seq
    return result


def run_protein_analysis(*args: str, site_of_interest=None) -> Union[List[str], str, list[float], float]:
    """
    Launch desired operation with proteins sequences. Pass comma-separated sequences,
    additional argument (if certain function requires it) and specify function name you want to apply to all sequences.
    Pass arguments strictly in this order, otherwise it won't be parsed.

    :param args:
    - seq (str): amino acids sequences for analysis in 1-letter or 3-letter code (as many as you wish)
    - curr_encoding (int): type of encoding of given protein sequences
    - operation name (str): specify procedure you want to apply
    :param site_of_interest: one letter encoding of desired site (for find_site function)
    :return: the result of procedure in list or str format
    """
    tool = args[-1]
    curr_encoding = args[-2]
    seqs = args[:-2]
    processed_result = []
    for seq in seqs:
        if not is_protein(seq, curr_encoding):
            raise ValueError('Please, use protein sequences!')
    seqs = change_encoding(seqs, curr_encoding)
    if tool == 'get_seq_characteristic':
        for seq in seqs:
            processed_result.append(get_seq_characteristic(seq))
    elif tool == 'find_site':
        for seq in seqs:
            processed_result.append(find_site(seq, site_of_interest))
    elif tool == 'calculate_protein_mass':
        for seq in seqs:
            processed_result.append(calculate_protein_mass(seq))
    elif tool == 'calculate_average_hydrophobicity':
        for seq in seqs:
            processed_result.append(calculate_average_hydrophobicity(seq))
    elif tool == 'calculate_isoelectric_point':
        for seq in seqs:
            processed_result.append(calculate_isoelectric_point(seq))
    elif tool == 'get_seq_characteristics':
        for seq in seqs:
            processed_result.append(get_seq_characteristic(seq))
    elif tool == 'get_mrna':
        for seq in seqs:
            processed_result.append(get_mrna(seq))
    else:
        raise ValueError(f'{tool} operation is not available!')
    return processed_result[0] if len(processed_result) == 1 else processed_result
