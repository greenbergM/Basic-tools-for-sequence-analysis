from typing import List, Union
import scripts.dna_rna_tools as nucl
import scripts.protein_tools as prot
import scripts.fastaq_filter as fasq


def run_dna_rna_tools(*args: str, seq_type='DNA') -> Union[list[str], str, list[float], float]:
    """
    Launch desired operation with nucleic acid sequences. Pass comma-separated sequences,
    additional argument (if certain function requires it) and specify function name you want to apply to all sequences.
    Pass arguments strictly in this order, otherwise it won't be parsed.

    :param args:
    - nucleic acid sequences for analysis (str)
    - operation name (str): specify procedure you want to apply (str)
    :param seq_type: type of desired complement sequence (for complement function) (str)

    :return: the result of procedure in list or str format
    """
    seqs = args[:-1]
    tool = args[-1]
    seqs_identity = nucl.nucl_acid_identity(seqs)
    if tool == 'nucl_acid_identity':
        return seqs_identity
    elif tool == 'complement':
        return nucl.complement(seqs, seq_type)
    elif tool == 'reverse_complement':
        return nucl.reverse_complement(seqs, seq_type)
    elif tool == 'transcribe':
        return nucl.transcribe(seqs, seqs_identity)
    elif tool == 'reverse':
        return nucl.reverse(seqs)
    elif tool == 'gc_content':
        return nucl.gc_content(seqs)
    else:
        raise ValueError(f'There is no {tool} tool available!')


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
        if not prot.is_protein(seq, curr_encoding):
            raise ValueError('Please, use protein sequences!')
    seqs = prot.change_encoding(seqs, curr_encoding)
    if tool == 'get_seq_characteristic':
        for seq in seqs:
            processed_result.append(prot.get_seq_characteristic(seq))
    elif tool == 'find_site':
        for seq in seqs:
            processed_result.append(prot.find_site(seq, site_of_interest))
    elif tool == 'calculate_protein_mass':
        for seq in seqs:
            processed_result.append(prot.calculate_protein_mass(seq))
    elif tool == 'calculate_average_hydrophobicity':
        for seq in seqs:
            processed_result.append(prot.calculate_average_hydrophobicity(seq))
    elif tool == 'calculate_isoelectric_point':
        for seq in seqs:
            processed_result.append(prot.calculate_isoelectric_point(seq))
    elif tool == 'get_mrna':
        for seq in seqs:
            processed_result.append(prot.get_mrna(seq))
    else:
        raise ValueError(f'{tool} operation is not available!')
    return processed_result[0] if len(processed_result) == 1 else processed_result


def run_filter_fastaq(seqs: dict[str:str], gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0) -> dict[str: str]:
    """
    Filter DNA sequences based on the GC-content, length and sequencing quality (phred33).
    :param seqs: sequences to be filtered with their names and quality (dict[str:str])
    :param gc_bounds: given threshold for GC-content (tuple/int)
    :param length_bounds: given threshold for length (tuple/int)
    :param quality_threshold: given threshold for quality (int)
    :return: given dict only with sequences that fit all the criteria.
    """
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    filtered_seqs = dict()
    for seq_name in seqs.keys():
        gc_result = fasq.gc_test(seqs[seq_name][0], gc_bounds)
        length_result = fasq.length_test(seqs[seq_name][0], length_bounds)
        quality_result = fasq.quality_test(seqs[seq_name][1], quality_threshold)
        if fasq.judge_seq(gc_result, length_result, quality_result):
            filtered_seqs[seq_name] = seqs[seq_name]
    return filtered_seqs

