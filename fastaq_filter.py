from typing import List, Union


def gc_test(seq: str, gc_bounds: tuple) -> bool:
    """
    Check sequence for GC-content.
    :param seq: given DNA sequence for analysis (str)
    :param gc_bounds: given threshold for GC-content (tuple)
    :return: True if sequence is within given threshold (bool)
    """
    gc_amount = (seq.count('G') + seq.count('C')) / len(seq) * 100
    if gc_amount < gc_bounds[0] or gc_amount > gc_bounds[1]:
        return False
    else:
        return True


def length_test(seq: str, length_bounds: tuple) -> bool:
    """
    Check sequence for length.
    :param seq: given DNA sequence for analysis (str)
    :param length_bounds: given threshold for sequence length (tuple)
    :return: True if sequence is within given threshold (bool)
    """
    if len(seq) < length_bounds[0] or len(seq) > length_bounds[1]:
        return False
    else:
        return True


def quality_test(seq_quality: str, quality_threshold: int) -> bool:
    """
    Check sequence for quality.
    :param seq_quality: string of quality (ASCII coding) for each nucleotide in given DNA sequence (str)
    :param quality_threshold: given threshold for sequence quality (int)
    :return: True if sequence is within given threshold (bool)
    """
    q_score_list = []
    for nucleotide_quality in seq_quality:
        q_score_list.append(ord(nucleotide_quality) - 33)
    if sum(q_score_list)/len(q_score_list) < quality_threshold:
        return False
    else:
        return True


def judge_seq(gc_result: bool, length_result: bool, quality_result: bool) -> bool:
    """
    Give verdict if DNA sequence fits all criteria or not
    :param gc_result: result of the GC-test (bool)
    :param length_result: result of length test (bool)
    :param quality_result: result of quality test (bool)
    :return: True if sequence fits all criteria (bool)
    """
    return gc_result and length_result and quality_result


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
        gc_result = gc_test(seqs[seq_name][0], gc_bounds)
        length_result = length_test(seqs[seq_name][0], length_bounds)
        quality_result = quality_test(seqs[seq_name][1], quality_threshold)
        if judge_seq(gc_result, length_result, quality_result):
            filtered_seqs[seq_name] = seqs[seq_name]
    return filtered_seqs
