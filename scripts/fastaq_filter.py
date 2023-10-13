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


def get_dict(fastaq_file_path: str) -> dict[str:tuple[str, str]]:
    """
    Converts fastaq file into dict where key is sequence name, and seqs, quality are value
    :param fastaq_file_path: path to the fastaq file (str)
    :return: dict where key is sequence name, and seqs, quality are value (dict[str:tuple[str, str])
    """
    with open(fastaq_file_path, mode='r') as fastaq:
        fastaq_dict = {}
        counter = 1
        for line in fastaq:
            if counter == 1:
                name = line.strip('\n')
                counter += 1
            elif counter == 2:
                seq = line.strip('\n')
                counter += 1
            elif counter == 3:
                commentary = line.strip('\n')
                counter += 1
            elif counter == 4:
                quality = line.strip('\n')
                fastaq_dict[name] = (seq, commentary, quality)
                counter = 1
    return fastaq_dict






