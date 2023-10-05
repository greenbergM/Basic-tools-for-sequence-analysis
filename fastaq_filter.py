def filter_gc(seq: str, gc_bounds: tuple) -> bool:
    gc_amount = (seq.count('G') + seq.count('C')) / len(seq) * 100
    if gc_amount < gc_bounds[0] or gc_amount > gc_bounds[1]:
        return False
    else:
        return True


def filter_length(seq: str, length_bounds: tuple) -> bool:
    if len(seq) < length_bounds[0] or len(seq) > length_bounds[1]:
        return False
    else:
        return True


def filter_quality(seq_quality: str, quality_threshold: int) -> bool:
    q_score_list = []
    for nucleotide_quality in seq_quality:
        q_score_list.append(ord(nucleotide_quality) - 33)
    if sum(q_score_list)/len(q_score_list) < quality_threshold:
        return False
    else:
        return True


def judge_seq(gc_result: bool, length_result: bool, quality_result: bool) -> bool:
    """
    :param gc_result:
    :param length_result:
    :param quality_result:
    :return:
    """
    return gc_result and length_result and quality_result


def run_filter_fastaq(seqs: dict, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0) -> dict:
    """
    :param seqs: sequences to be filtered with their names and quality (dict)
    :param gc_bounds:
    :param length_bounds:
    :param quality_threshold:
    :return: dict
    """
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    filtered_seqs = dict()
    for seq_name in seqs.keys():
        gc_result = filter_gc(seqs[seq_name][0], gc_bounds)
        length_result = filter_length(seqs[seq_name][0], length_bounds)
        quality_result = filter_quality(seqs[seq_name][1], quality_threshold)
        if judge_seq(gc_result, length_result, quality_result):
            filtered_seqs[seq_name] = seqs[seq_name]
    return filtered_seqs
