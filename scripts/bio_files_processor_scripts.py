def get_cds_list(input_gbk: str) -> list:
    """
    Creates list of all CDS from GBK file.
    :param input_gbk: path to GBK file (str)
    :return: list of all CDSs (list)
    """
    cds_list = []
    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith("     CDS"):
                cds_list.append(line.replace(' ', '').strip('\n'))
    return cds_list


def get_gene_to_cds_dict(input_gbk: str) -> dict:
    """
    Creates dict where genes are keys and CDSs are values from GBK file.
    :param input_gbk: path to GBK file (str)
    :return: dict where genes are keys and CDSs are values (dict)
    """
    gene_dict = {}

    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith('     CDS'):
                current_cds = line.replace(' ', '').strip('\n')
            if '/gene=' in line:
                current_gene = line.replace(' ', '').strip('\n').strip('/gene=').strip('"')
                gene_dict[current_gene] = current_cds
    return gene_dict


def get_cds_translation_dict(input_gbk: str) -> dict:
    """
    Creates dict where CDSs are keys and gene names (if they are presented) and translation sequences
    are values from GBK file.
    :param input_gbk: path to GBK file (str)
    :return: dict where CDSs are keys and gene names (if they are presented) and translation sequences
    are values from GBK file (dict).
    """
    cds_dict = {}
    current_gene = ''
    current_seq = ''
    translation_start = False
    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith('     CDS'):
                translation_start = False
                if current_seq:
                    cds_dict[current_cds] = (current_gene, current_seq)
                    current_gene = ''
                    current_seq = ''
                current_cds = line.replace(' ', '').strip('\n')
            elif '/gene=' in line:
                current_gene = line.replace(' ', '').strip('\n').strip('/gene=').strip('"')
            elif '/translation' in line:
                current_seq = line.replace(' ', '').strip('\n').strip('/translation=')
                translation_start = True
            elif 'ORIGIN' in line:
                translation_start = False
            elif translation_start:
                current_seq += line.replace(' ', '').strip('\n')
    cds_dict[current_cds] = (current_gene, current_seq)
    return cds_dict
