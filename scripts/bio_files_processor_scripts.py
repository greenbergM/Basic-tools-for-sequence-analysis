import os


def make_location(input_path: str, output_name: str, folder_name: str, commentary: str, file_type: str) -> str:
    """
    Generates location for output file and its name.
    :param input_path: path for input file (str).
    :param output_name: name for output file (str).
    :param folder_name: name for folder where output file will be saved (str)
    :param commentary: prefix for output file name (str)
    :param file_type: output file extension (str)
    :return: path to output file (str)
    """
    location = os.path.join(os.path.dirname(input_path), folder_name)
    os.makedirs(location, exist_ok=True)

    if output_name is None:
        output_name = f'{commentary}{os.path.basename(input_path)}'

    if not output_name.endswith(file_type):
        output_name = os.path.splitext(output_name)[0]
        output_name = f'{output_name}{file_type}'

    return os.path.join(location, output_name)


def get_cds_of_interest(gene_list: list, gene_cds_dict: dict, cds_list: list, n_before: int, n_after: int) -> list:
    """
    Finds neighbour CDSs for given genes.
    :param gene_list: genes of interest that are used for neighbor CDSs search (list)
    :param gene_cds_dict: dict where gene names are keys and CDSs names are values (dict)
    :param cds_list: list of CDSs (list)
    :param n_before: number of neighbor CDSs before gene of interest (int)
    :param n_after: number of neighbor CDSs after gene of interest (int)
    :return: list of chosen CDSs (list)
    """
    cds_of_interest = []
    for gene in gene_list:
        gene_cds = gene_cds_dict[gene]
        gene_position = cds_list.index(gene_cds)
        if gene_position - n_before < 0:
            raise ValueError(f'CDCs before {gene} are out of bounds! Use other value for n_before.')
        cds_before = cds_list[gene_position - n_before:gene_position]
        cds_after = cds_list[gene_position + 1:gene_position + n_after + 1]
        cds_of_interest = cds_of_interest + cds_before + cds_after
    return cds_of_interest


def get_fasta(output_fasta: str, cds_of_interest: list, translation_dict: dict):
    """
    Get FASTA file from given list of CDSs and their translations
    :param output_fasta: name for output FASTA file (str)
    :param cds_of_interest: list of cdc that is used for names
    :param translation_dict: dict where CDSs are keys and gene names (if they are presented) and translation sequences
    are values from GBK file (dict).
    """
    with open(output_fasta, mode='w') as fasta:
        for cds in cds_of_interest:
            fasta.write(f'>{cds} gene: {translation_dict[cds][0]}\n')
            fasta.write(translation_dict[cds][1].replace('"', '') + '\n')

    print('FASTA file for neighbour CDSs of given genes is created ')


def get_best_blast(input_file: str) -> list:
    """
    Gets descriptions of best blast results.
    :param input_file: path to blast results file (str)
    :return: list of descriptions of best blast results (list)
    """
    best_blast_results = []
    with open(input_file, mode='r') as blr:
        query_start = False
        description_start = False
        for line in blr:
            if line.startswith('Query #'):
                query_start = True
            elif query_start and line.startswith('Description  '):
                description_start = True
            elif query_start and description_start:
                current_result = line[:60].strip()
                best_blast_results.append(current_result)
                query_start = False
                description_start = False

    return best_blast_results
