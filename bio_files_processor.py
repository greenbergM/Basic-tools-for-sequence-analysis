import scripts.bio_files_processor_scripts as bfp


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta=None):
    """
    Creates fasta file with oneline sequences based on given fasta file with multiline sequences in the same directory.
    :param input_fasta: path to the multiline fasta file (str)
    :param output_fasta: name of oneline fasta file (str)
    """

    output_location = bfp.make_location(input_fasta, output_fasta, 'oneline_fasta', 'oneline_',
                                        '.fasta')
    seq = ''
    with open(input_fasta, mode='r') as fr, open(output_location, mode='w') as fw:
        for line in fr:
            if line.startswith('>'):
                if seq:
                    fw.write(seq + '\n')
                    seq = ''
                fw.write(line)
            else:
                seq += line.strip('\n')
        fw.write(seq)
    print('Conversion completed!')


def select_genes_from_gbk_to_fasta(*genes: str, input_gbk: str, n_before: int, n_after: int, output_fasta=None):
    """
    Creates fasta file with neighbor CDSs to given genes from GBK file and stores it
    in fasta_selected_from_gbk directory.
    :param genes: genes of interest that are used for neighbor CDSs search (str)
    :param input_gbk: path to GBK file (str)
    :param n_before: number of neighbor CDSs before gene of interest (int)
    :param n_after: number of neighbor CDSs after gene of interest (int)
    :param output_fasta: name of the output fasta file (str)
    """
    gene_list = [genes]

    # list of all CDSs in gbk file
    cds_list = []
    # dict where the keys are gene names and values are CDSs
    gene_cds_dict = {}
    # dict where the keys are CDSs and values are (gene, translation)
    cds_translation_dict = {}

    current_gene = ''
    current_seq = ''
    translation_start = False
    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith('     CDS'):
                translation_start = False
                if current_seq:
                    cds_translation_dict[current_cds] = (current_gene, current_seq)
                    current_gene = ''
                    current_seq = ''
                current_cds = line.replace(' ', '').strip('\n')
                cds_list.append(current_cds)
            elif '/gene=' in line:
                current_gene = line.replace(' ', '').strip('\n').strip('/gene=').strip('"')
                gene_cds_dict[current_gene] = current_cds
            elif '/translation' in line:
                current_seq = line.replace(' ', '').strip('\n').strip('/translation=')
                translation_start = True
            elif 'ORIGIN' in line:
                translation_start = False
            elif translation_start:
                current_seq += line.replace(' ', '').strip('\n')
    cds_translation_dict[current_cds] = (current_gene, current_seq)

    cds_of_interest = bfp.get_cds_of_interest(gene_list, gene_cds_dict, cds_list, n_before, n_after)

    output_location = bfp.make_location(input_gbk, output_fasta, 'fasta_selected_from_gbk',
                                        'CDS_selected_from_', '.fasta')

    bfp.get_fasta(output_location, cds_of_interest, cds_translation_dict)


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta=None):
    """
    Change the starting position of a DNA sequence in a FASTA file.
    :param input_fasta: Path to the input FASTA file (str).
    :param shift: The number of positions to shift the sequence (int).
    :param output_fasta: Name of the output FASTA file (str).
    """
    with open(input_fasta, mode='r') as fa:
        lines = fa.readlines()
        seq_name = lines[0]
        seq = lines[1].strip()
        shifted_seq = f'{seq[shift:]}{seq[:shift]}\n'

    output_location = bfp.make_location(input_fasta, output_fasta, 'shifted_fasta', 'shifted',
                                        '.fasta')

    with open(output_location, mode='w') as sfa:
        sfa.write(seq_name)
        sfa.write(shifted_seq)

    print('Starting position changed!')


def parse_blast_output(input_file: str, output_file=None, extension='.txt'):
    """
    Writes descriptions of best blast results from blast results file to a new file
    :param input_file: path to blast results file (str)
    :param output_file: name for output file (str)
    :param extension: extension for output file (str)
    """

    best_blast_results = sorted(bfp.get_best_blast(input_file))
    output_location = bfp.make_location(input_file, output_file, 'best_blast_results',
                                        'best_blast_', extension)

    bfp.write_blast_results(output_location, best_blast_results)
