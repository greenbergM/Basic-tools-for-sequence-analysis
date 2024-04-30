from dataclasses import dataclass
from typing import Optional
from additional_scripts.bio_files_processor_scripts import (
    make_location,
    get_fasta,
    get_best_blast,
    get_cds_of_interest,
    write_blast_results,
)


### HW15


@dataclass
class FastaRecord:
    """
    Dataclass for storing FASTA record (ID, description and sequence).
    """

    id: str
    description: str
    seq: str

    def __repr__(self):
        return f">{self.id} {self.description}\n{self.seq}"


class FastaIter:
    """
    Iterator for FASTA files. Returns records in FastaRecord dataclass format.
    Allows for multiline FASTA records.
    """

    def __init__(self, fasta_handler):
        self._fasta = fasta_handler
        self._records = []
        self.current_pos = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_pos == len(self._records):
            seq = ""
            id_and_description = None
            for line in self._fasta:
                line = line.strip()
                if line.startswith(">"):
                    if seq:
                        self._records.append(
                            (id_and_description[0], id_and_description[1], seq)
                        )
                        seq = ""
                    id_and_description = line[1:].split(maxsplit=1)
                else:
                    seq += line
            if seq:
                self._records.append(
                    (id_and_description[0], id_and_description[1], seq)
                )

            if self.current_pos == len(self._records):
                raise StopIteration

        current_pos = self.current_pos
        self.current_pos += 1
        return FastaRecord(*self._records[current_pos])


class OpenFasta:
    """
    Context manager for opening and reading a FASTA file.
    Returns records in FastaRecord dataclass format.
    Allows for multiline FASTA records.
    """

    def __init__(self, file_path, mode="r"):
        self.file_path = file_path
        self.mode = mode
        self.iterator = None

    def __enter__(self):
        self.file_handler = open(self.file_path, self.mode)
        self.iterator = FastaIter(self.file_handler)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_handler:
            self.file_handler.close()

    def __iter__(self):
        return self.iterator

    def read_record(self):
        """
        Return the next FastaRecord from the FASTA file.
        """

        return next(iter(self))

    def read_records(self):
        """
        Return all FastaRecords from the FASTA file.
        """
        all_records = []
        for record in self.iterator:
            all_records.append(record)
        return all_records


## HW6


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: Optional[str] = None
) -> None:
    """
    Creates FASTA file with oneline sequences based on given FASTA file with multiline sequences and stores it
    in oneline_fasta folder in the same directory.

    :param input_fasta: path to the multiline FASTA file (str)
    :param output_fasta: name of oneline FASTA file (str);
    if not given the output file name will be oneline_*input file name*.fasta
    """

    output_location = make_location(
        input_fasta, output_fasta, "oneline_fasta", "oneline_", ".fasta"
    )

    with open(input_fasta) as inp_fa, open(output_location, "w") as opt_fa:
        first_line = True
        for line in inp_fa:
            if line.startswith(">"):
                if not first_line:
                    opt_fa.write("\n")
                opt_fa.write(line)
                first_line = False
            else:
                opt_fa.write(line.strip())

    print("Conversion completed!")


def select_genes_from_gbk_to_fasta(
    *genes: str,
    input_gbk: str,
    n_before: int,
    n_after: int,
    output_fasta: Optional[str] = None,
) -> None:
    """
    Creates FASTA file with neighbor CDSs to given genes from GBK file and stores it
    in fasta_selected_from_gbk folder in the same directory.

    :param genes: genes of interest that are used for neighbor CDSs search (str)
    :param input_gbk: path to GBK file (str)
    :param n_before: number of neighbor CDSs before gene of interest (int)
    :param n_after: number of neighbor CDSs after gene of interest (int)
    :param output_fasta: name of the output FASTA file (str);
    if not given the output file name will be CDS_selected_from_*input file name*.fasta
    """

    gene_list = [genes]

    # list of all CDSs in gbk file
    cds_list = []

    # dict where the keys are gene names and values are CDSs
    gene_cds_dict = {}

    # dict where the keys are CDSs and values are (gene, translation)
    cds_translation_dict = {}

    current_gene = ""
    current_seq = ""
    translation_start = False
    with open(input_gbk, mode="r") as gbk:
        for line in gbk:
            if line.startswith("     CDS"):
                translation_start = False
                if current_seq:
                    cds_translation_dict[current_cds] = (current_gene, current_seq)
                    current_gene = ""
                    current_seq = ""
                current_cds = line.replace(" ", "").strip("\n")
                cds_list.append(current_cds)
            elif "/gene=" in line:
                current_gene = (
                    line.replace(" ", "").strip("\n").strip("/gene=").strip('"')
                )
                gene_cds_dict[current_gene] = current_cds
            elif "/translation" in line:
                current_seq = line.replace(" ", "").strip("\n").strip("/translation=")
                translation_start = True
            elif "ORIGIN" in line:
                translation_start = False
            elif translation_start:
                current_seq += line.replace(" ", "").strip("\n")
    cds_translation_dict[current_cds] = (current_gene, current_seq)

    cds_of_interest = get_cds_of_interest(
        gene_list, gene_cds_dict, cds_list, n_before, n_after
    )

    output_location = make_location(
        input_gbk,
        output_fasta,
        "fasta_selected_from_gbk",
        "CDS_selected_from_",
        ".fasta",
    )

    get_fasta(output_location, cds_of_interest, cds_translation_dict)


def change_fasta_start_pos(
    input_fasta: str, shift: int, output_fasta: Optional[str] = None
) -> None:
    """
    Change the starting position of a DNA sequence in a FASTA file and stores it
    in shifted_fasta folder in the same directory.

    :param input_fasta: path to the input FASTA file (str).
    :param shift: the number of positions to shift the sequence (int).
    :param output_fasta: name of the output FASTA file (str);
    if not given the output file name will be shifted_*input file name*.fasta
    """
    with open(input_fasta, mode="r") as fa:
        lines = fa.readlines()
        seq_name = lines[0]
        seq = lines[1].strip()
        shifted_seq = f"{seq[shift:]}{seq[:shift]}\n"

    output_location = make_location(
        input_fasta, output_fasta, "shifted_fasta", "shifted_", ".fasta"
    )

    with open(output_location, mode="w") as sfa:
        sfa.write(seq_name)
        sfa.write(shifted_seq)

    print("Starting position changed!")


def parse_blast_output(
    input_file: str, output_file: Optional[str] = None, extension: str = ".txt"
) -> None:
    """
    Writes descriptions of best blast results from blast results file to a new file and stores it
    in best_blast_results folder in the same directory.

    :param input_file: path to blast results file (str)
    :param output_file: name for output file (str);
    if not given the output file name will be best_blast_*input file name*.*extension*
    :param extension: extension for output file (str); if not changed will be '.txt'
    """

    best_blast_results = sorted(get_best_blast(input_file))
    output_location = make_location(
        input_file, output_file, "best_blast_results", "best_blast_", extension
    )

    write_blast_results(output_location, best_blast_results)
