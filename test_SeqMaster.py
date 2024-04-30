import pytest
import os

from SeqMaster import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    filter_fastq,
    IncorrectNucleotideError,
)

from bio_files_processor import convert_multiline_fasta_to_oneline


class TestBiologicalSeq:
    """
    Test cases for Biological Sequences
    """

    @pytest.fixture
    def input_dna(self):
        """
        Example DNA sequence
        """
        dna = "ACttCcGT"
        return dna

    def test_identity_dna(self, input_dna):
        """
        Test DNA sequence identity (given sequence is DNA)
        """
        inp = input_dna
        target = True
        res = DNASequence(inp).check_alphabet()
        assert target == res

    def test_identity_rna(self, input_dna):
        """
        Test RNA sequence identity (given sequence is DNA)
        """
        inp = input_dna
        with pytest.raises(IncorrectNucleotideError):
            RNASequence(inp).check_alphabet()

    def test_identity_hybrid(self):
        """
        Test DNA sequence identity (given sequence is not DNA or RNA)
        """
        inp = "ACttCcGTUUUUu"
        with pytest.raises(IncorrectNucleotideError):
            DNASequence(inp).check_alphabet()

    def test_transcribe_dna(self, input_dna):
        """
        Test DNA transcription to RNA
        """
        inp = input_dna
        target = RNASequence("ACuuCcGU").sequence
        res = DNASequence(inp).transcribe().sequence
        assert target == res

    def test_prot_len(self):
        """
        Test amino acid sequence length
        """
        inp = "AlaValSer"
        target = 3
        res = len(AminoAcidSequence(inp, 3))
        assert target == res

    def test_mol_mass(self):
        """
        Test amino acid sequence molecular mass
        """
        inp = "ARN"
        target = 89 + 174 + 132
        res = AminoAcidSequence(inp, 1).get_molecular_mass()
        assert target == res


class TestFastqFilter:
    """
    Test cases for FASTQ filtrator
    """

    @pytest.fixture
    def tmp_fastq(self):
        """
        Temporary FASTQ file
        """
        file_path = "./tmp.fastq"
        with open(file_path, "w") as f:
            # Bad seq
            f.write("@Seq4\nGATA\n+\n!!!!\n")

            # Good seq
            f.write("@Seq1\nACAT\n+\nIIII\n")

            # Too long seq
            f.write("@Seq2\nACTATA\n+\nIIIIII\n")

            # Too much GC seq
            f.write("@Seq3\nGGGA\n+\nIIII\n")

        yield file_path
        os.remove(file_path)

    def test_fastq_filter(self, tmp_fastq):
        """
        Test FASTQ file filtering
        """
        inp = tmp_fastq
        target = "@Seq1\nACAT\n+\nIIII"
        filter_fastq(inp, None, 50, 5, 10)
        with open("./fastq_filtrator_results/filtered_tmp.fastq", "r") as f:
            res = f.read().rstrip("\n")
        os.remove("./fastq_filtrator_results/filtered_tmp.fastq")
        os.rmdir("./fastq_filtrator_results")
        assert target == res


class TestConvertMultilineOneline:
    """
    Test cases for converter of multiline FASTA files to oneline
    """

    @pytest.fixture
    def tmp_fasta(self):
        """
        Temporary multiline FASTA file
        """
        file_path = "./tmp.fasta"
        with open(file_path, "w") as f:
            f.write(">Seq\nAGg\nCaT")

        yield file_path
        os.remove(file_path)

    def test_fasta_converter(self, tmp_fasta):
        """
        Test conversion of multiline FASTA to oneline
        """
        inp = tmp_fasta
        target = ">Seq\nAGgCaT"
        convert_multiline_fasta_to_oneline(inp)
        with open("./oneline_fasta/oneline_tmp.fasta", "r") as f:
            res = f.read().rstrip("\n")
        os.remove("./oneline_fasta/oneline_tmp.fasta")
        os.rmdir("./oneline_fasta")
        assert target == res
