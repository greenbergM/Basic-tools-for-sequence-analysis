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
    @pytest.fixture
    def input_dna(self):
        dna = "ACttCcGT"
        return dna

    def test_identity_dna(self, input_dna):
        inp = input_dna
        target = True
        res = DNASequence(inp).check_alphabet()
        assert target == res

    def test_identity_rna(self, input_dna):
        inp = input_dna
        with pytest.raises(IncorrectNucleotideError):
            RNASequence(inp).check_alphabet()

    def test_identity_hybrid(self):
        inp = "ACttCcGTUUUUu"
        with pytest.raises(IncorrectNucleotideError):
            DNASequence(inp).check_alphabet()

    def test_transcribe_dna(self, input_dna):
        inp = input_dna
        target = RNASequence("ACuuCcGU").sequence
        res = DNASequence(inp).transcribe().sequence
        assert target == res

    def test_prot_len(self):
        inp = "AlaValSer"
        target = 3
        res = len(AminoAcidSequence(inp, 3))
        assert target == res

    def test_mol_mass(self):
        inp = "ARN"
        target = 89 + 174 + 132
        res = AminoAcidSequence(inp, 1).get_molecular_mass()
        assert target == res


class TestFastqFilter:
    @pytest.fixture
    def tmp_fastq(self):
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
        inp = tmp_fastq
        target = "@Seq1\nACAT\n+\nIIII"
        filter_fastq(inp, None, 50, 5, 10)
        with open("./fastq_filtrator_results/filtered_tmp.fastq", "r") as f:
            res = f.read().rstrip("\n")
        os.remove("./fastq_filtrator_results/filtered_tmp.fastq")
        os.rmdir("./fastq_filtrator_results")
        assert target == res


class TestConvertMultilineOneline:
    @pytest.fixture
    def tmp_fasta(self):
        file_path = "./tmp.fasta"
        with open(file_path, "w") as f:
            f.write(">Seq\nAGg\nCaT")

        yield file_path
        os.remove(file_path)

    def test_fasta_converter(self, tmp_fasta):
        inp = tmp_fasta
        target = ">Seq\nAGgCaT"
        convert_multiline_fasta_to_oneline(inp)
        with open("./oneline_fasta/oneline_tmp.fasta", "r") as f:
            res = f.read().rstrip("\n")
        os.remove("./oneline_fasta/oneline_tmp.fasta")
        os.rmdir("./oneline_fasta")
        assert target == res
