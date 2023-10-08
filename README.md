# SeqMaster

This repository contains an open-source library which make work with protein, DNA and RNA. For proteins it is capable to process multiple proteins and peptides sequences, calculate physical features, find specific sites. For nucleic acids it can transcribe, find complement sequences, reverse, calculate GC-content. Also SeqMaster supports FASTAQ-format and can be used to filter DNA sequences based on their length, GC-content and sequence quality. 

## Installation

To use this toolbox one need to clone repository

```shell
git clone https://github.com/greenbergM/SeqMaster
cd SeqMaster
```

### System requirements:

Key packages and programs:
- [Python](https://www.python.org/downloads/) (version >= 3.9)

## Usage

```python
# import main function
from SeqMaster import run_dna_rna_tools, run_protein analysis, run_filter_fastaq

```

## Sequence tools

This section contains description of functions you can find in our library.

### DNA and RNA tools

Main function: `run_dna_rna_tools(*args, seq_type = 'DNA')`

#### Arguments

**args**:
- DNA or RNA sequnces in `str` format
- operation name in str format

    **Note**
	- There can be any number of sequnces
	- Sequences can contain upper case and lower case letters
	- Operation always must be the last argument!

**seq_type**:
- type of desired complementary sequence
- by default has value 'DNA'
- only supported values are 'DNA' and 'RNA'

#### Operations

- `nucl_acid_identity` - shows the type of given sequences; if sequence cannot be identified as `'DNA'` or `'RNA'` it will be given name `'uncertain'`
	*Returns* 
	`result: dict` 

- `reverse` - gives the reversed copy of a sequence
	*Returns*
	`result: Union[list[str], str]` 

- `complement` - gives complement `DNA` or `RNA` sequence depending on the `seq_type` argument. 
	*Returns*
	`result: Union[list[str], str]` 

- `reverse_complement` - gives reversed complement `DNA` or `RNA` sequence depending on the `seq_type` argument. 
	*Returns*
	`result: Union[list[str], str]` 

- `transcribe` - transcribes given sequences; DNA sequences are transcribed to RNA, RNA sequences are reversed transcribed to DNA.
	*Returns*
	`result: Union[list[str], str]` 


- `gc_content` - gives the GC-content in percent for given sequences
	*Returns*
	`result: Union[list[float], float]`

**Example**

```python
seq1= 'AGTC'
seq2 = 'AUGGG'

run_dna_rna_tools(seq1, 'nucl_acid_identity')
run_dna_rna_tools(seq1, 'complement', seq_type = 'RNA')
run_dna_rna_tools(seq1, seq2, 'transcribe')

```

### Protein tools

Main function: `run_protein_analysis(*args, site_of_interest = None)`

#### Arguments

**args**:
- protein sequnces in `str` format
- sequence encoding type in `int` format (`1` - one-letter encoding or `3` - three-letter encoding)
- operation on `str` format

    **Note**
	- There can be any number of sequnces
	- Sequences can contain upper case and lower case letters
	- Sequences can be given in one-letter format or three-letter format (it must be clarified by sequence encoding type parameter); in three-letter format they can contain spaces between amino acids. 
	- Arguments always must be in the strict order: sequences, sequence encoding type, operation!

**site_of_interest** 
 - sequence of site to find in given protein sequences via `find_site` operation
 - must be in one-letter encoding
 - by default set to be `None`

 #### Operations

- `get_seq_characteristic` - counts entry of each residue type in given sequences.
 	*Returns*
	`result: Union[list[dict], dict]` 

- `find_site` - find if sequence contains certain site (given by `site_of_interest` argument) and get positions of this site
	*Returns*
	`result: Union[list[str], str]` 

- `calculate_protein_mass` - calculates mass of residues in given sequences in Da
	*Returns*
	`result: Union[list[float], float]` 

- `calculate_average_hydrophobicity` - calculates hydrophobicity index for each given sequence as sum of index for each residue in sequence divided by its length
	*Returns*
	`result: Union[list[float], float]` 

- `calculate_isoelectric_point` - calculates soelectrinc point for each given sequence as sum of known pI for residues in sequence
	*Returns*
	`result: Union[list[float], float]` 

- `get_mrna` - gives encoding mRNA nucleotides in degenerate form for each given sequence
	*Returns*
	`result: Union[list[str], str]` 

**Example**

```python
seq1= 'His Val Leu Ile'
seq2 = 'HVGGLLA'

run_protein_analysis(seq1, 3,  'get_mrna')
run_protein_analysis(seq2, 1, 'find_site', site_of_interest = 'GGL')
run_protein_analysis(seq1, seq1, 3, 'calculate_average_hydrophobicity')

```


### FASTAQ filter

Main function: `run_filter_fastaq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0)`


#### Arguments

**seqs**:
- FASTAQ sequences organised in `dict`: key = sequence name, value = sequence (`str`) and sequence quality in ASCII code (`str`)

**gc_bounds**
- lower and upper boundaries of GC-content in percent orginised in `tuple`
- by default is set to `(0, 100)`
- if one number given it will be considered the upper boundary with the lover boundary = 0

**length_bounds**
- lower and upper boundaries of sequence length orginised in `tuple`
- by default is set to `(0, 2 ** 32)`
- if one number given it will be considered the upper boundary with the lover boundary = 0

**quality_threshold**
- threshold average Q-score (`int`) - lower boundary 
- by default is set to 0

#### Operations

Function filters sequences based on given arguments. Sequences will stay in output only if they fit all the given criteria. 
	*Returns*
	`result: dict` 


**Example**
```python
seqs = {'name1': ('ACGTGTAAAGC', 'jjG#HHFddd@'), 'name2': ('GGGTCA', '!@jHHj')}

run_filter_fastaq(seqs, gc_bounds=(20, 70), length_bounds=(4, 100), quality_threshold=27)

```











