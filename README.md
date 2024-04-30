# Owerview

This repository contains several Python scripts developed for Python coursework during my studies at the Bioinformatics Institute. There are three main Python scripts (see below). Example usage for each can be found in the `Showcases.ipynb` notebook. Additionally, the repository includes a testing script containing several assert tests.


## SeqMaster.py 

This script represents my work on object-oriented programming (OOP) assignments (hw14). It includes classes designed for handling different biological sequences such as DNA, RNA, and Amino Acid sequences. These classes offer functionalities for basic sequence analysis, such as finding complementary sequences for nucleic acids and calculating molecular masses for amino acids. Additionally, there's a function called `filter_fastq`, leveraging Biopython, enabling FASTQ file filtering based on user-defined parameters.

Additionally, SeqMaster.py contains functions related to API homework (hw17). It comprises two components: an API for a Telegram bot capable of calculating function working time and sending log files, and an API for utilizing the GENSCAN web-page.


## bio_files_processor.py

This script contains my work on iterators (hw15). There is code for iterator, that allow for reading by records in FASTA file. 
Additionally, there are some functions from earlier assignment (hw6) focusing on working with various biological files, like parsing GenBank and BLAST outputs, converting multi-line FASTA to one-line format and altering the starting position in FASTA files. These functions rely on several helper functions located in the `additional_scripts` folder within `bio_files_processor_scripts.py` (That is the demonstartion that I know how to make modules :) ). 


## custom_random_forest.py

Script that contains my realisation of RandomForestClassifier with integrated multithreading.


## Showcase

`Showcases.ipynb` notebook contains some examples for `SeqMaster.py`, `bio_files_processor.py` and `custom_random_forest.py`. `Data` folder contains some files for demonstration.

## Testing 

`test_SeqMaster.py` script contains assert tests for `SeqMaster.py` and `bio_files_processor.py`



# Dependencies

All the dependencies can be found in `evironment.yml` file 

