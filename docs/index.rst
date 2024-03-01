.. SNVGuru documentation master file, created by
   sphinx-quickstart on Mon Apr 17 02:40:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SNVGuru's documentation!
=======================================

============
Introduction
============

SNVGuru is an RNA-seq analysis tool made in Python that downloads and filters high-quality reads, discarding reads that align with a host genome, and it calls and analyzes the single nucleotide variants found. It supports multiple aligning tools, and uses JACUSA, REDItools2 and SnpEff for calling the SNVs. At the end, you will get an HTML report with the basic parameters used and the explanation of every generated figure.

===============
How to install?
===============

* Download SNVGuru from GitHub running `git clone <URL>`.
* Run `cd snvguru`.
* Run `conda env create -f pipeline_environment.yml`. Be aware that you must have Miniconda or Anaconda installed (see :ref:`miniconda-installation`).
* Run `pip install -r requirements.txt`.
* SNVGuru can install the required tools for you. For that, run `python3 src/main.py -d`. This is only needed the first time. These tools will be located in the `tools` folder.

===========
How to run?
===========

For running SNVGuru, the command is `python3 src/main.py`. It can read the configuration (including the input files) from the `config/main.config` file, or you can use the multiple arguments to customize your execution. Use `python3 src/main.py -h` for a description of all available arguments.

=================
How to configure?
=================
  
In the configuration folder (`config/` or your folder of choice using the `-c` argument) you will find 12 different `.config` files. 11 of these are for the tool the name refers to, and, in general, you will not have to modify these, unless you are using Minimap2 or DNA sequences with Magic-BLAST. For example, `bwa.config` is for BWA. The most important configuration file, and the one you might want to check and modify to suit your needs, is `main.config`. The most important parameters you might want to configure here are `source`, `inputType` (if `source` is `file`), `inputFastqDir` (if `source` is `file`), `workPath`, `hostReferencePath`, `pathogenReferenceGenomePaths`, `pathogenReferenceProteinPaths`, `pathogenReferenceGenesPaths`, `alignmentSoftwareHost` (if you want to eliminate the host-contaminated reads first) and `alignmentSoftwarePathogen`.

* `source` (`-s`): It can be 'project', 'file' or 'sra'.
   * `project`: It will read a list of BioProject IDs from `projects.txt`.
   * `sra`: It will read a list of SRA IDs from sras.txt.
   * `file`: It will read a a list of files from `singleInput.txt`, `pairedInput.txt` or `mixedInput.txt`, depending on the `inputType` value.
      * `inputType` (`-it`): It can be either 'single', 'paired' or 'mixed'. Will only work if `source` is `file`. All files read must be located at `inputFastqDir`.
         * `single`: It will read the files in `singleInput.txt`. It has three columns: Run (the sample long ID), ID (the sample short ID for pipeline use) and File (the file name). All reads must be single end.
         * `paired`: It will read the files in `pairedInput.txt`. It has four columns: Run (the sample long ID), ID (the sample short ID for pipeline use), and File1 and File2 (the file names of the main and the mate reads). All reads must be paired end.
         * `mixed`: It will read the files in `mixedInput.txt`. It has five columns: Run (the sample long ID), ID (the sample short ID for pipeline use), Type (either `single` or `paired`) and File1 and File2 (the file names of the main and the mate reads). File2 is not required if the sample is single-end.
      * `inputFastqDir` (`-if`): Directory where all input FASTQ files from the samples are located. Will only work if `source` is `file`.
* `workPath` (`-w`): When you download SNVGuru, it will have the value `workspace`, which means that your results will be located at `workspace/`. If you want to run the pipeline with different configurations, you might want to have a different `workPath` for every configuration.
* `hostReferencePath` (`-hr`): Location of the host reference genome FASTA file.
* Pathogen reference files: Each pathogen reference has three files needed: The genome FASTA file, the proteome FASTA file and the genes file. If you are running the samples against multiple genomes, make sure that they are input in the same order for the three following parameters. 
   * `pathogenReferenceGenomePaths` (`-prf`): Location of the pathogen reference genome FASTA files. If you are running the samples against multiple genomes, they must be separated by comma.
   * `pathogenReferenceProteinPaths` (`-prp`): Location of the pathogen reference proteome FASTA files. If you are running the samples against multiple genomes, they must be separated by comma.
   * `pathogenReferenceGenesPaths` (`-prg`): Location of the pathogen reference genes file. Accepted formats are GFF (`.gff`, `.gff3`), GTF (`.gtf`), GenBank (`.gbk`, `.gbff`, `.gb`) or RefSeq (`.refseq`).
* Alignment tools: There are two parameters for setting the tools used for the alignment steps:
   * `alignmentSoftwareHost` (`-ah`): Selected tool for running the alignment against the host.
   * `alignmentSoftwarePathogen` (`-ap`): Selected tool for running the alignment against the pathogens.
   * These tools can be:
      * `hisat2`: Hisat2 is suggested for short RNA-seq reads. 
      * `star`: STAR is suggested for short RNA-seq reads. 
      * `bwa`: BWA is suggested for short DNA reads.
      * `minimap2`: Minimap2 is suggested for long DNA or RNA-seq reads. 
      * `gmap`: GMAP is suggested for long cDNA reads.
      * `magicblast`: Magic-BLAST can be used for any type of read. 


==========================================================
Do you have a sample report? How to interpret the figures?
==========================================================

You can check `this sample report <influenzaA/analysis_report.html>`_ for *influenza A*, or `this one <mtuberculosis/analysis_report.html>`_ for *Mycobacterium tuberculosis*, or `this other one <hcapsulatum/analysis_report.html>`_ for *Histoplasma capsulatum*. 

.. _miniconda-installation:

=========================
How to install Miniconda?
=========================

* Download the installer from https://docs.conda.io/en/latest/miniconda.html#linux-installers.
* Run `bash Miniconda3-latest-Linux-x86_64.sh`. The filename can change.
* Accept all the default configuration.
* Close and reopen the terminal window.
* You can test that it is installed by running `conda list`. It should display a list of installed packages.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
