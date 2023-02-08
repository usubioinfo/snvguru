# VirnaSNV

VirnaSNV is a Python pipeline for finding and analyzing single nucleotide variants (SNVs) from viral RNA-seq data. It can be used on a mix of single-end and paired-end RNA-seq runs given a list of projects, SRAs or file names, and the alignment can be executed against multiple viral reference genomes.

## Before installing

### Setting up Reditools2

Since Reditools2 uses Python 2, and VirnaSNV uses Python 3, and you do not have Reditools2 available as an HPC module, you have the following options:

* Install Reditools2 on Python 2 default environment, and on the config/main.config file set up reditoolsCommand as <pre>python2 /path/to/reditools.py</pre>.

* Install VirnaSNV on Python 3 default environment, and if Reditools2 runs on a virtual environment, activate Reditools virtual environment and run VirnaSNV...

* On the config/main.config file, set up reditoolsCommand as <pre>~/path/to/environment/for/reditools2/bin/python /path/to/reditools.py</pre>

## Installation

You can install the package through PyPI running the following command:

<pre>
pip3 install virnasnv
</pre>

## Configuration Files

### main.config

The main parameters you must be aware of are:

* `source`, where you can specify whether the list of runs will be extracted from a list of projects (`project`), a list of SRAs (`sra`) or file names (`file`).
* `workPath`, where you specify the working directory for your temporary files and result files.
* `hostReferencePath`, where you say what is the path to the host reference genome.
* `virusReferencePath`, a directory where all the viral genomes to be evaluated are found.
* `cropSoftware`, where you choose the tool for cropping the reads, either `trimmomatic` or `trimgalore`.
* `alignmentSoftware`, where you specify the tool for aligning the reads against the human and the viral genomes. You can choose `hisat2`, `star`, `bwa`, `magicblast`, `minimap2` or `gmap`.
* `reditoolsCommand`, where you put the command to run Reditools2 as instructed in the section **Setting up Reditools2**.
* All the options that end with `Path`, such as `sratoolkitPath` or `fastqcPath` in order to configure the path to the executable or directory where the executables are for each tool.

### projects.config

If you chose `source project`, the pipeline will read this configuration file that lists the sequencing projects, each one stating whether they are single-end or paired-end.

### sra.config

If you chose `source sra`, the pipeline will read this configuration file that lists the SRAs, each one stating whether they are single-end or paired-end.

### files.config

If you chose `source file`, the pipeline will read this configuration file that lists the file names/paths without extension (.fastq), each one stating whether they are single-end or paired-end, and it will work with all files with the .fastq extension. For example, if the file path is ../test/run_1.fastq and ../test/run_2.fastq, and file.config has ../test/run, both files will be considered.

### align\[Tool].config

These are configuration files exclusive to each alignment tool. The most important parameter to modify are the ones that end in `Path` in order to configure the path to the executable or directory where the executables are for each tool. Most of the parameters have two columns: one for the human reference file and one for the viral reference files.

## Results

The output files are found in the `results` directory of the working directory. Each viral reference genome will contain:

* Common SNVs count graph
* Recurrent SNVs count graph
* SNV frequency graph
* SNV distribution graph
* Common SNVs Excel file
* Recurrent SNVs Excel file
* Each run will also contain:
    * Common SNVs count graph
    * Recurrent SNVs count graph
    * Common SNVs Excel file
    * Recurrent SNVs Excel file
    * JACUSA Excel file
    * Reditools2 Excel file
