"""This module creates the argument parser and sets up the 
configuration given the command line arguments.
"""

import config
import logger
import util
import pathlib
import software
import os
from argparse import ArgumentParser

log = None
steps = [ # PER RUN
    'runFastQC', 'cropLowQuality', #1-Quality (FastQC, trimming)
    'hostAlignment', 'removeAligned', 'pathogenAlignment', 'sort', #2-Alignment (host, pathogen)
    'runQualimap', 'discardHighError', #3-Qualimap 
    'callReditools', 'callJacusa', 'filterSOR', #4-SNV_Calling
    "runSnpEff", # 5-SnpEff 
    'mergeCalling', 'generateGraphs', # 6-Visualization (files, plots)
    'generateReport'
]
tools = ["hisat2", "bwa", "star", "magicblast", "minimap2", "gmap"]

def _getParser():
    """Creates the argument options for the pipeline.

    Returns:
        ArgumentParser: Command line argument parser.
    """
    parser = ArgumentParser()
    parser.add_argument("-s", "--step", metavar="OPTION", choices=steps, dest="step", type=str, help="only execute a specific step. Available options: 'runFastQC', 'cropLowQuality', 'hostAlignment', 'removeAligned', 'pathogenAlignment', 'sort', 'runQualimap', 'discardHighError', 'callReditools', 'callJacusa', 'filterSOR', 'runSnpEff', 'mergeCalling', 'generateGraphs' and 'generateReport'")
    parser.add_argument("-r", "--resume", metavar="OPTION", choices=steps, dest="resume", type=str, help="resume the process from a specific step. Available options: 'runFastQC', 'cropLowQuality', 'hostAlignment', 'removeAligned', 'pathogenAlignment', 'sort', 'runQualimap', 'discardHighError', 'callReditools', 'callJacusa', 'filterSOR', 'runSnpEff', 'mergeCalling', 'generateGraphs' and 'generateReport'")
    parser.add_argument("-w", "--workingDir", metavar="PATH", dest="workingDir", help="path to working directory")
    parser.add_argument("-c", "--configDir", metavar="CONFIG", dest="configDir", help="path to configuration directory")
    parser.add_argument("-sr", "--sourceType", metavar="OPTION", choices=["project", "sra", "file"], dest="sourceType", help="source type for the list of reads. Available options: 'project' (will read projects.txt), 'sra' (will read sras.txt) and 'file' (will read singleInput.txt, mixedInput.txt or pairedInput.txt depending on the input type)")
    parser.add_argument("-it", "--inputType", metavar="OPTION", choices=["single", "paired", "mixed"], dest="inputType", help="input type for the read files if sourceType = 'file'. Available options: 'single', 'paired' and 'mixed'")
    parser.add_argument("-if", "--inputFastqDir", metavar="PATH", dest="inputFastqDir", help="directory where all the input FASTQ files are located. Only for sourceType = 'file'")
    parser.add_argument("-d", "--download", dest="download", help="download all required tools", action="store_const", const=1)
    parser.add_argument("-hs", "--slurm", dest="slurm", help="run SLURM jobs", action="store_const", const=1)
    parser.add_argument("-ht", "--slurmTime", metavar="HOURS", dest="slurmTime", type=int, help="SLURM time allocation in hours per job")
    parser.add_argument("-hm", "--slurmMem", metavar="MEMORY", dest="slurmMem", type=str, help="SLURM memory allocation per job. For example, '8G' or '8000M'")
    parser.add_argument("-hc", "--slurmCpus", metavar="CPUS", dest="slurmCpus", type=int, help="SLURM CPUs per job")
    parser.add_argument("-hr", "--referenceHostPath", metavar="PATH", dest="referenceHostPath", help="path to reference host genome file")
    parser.add_argument("-prf", "--referencePathogenGenomePaths", metavar="PATH", dest="referencePathogenGenomePaths", help="paths to reference pathogen genome FASTA files (comma-separated)")
    parser.add_argument("-prg", "--referencePathogenGenesPaths", metavar="PATH", dest="referencePathogenGenesPaths", help="paths to reference pathogen genes files (comma-separated). Accepted extensions: gff, gff3, gtf, gbk, gbff, gb, refseq.")
    parser.add_argument("-prp", "--referencePathogenProteinPaths", metavar="PATH", dest="referencePathogenProteinPaths", help="paths to reference pathogen protein FASTA files (comma-separated)")
    parser.add_argument("-ah", "--alignmentSoftwareHost", metavar="OPTION", choices=tools, dest="alignmentSoftwareHost", help="tool for aligning against the reference host genome. Available options: 'bwa', 'hisat2', 'star', 'magicblast', 'minimap2' and 'gmap'")
    parser.add_argument("-ap", "--alignmentSoftwarePathogen", metavar="OPTION", choices=tools, dest="alignmentSoftwarePathogen", help="tool for aligning against the reference pathogen genome. Available options: 'bwa', 'hisat2', 'star', 'magicblast', 'minimap2' and 'gmap'")
    parser.add_argument("-cs", "--cropSoftware", metavar="OPTION", choices=['trimmomatic', 'trimgalore'], dest="cropSoftware", help="tool for cropping the reads. Available options: 'trimmomatic' and 'trimgalore'")
    parser.add_argument("-cas", "--callingSoftware", metavar='OPTION', choices=["jacusa", "reditools", "both"], dest="callingSoftware", help="tool(s) for calling the SNVs, Available options: 'jacusa', 'reditools' and 'both'. If 'both', the common SNVs to both tools will be used for the analysis.")
    parser.add_argument("-skf", "--skipFastQC", dest="skipFastQC", help="Disable the FastQC execution", action="store_const", const=1)
    parser.add_argument("-ske", "--skipDiscardHighError", dest="skipDiscardHighError", help="disable the step that discards the files that got a high alignment error on Qualimap", action="store_const", const=1)
    parser.add_argument("-skc", "--skipCropLowQuality", dest="skipCropLowQuality", help="disable the step that crops the reads based on the results of FastQC", action="store_const", const=1)
    parser.add_argument("-skh", "--skipHostAlignment", dest="skipHostAlignment", help="disable the alignment against the host genome", action="store_const", const=1)
    parser.add_argument("-skr", "--skipRemoveAligned", dest="skipRemoveAligned", help="disable the step that removes the sequences that aligned with the host genome", action="store_const", const=1)
    parser.add_argument("-skq", "--skipQualimap", dest="skipQualimap", help="disable the Qualimap execution", action="store_const", const=1)


    return parser
    
def parseArgs():
    """Parses the arguments and finishes setting up the configuration.

    Returns:
        str: Returns the step that is going to be executed, according to the -s argument. If no -s argument was provided, it returns an empty string.
    """

    args = parser.parse_args()

    if args.workingDir != None:
        config.setWorkPath(args.workingDir)
        
    if args.configDir:
        config.loadConfig(configDir=os.path.join(args.configDir, ''), stepP=args.step, resumeFromP=args.resume)
    else:
        config.loadConfig()

    log = logger.logger
    logger.loadConfig()

    step = ""
    if args.step:
        step = args.step
        if args.step not in steps:
            log.error("Invalid execution step. Accepted values are 'runFastQC', 'cropLowQuality', 'hostAlignment', 'removeAligned', 'pathogenAlignment', 'sort', 'runQualimap', 'removeHighError', 'callReditools', 'callJacusa', 'filterSOR', 'runSnpEff', 'mergeCalling', 'generateGraphs', 'generateReport'.")
            util.stopProgram()

    if args.resume:
        if args.resume not in steps:
            log.error("Invalid execution step. Accepted values are 'runFastQC', 'cropLowQuality', 'hostAlignment', 'removeAligned', 'pathogenAlignment', 'sort', 'runQualimap', 'removeHighError', 'callReditools', 'callJacusa', 'filterSOR', 'runSnpEff', 'mergeCalling', 'generateGraphs', 'generateReport'.")
            util.stopProgram()
    
    if step != "":
        log.info(f"Executing step {step} only...")
        
    if config.resumeFrom != "None":
        log.info(f"Resuming from step {config.resumeFrom}...")
    
    if args.download == 1:
        software.downloadTools()
        log.info("All required software installed correctly.")
        util.stopProgram()

    if args.referenceHostPath != None:
        config.originalHostPath = [args.referenceHostPath]
        
    if args.referencePathogenGenomePaths != None:
        config.originalPathogenGenomePaths = args.referencePathogenGenomePaths.split(",")
    if args.referencePathogenGenesPaths != None:
        config.originalPathogenGenesPaths = args.referencePathogenGenesPaths.split(",")
    if args.referencePathogenProteinPaths != None:
        config.originalPathogenProteinPaths = args.referencePathogenProteinPaths.split(",")
    if (
        len(config.originalPathogenGenesPaths) != len(config.originalPathogenProteinPaths) or
        len(config.originalPathogenGenesPaths) != len(config.originalPathogenGenomePaths) or
        len(config.originalPathogenGenomePaths) != len(config.originalPathogenProteinPaths)
    ):
        log.info("The number of pathogen reference genome, protein and genes files are not the same.")
        util.stopProgram()
    if (
        len(config.originalPathogenGenesPaths) == 0 or
        len(config.originalPathogenProteinPaths) == 0 or
        len(config.originalPathogenGenomePaths) == 0
    ):
        log.info("You must provide at least one genome file, one protein file and one genes file.")
        util.stopProgram()

    for genome, genes, protein in zip(config.originalPathogenGenomePaths, config.originalPathogenGenesPaths, config.originalPathogenProteinPaths):
        path = pathlib.Path(genome)
        refFile = path.name.split(".")
        ref = refFile[0]
        genesExtension = pathlib.Path(genes).name.split(".")[1].lower()
        genesFormat = ""
        if genesExtension == "gtf":
            genesFormat = "gtf22"
        elif genesExtension in ["gff", "gff3"]:
            genesFormat = "gff3"
            genesExtension = "gff"
        elif genesExtension in ["gbk", "gbff", "gb"]:
            genesFormat = "genbank"
            genesExtension = "gbk"
        elif genesExtension == "refseq":
            genesFormat = "refSeq"
        else:
            log.info(f"The genes file {genes} does not have a valid extension. The accepted extensions are gff, gtf, gbk, gbff and refseq.")
            util.stopProgram()
        config.pathogenReferenceGenesFormats.append(genesFormat)
        util.makeDirectory(f"{config.workPath}/data/{ref}/")
        if not os.path.isfile(genome):
            log.info(f"The genome file {refFile} was not found.")
            util.stopProgram()
        else:
            path = f"{config.workPath}/data/{ref}/genome.fa"
            if not os.path.isfile(path):
                util.execCmd(f"cp {genome} {config.workPath}/data/{ref}/genome.fa")
            config.pathogenReferenceGenomePaths.append(f"{config.workPath}/data/{ref}/genome.fa")
        if not os.path.isfile(genes):
            log.info(f"The genes file {refFile} was not found.")
            util.stopProgram()
        else:
            path = f"{config.workPath}/data/{ref}/genes.{genesExtension}"
            if not os.path.isfile(path):
                util.execCmd(f"cp {genes} {config.workPath}/data/{ref}/genes.{genesExtension}")
            config.pathogenReferenceGenesPaths.append(f"{config.workPath}/data/{ref}/genes.{genesExtension}")
        if not os.path.isfile(protein):
            refFile = protein.split("/")[-1]
            log.info(f"The protein file {refFile} was not found.")
            util.stopProgram()
        else:
            path = f"{config.workPath}/data/{ref}/protein.fa"
            if not os.path.isfile(path):
                util.execCmd(f"cp {protein} {config.workPath}/data/{ref}/protein.fa")
            config.pathogenReferenceProteinPaths.append(f"{config.workPath}/data/{ref}/protein.fa")

    if args.sourceType != None:
        config.source = args.sourceType
    if args.inputType != None:
        config.inputType = args.inputType
    if args.inputFastqDir != None:
        config.inputFastqDir = args.inputFastqDir

    if args.slurm == 1:
        config.slurm = True
    if args.slurmTime != None:
        config.slurmTime = args.slurmTime
    if args.slurmMem != None:
        config.slurmMem = args.slurmMem
    if args.slurmCpus != None:
        config.slurmCpus = args.slurmCpus
    if args.alignmentSoftwareHost != None:
        config.alignmentSoftwareHost = args.alignmentSoftwareHost
    if args.alignmentSoftwarePathogen != None:
        config.alignmentSoftwarePathogen = args.alignmentSoftwarePathogen
    if args.cropSoftware != None:
        config.cropSoftware = args.cropSoftware
    if args.callingSoftware != None:
        config.callingSoftware = args.callingSoftware
    util.makeDirectory(f"{config.workPath}/logs/")
    if config.slurm:
        util.makeDirectory(f"{config.workPath}/logs/slurm/")

    if args.skipFastQC == 1:
        config.runFastQC = "Off"
    if args.skipCropLowQuality == 1:
        config.cropLowQualityRuns = "Off"
    if args.skipHostAlignment == 1:
        config.runHostAlignment = "Off"
    if args.skipRemoveAligned == 1:
        config.removeAlignedWithHost = "Off"
    if args.skipQualimap == 1:
        config.runQualimap = "Off"
    if args.skipDiscardHighError == 1:
        config.removeHighErrorRuns = "Off"

    if config.source not in ["project", "file", "sra"]:
        log.error("Invalid source in configuration file. Accepted values are 'project', 'file' and 'sra'.")
        util.stopProgram()

    config.step = step

parser = _getParser()
parseArgs()