"""This module handles everything about quality analysis.
"""

import util
import config
import logger
import glob
import pathlib

log = logger.logger


def runFastQC(sras):
    """Runs a quality analysis on the reads using FastQC.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    qualityDir = config.workPath + "/1-quality"
    fastqcDir = qualityDir + "/fastqc"
    util.makeDirectory(fastqcDir)
    fastqcDirTemp = qualityDir + "/fastqctmp"
    util.makeDirectory(fastqcDirTemp)
    jobs = []
    for sra in sras:
        files = []
        filePath = pathlib.Path(sra[0][0])
        files = f" ".join(sra[0])
        filesLog = " ".join(sra[0])
        log.info(f"Analyzing the quality of {filesLog}...")
        cmd = f"{config.fastqcPath} --extract -t {config.threads} -o {fastqcDir}/ {config.fastqc} {files}"
        util.runCommand(cmd, jobName="fastqc", jobs=jobs)
    util.waitForJobs(jobs)
    cmd = f"rm -r {fastqcDirTemp}"
    util.execCmd(cmd)

def runQualimap(sras):
    """Runs a quality analysis on the alignments using Qualimap.

    It analyzes the quality of the alignment of each run against
    each viral genome.  

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    qualityDir = config.workPath + "/3-qualimap"
    util.makeDirectory(qualityDir)
    bamDir = config.workPath + "/2-alignment/pathogen/bam"
    jobs = []
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name 
        util.makeDirectory(f"{qualityDir}/{ref}")
        for sra in sras:
            run = sra[2]
            if sra[1] == "paired":
                run = run.replace("_1", "")
            util.makeDirectory(f"{qualityDir}/{ref}/{run}")
            files = glob.glob(f"{bamDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.bam")
            if len(files) == 0:
                log.error(f"Sorted BAM file with the alignment for run with ID {run} against pathogen reference {ref} not found.")
                util.stopProgram()
            log.info(f"Analyzing the alignment of run with ID {run} vs {ref}...")
            cmd = f"{config.qualimapPath} bamqc -bam {bamDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.bam -outdir {qualityDir}/{ref}/{run} -nt {config.threads}{config.qualimap}"
            util.runCommand(cmd, jobName="qualimap", jobs=jobs)
    util.waitForJobs(jobs)

def discardHighError(sras):
    """Reads the quality analysis made by Qualimap and discards
    the runs that yielded a low quality alignment.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  

    Returns:
        list: List of filtered tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end)
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID
    """
    log.info(f"Discarding high error runs...")
    qualityDir = config.workPath + "/3-qualimap"
    runs = []
    newSras = []
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name 
        for sra in sras:
            run = sra[2]
            files = glob.glob(f"{qualityDir}/{ref}/{run}/genome_results.txt")
            if len(files) == 0:
                log.error(f"Qualimap results for run with ID {run} against pathogen reference {ref} not found.")
                util.stopProgram()
            with open(f"{qualityDir}/{ref}/{run}/genome_results.txt") as f:
                for line in f:
                    if line.lstrip().startswith("general error rate"):
                        value = float(line.strip().split("=")[1].strip())
                        if value < config.qualimapMaxError:
                            if run not in runs:
                                runs.append(run)
                                newSras.append(sra)
    
    if len(sras) == len(newSras):
        log.info("No runs were discarded.")
    else:
        for sra in sras:
            run = sra[2]
            if run not in runs:
                log.info(f"Discarded run {run}.")
    return newSras
