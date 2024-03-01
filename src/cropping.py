"""This module handles everything about the cropping of the reads.
"""

import util
import config
import logger
import glob
import pathlib

log = logger.logger
qualityDir = config.workPath + "/1-quality"
fastqcDir = qualityDir + "/fastqc"
fastqDir = qualityDir + "/fastq"

def detectFilesToCrop(sras):
    """It reads the FastQC results for each run and decides whether
    they should be cropped or not depending on the quality of the 
    reads.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  

    Returns:
        list: List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end)
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID
    """
    toCrop = []
    for sra in sras:
        for f in sra[0]:
            added = False
            runDir = pathlib.Path(f).name.replace(".fastq", "_fastqc")
            files = glob.glob(f"{fastqcDir}/{runDir}/fastqc_data.txt")
            if len(files) == 0:
                print(f"{fastqcDir}/{runDir}/fastqc_data.txt")
                log.error(f"FastQC analysis for file {f} not found.")
                util.stopProgram()
            with open(f"{fastqcDir}/{runDir}/fastqc_data.txt") as f:
                values = []
                for _ in range(13):
                    f.readline()
                for line in f:
                    if line.startswith(">>END_MODULE"):
                        break
                    line = line.strip()
                    val = float(line.split()[1])
                    values.append(val)
                maxValue = -1
                for val in values:
                    if val < config.cropMinMeanQuality:
                        toCrop.append(sra)
                        added = True
                        break
                    elif val > maxValue:
                        maxValue = val
                    elif maxValue - config.cropMaxDecay > val:
                        toCrop.append(sra)
                        added = True
                        break
                if added:
                    break
    return toCrop

def runTrimmomatic(sras):
    """Crops the given runs using Trimmomatic.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    util.makeDirectory(fastqDir)
    jobs = []
    for f in sras:
        filePaths = f[0]
        fileNames = []
        for fp in filePaths:
            fileNames.append(pathlib.Path(fp).name)
        runType = f[1]
        runId = f[2]
        filesString = " and ".join(fileNames)
        log.info(f"Cropping {filesString}...")
        if runType == "single":
            cmd = f"{config.javaPath} -jar {config.trimmomaticPath} SE {filePaths[0]} {fastqDir}/{runId}.fastq CROP:{config.cropSize}"
        else:
            cmd = f"{config.javaPath} -jar {config.trimmomaticPath} PE {filePaths[0]} {filePaths[1]} {fastqDir}/{runId}_1.fastq {fastqDir}/{runId}_1.fastq.unpaired {fastqDir}/{runId}_2.fastq {fastqDir}/{runId}_2.fastq.unpaired CROP:{config.cropSize}"
        util.runCommand(cmd, jobName="trimmomatic", jobs=jobs)
    util.waitForJobs(jobs)

def runTrimGalore(sras):
    """Crops the given runs using Trim Galore.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    util.makeDirectory(fastqDir)
    jobs = []
    for f in sras:
        filePaths = f[0]
        fileNames = []
        for fp in filePaths:
            fileNames.append(pathlib.Path(fp).name)
        fileStems = []
        for fp in filePaths:
            fileStems.append(pathlib.Path(fp).stem)
        runType = f[1]
        runId = f[2]
        filesString = " and ".join(fileNames)
        log.info(f"Cropping {filesString}...")
        if runType == "single":
            cmd = f"{config.trimGalorePath} -o {fastqDir} --hardtrim5 {config.cropSize} {filePaths[0]}"
            cmd += f"; mv {fastqDir}/{fileStems[0]}.{config.cropSize}bp_5prime.fq {fastqDir}/{runId}.fastq"
            util.runCommand(cmd, jobName="trimgalore", jobs=jobs)
        else:
            cmd = f"{config.trimGalorePath} -o {fastqDir} --paired --hardtrim5 {config.cropSize} {filePaths[0]} {filePaths[1]}"
            cmd = f"; mv {fastqDir}/{fileStems[0]}.{config.cropSize}bp_5prime.fq {fastqDir}/{runId}_1.fastq"
            cmd = f"; mv {fastqDir}/{fileStems[1]}.{config.cropSize}bp_5prime.fq {fastqDir}/{runId}_2.fastq"
            util.runCommand(cmd, jobName="trimgalore", jobs=jobs)
    util.waitForJobs(jobs)
