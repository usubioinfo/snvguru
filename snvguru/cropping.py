"""This module handles everything about the cropping of the reads.
"""

from snvguru import util
from snvguru import config
from snvguru import logger
import glob
import pathlib
import os
import shutil

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
                in_module = False
                for line in f:
                    line = line.strip()
                    if line.startswith(">>Per base sequence quality"):
                        in_module = True
                        continue
                    if in_module:
                        if line.startswith(">>END_MODULE"):
                            break
                        if line.startswith("#"):
                            continue
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                val = float(parts[1])
                                values.append(val)
                            except ValueError:
                                pass
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
            cmd = f"{config.trimmomaticPath} SE {filePaths[0]} {fastqDir}/{runId}.fastq CROP:{config.cropSize}"
        else:
            cmd = f"{config.trimmomaticPath} PE {filePaths[0]} {filePaths[1]} {fastqDir}/{runId}_1.fastq {fastqDir}/{runId}_1.fastq.unpaired {fastqDir}/{runId}_2.fastq {fastqDir}/{runId}_2.fastq.unpaired CROP:{config.cropSize}"
        util.runCommand(cmd, jobName="trimmomatic", jobs=jobs)
    util.waitForJobs(jobs)

import subprocess
import urllib.request

def _ensureTrimGalore():
    """Checks if trim_galore works on the host system.
    If it fails due to GLIBC incompatibility (e.g. on CentOS 7), it automatically
    downloads and configures the standalone TrimGalore perl script from GitHub.
    """
    cmd = f"{config.trimGalorePath} --version"
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    err = res.stderr.decode("utf-8")
    out = res.stdout.decode("utf-8")
    
    if res.returncode != 0 and ("GLIBC" in err or "GLIBC" in out or "not found" in err or "not found" in out):
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        tools_dir = os.path.join(repo_root, "tools")
        standalone_tg = os.path.join(tools_dir, "trim_galore")
        if not os.path.exists(standalone_tg):
            log.info("Host GLIBC compatibility issue detected with trim_galore. Automatically downloading standalone Trim Galore script...")
            os.makedirs(tools_dir, exist_ok=True)
            url = "https://raw.githubusercontent.com/FelixKrueger/TrimGalore/0.6.10/trim_galore"
            try:
                urllib.request.urlretrieve(url, standalone_tg)
                os.chmod(standalone_tg, 0o755)
            except Exception as e:
                log.error(f"Failed to auto-download standalone trim_galore: {e}")
                return
        config.trimGalorePath = standalone_tg
        log.info(f"Using standalone Trim Galore from {config.trimGalorePath}")

def runTrimGalore(sras):
    """Crops the given runs using Trim Galore.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    _ensureTrimGalore()
    util.makeDirectory(fastqDir)
    jobs = []
    for f in sras:
        filePaths = f[0]
        fileNames = [pathlib.Path(fp).name for fp in filePaths]
        runType = f[1]
        filesString = " and ".join(fileNames)
        log.info(f"Cropping {filesString}...")
        if runType == "single":
            cmd = f"{config.trimGalorePath} -o {fastqDir} --hardtrim5 {config.cropSize} {filePaths[0]}"
            util.runCommand(cmd, jobName="trimgalore", jobs=jobs)
        else:
            cmd = f"{config.trimGalorePath} -o {fastqDir} --paired --hardtrim5 {config.cropSize} {filePaths[0]} {filePaths[1]}"
            util.runCommand(cmd, jobName="trimgalore", jobs=jobs)
    util.waitForJobs(jobs)

    # Rename Trim Galore output files to standard {runId}.fastq / {runId}_1.fastq naming convention
    for f in sras:
        filePaths = f[0]
        fileStems = [pathlib.Path(fp).stem for fp in filePaths]
        runType = f[1]
        runId = f[2]
        if runType == "single":
            candidates = [
                f"{fastqDir}/{fileStems[0]}.{config.cropSize}bp_5prime_trimmed.fq",
                f"{fastqDir}/{fileStems[0]}.{config.cropSize}bp_5prime.fq",
                f"{fastqDir}/{fileStems[0]}_trimmed.fq",
                f"{fastqDir}/{fileStems[0]}.fq",
            ]
            dest = f"{fastqDir}/{runId}.fastq"
            for cand in candidates:
                if os.path.exists(cand):
                    shutil.move(cand, dest)
                    break
        else:
            candidates_1 = [
                f"{fastqDir}/{fileStems[0]}.{config.cropSize}bp_5prime_val_1.fq",
                f"{fastqDir}/{fileStems[0]}.{config.cropSize}bp_5prime.fq",
                f"{fastqDir}/{fileStems[0]}_val_1.fq",
            ]
            dest_1 = f"{fastqDir}/{runId}_1.fastq"
            for cand in candidates_1:
                if os.path.exists(cand):
                    shutil.move(cand, dest_1)
                    break

            candidates_2 = [
                f"{fastqDir}/{fileStems[1]}.{config.cropSize}bp_5prime_val_2.fq",
                f"{fastqDir}/{fileStems[1]}.{config.cropSize}bp_5prime.fq",
                f"{fastqDir}/{fileStems[1]}_val_2.fq",
            ]
            dest_2 = f"{fastqDir}/{runId}_2.fastq"
            for cand in candidates_2:
                if os.path.exists(cand):
                    shutil.move(cand, dest_2)
                    break