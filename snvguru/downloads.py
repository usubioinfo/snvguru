"""This module handles everything about the download of the runs.
"""

from snvguru import util
import xml.etree.ElementTree as ET
from snvguru import logger
from snvguru import config
import glob
import re
import os
from Bio import Entrez

Entrez.email = "dummy@dummy.com"

log = logger.logger

def getExperimentsList():
    """Gets the accession numbers (SRAs) of all the runs within the
    projects found in the projects.txt file from the NCBI database.

    Returns:
        sras (list): List of tuples with the following data:
            0 - A list of the paths for the input run (one file if single-end, two if paired-end)
            1 - Run type. "single" if single-end, "paired" if paired-end
            2 - Run ID
    """
    with open("projects.txt", "r") as f:
        # Skip the first line
        line = f.readline()
        line = f.readline().strip()
        types = []
        ids = []
        while line != None and line != "":
            if line.startswith("#"):
                line = f.readline().strip()
                continue
            project = line.split()[0]
            # Search each project and get the list of experiment IDs
            search = Entrez.esearch(db="sra", term=f"{project}[BioProject]", rettype="gb", retmode="text", retmax="10000")
            read = Entrez.read(search)
            search.close()
            ids += read["IdList"]
            count = int(read["Count"])
            if count == 0:
                log.info(f"No experiments found for project {project}.")
            else:
                log.info(f"{count} experiments found for project {project}.")
            line = f.readline().strip()

        # Get sras for the experiments
        search = Entrez.epost(db='sra', id=",".join(ids), retmode='text', rettype='gb')
        read = Entrez.read(search)
        search.close()
        webenv = read["WebEnv"]
        key = read["QueryKey"]
        batchSize = 200
        sras = []
        for start in range(0, count, batchSize):
            search = Entrez.efetch(
                db="sra",
                retmode="xml",
                retstart=start,
                retmax=batchSize,
                webenv=webenv,
                query_key=key,
                idtype="acc",
            )
            data = search.read()
            search.close()
            root = ET.fromstring(data)
            typeCounter = 0
            for el in root.iter("RUN_SET"):
                for el2 in el.iter("RUN"):
                    sras.append(el2.attrib["accession"])
                typeCounter += 1
        return sras

import subprocess
import shutil
import urllib.request
import tarfile

def _ensureSraToolkit():
    """Checks if prefetch works on the host system.
    If it fails due to GLIBC incompatibility (e.g. on CentOS 7), it automatically
    downloads and configures NCBI's static CentOS build.
    """
    cmd = f"{config.sratoolkitPath}prefetch --version"
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    err = res.stderr.decode("utf-8")
    out = res.stdout.decode("utf-8")
    
    if res.returncode != 0 and ("GLIBC" in err or "GLIBC" in out or "not found" in err or "not found" in out):
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        tools_dir = os.path.join(repo_root, "tools")
        sra_bin_dir = os.path.join(tools_dir, "sratoolkit", "bin")
        if not os.path.exists(os.path.join(sra_bin_dir, "prefetch")):
            log.info("Host GLIBC compatibility issue detected with sra-tools. Automatically downloading NCBI static binaries for CentOS/Linux...")
            os.makedirs(tools_dir, exist_ok=True)
            tar_path = os.path.join(tools_dir, "sratoolkit.tar.gz")
            url = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
            try:
                urllib.request.urlretrieve(url, tar_path)
                with tarfile.open(tar_path, "r:gz") as tar:
                    tar.extractall(path=tools_dir)
                extracted_dirs = glob.glob(os.path.join(tools_dir, "sratoolkit*centos*"))
                if extracted_dirs:
                    target_dir = os.path.join(tools_dir, "sratoolkit")
                    if os.path.exists(target_dir):
                        shutil.rmtree(target_dir)
                    os.rename(extracted_dirs[0], target_dir)
                if os.path.exists(tar_path):
                    os.remove(tar_path)
            except Exception as e:
                log.error(f"Failed to auto-download static SRA toolkit: {e}")
                return
        config.sratoolkitPath = sra_bin_dir + "/"
        log.info(f"Using NCBI static SRA toolkit from {config.sratoolkitPath}")

def downloadSRAs(sras):
    """Downloads the FASTQ files for each run using the prefetch and
    fasterq-dump tools from SRAtoolkit.

    Args:
        sras (list): List of accession numbers

    Returns:
        list: List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end)
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID
    """
    _ensureSraToolkit()
    fastqDir = config.workPath + "/data/fastq"
    util.makeDirectory(fastqDir)
    sraDir = config.workPath + "/data/sra"
    util.makeDirectory(sraDir)
    auxSras = []
    runType = ""
    for sra in sras:
        files = None
        if os.path.exists(f"{fastqDir}/{sra}_1.fastq") or os.path.exists(f"{fastqDir}/{sra}_2.fastq") or os.path.exists(f"{fastqDir}/{sra}.fastq"):
            log.info(f"{sra} already downloaded. Skipping...")
            files = glob.glob(f"{fastqDir}/{sra}*.fastq")
            files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
            for i in range(len(files)):
                files[i] = files[i].split("/")[-1]
            if len(files) > 1: 
                runType = "paired"
            elif len(files) == 1:
                runType = "single"
        else:
            log.info(f"Downloading {sra}...")
            cmd = f"{config.sratoolkitPath}prefetch {sra} --output-directory {sraDir}"
            out, err = util.execCmd(cmd)
            cmd = f"{config.sratoolkitPath}fasterq-dump {sraDir}/{sra} --outdir {fastqDir} --split-3"
            util.execCmd(cmd)
            files = glob.glob(f"{fastqDir}/{sra}*.fastq")
            if len(files) == 1:
                for f in files:
                    cmd = f"cat {f}"
                    util.execCmd(cmd, f"{fastqDir}/{sra}.fastq", "a")
                runType = "single"
            else:
                files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
                runType = "paired"
            
        if runType == "single":
            auxSras.append(([f"{fastqDir}/{sra}.fastq"], "single", sra))
        else:
            auxSras.append(([f"{fastqDir}/{sra}_1.fastq", f"{fastqDir}/{sra}_2.fastq"], "paired", sra))
    cmd = f"rm -r {sraDir}"
    util.execCmd(cmd)

    return auxSras
                