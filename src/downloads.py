import util
import xml.etree.ElementTree as ET
import logger
import config
import glob
import re
import os
from Bio import Entrez

Entrez.email = "dummy@dummy.com"

# log = util.getLogger()
log = logger.logger

def getExperimentsList():
    with open(f"../config/projects.config", "r") as f:
        # Skip the first line
        line = f.readline()
        line = f.readline().strip()
        types = []
        ids = []
        while line != None and line != "":
            if line.startswith("#"):
                continue
            project = line.split()[0]
            type = line.split()[1]
            # Search each project and get the list of experiment IDs
            search = Entrez.esearch(db="sra", term=f"{project}[BioProject]", rettype="gb", retmode="text")
            read = Entrez.read(search)
            search.close()
            ids += read["IdList"]
            for _ in read["IdList"]:
                types.append(type)
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
            for el in root.iter("RUN"):
                sras.append((el.attrib["accession"], types[typeCounter]))
                typeCounter += 1
        return sras

def downloadSRAs(sras):
    fastqDir = config.workPath + "/fastq"
    util.makeDirectory(fastqDir)
    auxSras = []
    for sra in sras:
        files = None
        if os.path.exists(f"{fastqDir}/{sra[0]}_1.fastq") or os.path.exists(f"{fastqDir}/{sra[0]}_2.fastq") or os.path.exists(f"{fastqDir}/{sra[0]}.fastq"):
            log.info(f"{sra[0]} already downloaded. Skipping...")
            files = glob.glob(f"{fastqDir}/{sra[0]}*.fastq")
            files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
            for i in range(len(files)):
                files[i] = files[i].split("/")[-1]
        else:
            log.info(f"Downloading {sra[0]}...")
            cmd = f"{config.sratoolkitPath}fastq-dump --outdir {fastqDir} --split-files {sra[0]}"
            util.execCmd(cmd)
            if sra[1] == "single":
                files = glob.glob(f"{fastqDir}/{sra[0]}*.fastq")
                files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
                for file in files:
                    cmd = f"cat {file}"
                    util.execCmd(cmd, f"{fastqDir}/{sra[0]}.fastq", "a")
            else:
                files = glob.glob(f"{fastqDir}/{sra[0]}*.fastq")
                files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
                for i in range(len(files)):
                    files[i] = files[i].split("/")[-1]
        if sra[1] == "single":
            auxSras.append(([f"{sra[0]}.fastq"], "single"))
        else:
            auxSras.append((files, "paired"))
    return auxSras
                