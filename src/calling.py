import util
import config
import logger
import io
import pandas as pd
import numpy as np
import glob 

log = logger.logger

def runReditools2(sras):
    indexDir = config.workPath + "/indices"
    homopolDir = config.workPath + "/homopol"
    util.makeDirectory(homopolDir)
    callingDir = config.workPath + "/calling"
    util.makeDirectory(callingDir)
    bamDir = config.workPath + "/bam"
    for fullRef in config.virusReferencePath:
        cmd = f"cp {fullRef} {indexDir}"
        util.execCmd(cmd)
        refFile = fullRef[fullRef.rfind("/") + 1:]
        ref = refFile[:refFile.rfind(".")]        
        util.makeDirectory(f"{homopolDir}/{ref}")
        util.makeDirectory(f"{callingDir}/{ref}")
        cmd = f"{config.samtoolsPath} faidx {indexDir}/{refFile}"
        util.execCmd(cmd)
        for sra in sras:
            run = sra[0][0]
            run = run[:run.rfind(".")]
            type = sra[1]
            if type == "paired":
                run = run.replace("_1", "")
            file = f"{run}_virus_{config.alignmentSoftwareVirus}"
            files = glob.glob(f"{bamDir}/{ref}/{file}.bam")
            if len(files) == 0:
                log.error(f"Sorted BAM file with the alignment for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            log.info(f"Running SNV calling with Reditools2: {run} vs {ref}...")
            cmd = f"{config.samtoolsPath} index {bamDir}/{ref}/{file}.bam"
            util.execCmd(cmd)
            cmd = f"python {config.reditoolsCommand} -c -f {bamDir}/{ref}/{file}.bam -r {indexDir}/{refFile} -m {homopolDir}/{ref}/{run}.homopol.txt -o {callingDir}/{ref}/{run}.calling.txt -S -s {config.reditoolsStrand} -os {config.reditoolsHomopolSpan}{config.reditools}"
            util.execCmd(cmd)

def runJacusa(sras):
    indexDir = config.workPath + "/indices"
    callingDir = config.workPath + "/calling"
    bamDir = config.workPath + "/bam"
    for fullRef in config.virusReferencePath:
        cmd = f"cp {fullRef} {indexDir}"
        util.execCmd(cmd)
        ref = fullRef[fullRef.rfind("/") + 1:]
        ref = ref[:ref.rfind(".")]
        for sra in sras:
            run = sra[0][0]
            run = run[:run.rfind(".")]
            type = sra[1]
            if type == "paired":
                run = run.replace("_1", "")
            file = f"{run}_virus_{config.alignmentSoftwareVirus}"
            files = glob.glob(f"{bamDir}/{ref}/{file}.bam")
            if len(files) == 0:
                log.error(f"Sorted BAM file with the alignment for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            log.info(f"Running SNV calling with JACUSA: {run} vs {ref}...")
            cmd = f"java -jar {config.jacusaPath} call-1 -p {config.threads} -r {callingDir}/{ref}/{run}.vcf -a {config.jacusaFeatureFilter} -s -f V{config.jacusa} {bamDir}/{ref}/{file}.bam"
            util.execCmd(cmd)

def filterAS_StrandOddsRatio(sras):
    indexDir = config.workPath + "/indices"
    depthsDir = config.workPath + "/depths"
    util.makeDirectory(depthsDir)
    bamDir = config.workPath + "/bam"
    for fullRef in config.virusReferencePath:
        cmd = f"cp {fullRef} {indexDir}"
        util.execCmd(cmd)
        ref = fullRef[fullRef.rfind("/") + 1:]
        ref = ref[:ref.rfind(".")]
        util.makeDirectory(f"{depthsDir}/{ref}")
        for sra in sras:
            run = sra[0][0]
            run = run[:run.rfind(".")]
            type = sra[1]
            if type == "paired":
                run = run.replace("_1", "")
            file = f"{run}_virus_{config.alignmentSoftwareVirus}"
            files = glob.glob(f"{bamDir}/{ref}/{file}.bam")
            if len(files) == 0:
                log.error(f"Sorted BAM file with the alignment for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            cmd = f"{config.bcftoolsPath} mpileup -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP -O v -A -C -I -d 1000000{config.bcftools} -f {fullRef} -o {depthsDir}/{ref}/{run}.vcf {bamDir}/{ref}/{file}.bam"
            util.execCmd(cmd)
            mpileup = _readMpileupVcf(f"{depthsDir}/{ref}/{run}.vcf")
            mpileup[["PL", "DP", "SP", "ADF", "ADR", "AD"]] = mpileup["FORMAT2"].str.split(":", expand=True)
            mpileup[["fwdRefDepth", "fwdAltDepth"]] = mpileup["ADF"].str.split(",", expand=True).apply(lambda x: x[:2], axis=1).apply(pd.to_numeric)
            mpileup[["revRefDepth", "revAltDepth"]] = mpileup["ADR"].str.split(",", expand=True).apply(lambda x: x[:2], axis=1).apply(pd.to_numeric)
            mpileup["fwdRefDepth"] = mpileup["fwdRefDepth"] + 1
            mpileup["fwdAltDepth"] = mpileup["fwdAltDepth"] + 1
            mpileup["revRefDepth"] = mpileup["revRefDepth"] + 1
            mpileup["revAltDepth"] = mpileup["revAltDepth"] + 1
            mpileup["R"] = (mpileup["fwdRefDepth"] * mpileup["revAltDepth"]) / (mpileup["fwdAltDepth"] * mpileup["revRefDepth"])
            mpileup["sym"] =  mpileup["R"] + (1 / mpileup["R"])
            mpileup["refRatio"] = mpileup[["fwdRefDepth", "revRefDepth"]].min(axis=1) / mpileup[["fwdRefDepth", "revRefDepth"]].max(axis=1)
            mpileup["altRatio"] = mpileup[["fwdAltDepth", "revAltDepth"]].min(axis=1) / mpileup[["fwdAltDepth", "revAltDepth"]].max(axis=1)
            mpileup["SOR"] =  np.log(mpileup["sym"]) + np.log(mpileup["refRatio"]) - np.log(mpileup["altRatio"])
            mpileup.to_csv(f"{depthsDir}/{ref}/{run}.csv", index = False)
            mpileup = mpileup[mpileup["SOR"] <= config.maxAS_StrandOddsRatio]
            mpileup.to_csv(f"{depthsDir}/{ref}/{run}_filtered.csv", index = False)

def _readMpileupVcf(path):
    lines = []
    with open(path, 'r') as f:
        start = False
        for line in f:
            if start:
                lines.append(line.strip())
            if line.startswith("#CHROM"):
                lines.append(line.strip())
                start = True
    csv = pd.read_csv(
        io.StringIO('\n'.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM', "REF": "Reference", "POS": "Position"})
    csv.columns = [*csv.columns[:-1], 'FORMAT2']
    return csv