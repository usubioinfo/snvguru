"""This module handles everything about the SNV calling.
"""

import util
import config
import logger
import io
import pandas as pd
import numpy as np
import glob 
import pathlib
import dask.dataframe as dd

log = logger.logger
baseDir = config.workPath + "/4-snvCalling"

def runReditools2(sras):
    """Runs REDItools2.

    It creates an index for each reference genome using samtools faidx 
    and runs the calling process for each run BAM file. The resulting 
    files' names end with .calling.txt, which are saved in the 
    4-snvCalling/calling directory. Finally, these files are 
    transformed into VCF files. These files' names end with 
    .reditools.vcf, which are saved in the 
    4-snvCalling/calling directory

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    indexDir = baseDir + "/indices"
    homopolDir = baseDir + "/homopol"
    util.makeDirectory(homopolDir)
    callingDir = baseDir + "/calling/reditools"
    util.makeDirectory(callingDir)
    bamDir = config.workPath + "/2-alignment/pathogen/bam"
    jobs = []
    jobsRef = {}
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name      
        util.makeDirectory(f"{indexDir}/{ref}")
        cmd = f"cp {fullRef} {indexDir}/{ref}/"
        util.execCmd(cmd)
        util.makeDirectory(f"{homopolDir}/{ref}")
        util.makeDirectory(f"{callingDir}/{ref}")
        log.info(f"Building Samtools index file for {ref}...")
        files = glob.glob(f"{indexDir}/{ref}/{ref}_faidx.done")
        if len(files) > 0:
            log.info(f"Samtools index for {ref} already built.")
            jobsRef[ref] = ""
        else:
            cmd = f"{config.samtoolsPath} faidx {indexDir}/{ref}/genome.fa"
            jobsRef[ref] = util.runCommand(cmd, jobName="faidx", outFile=f"{indexDir}/{ref}/{ref}_faidx.done")
        for sra in sras:
            run = sra[2]
            f = f"{run}_{config.alignmentSoftwarePathogen}"
            files = glob.glob(f"{bamDir}/{ref}/{f}.bam")
            if len(files) == 0:
                log.error(f"Sorted BAM file with the alignment for run with ID {run} against pathogen reference {ref} not found.")
                util.stopProgram()
            log.info(f"Running SNV calling with Reditools2: {run} vs {ref}...")
            cmd = f"{config.samtoolsPath} index {bamDir}/{ref}/{f}.bam"
            util.runCommand(cmd, jobName="index", jobs=jobs, dep=jobsRef[ref])
    util.waitForJobs(jobs)
    jobs = []
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        for sra in sras:
            run = sra[2]
            cmd = f"{config.reditoolsCommand} -c -f {bamDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.bam -r {indexDir}/{ref}/genome.fa -m {homopolDir}/{ref}/{run}.homopol.txt -o {callingDir}/{ref}/{run}.reditools.txt{config.reditools}"
            util.runCommand(cmd, jobName="reditools2", jobs=jobs)
    util.waitForJobs(jobs)
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        for sra in sras:
            run = sra[2]
            _reditoolsToVcf(f"{callingDir}/{ref}/{run}.reditools.txt", run)

def runJacusa(sras):
    """Runs JACUSA.

    It runs the calling process for each run BAM file. The resulting 
    files' names end with .jacusa.vcf, and are found in the 
    4-snvCalling/calling directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID   
    """
    indexDir = baseDir + "/indices"
    callingDir = baseDir + "/calling/jacusa"
    util.makeDirectory(callingDir)
    bamDir = config.workPath + "/2-alignment/pathogen/bam"
    jobs = []
    for fullRef in config.pathogenReferenceGenomePaths:
        cmd = f"cp {fullRef} {indexDir}"
        util.execCmd(cmd)
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        util.makeDirectory(f"{callingDir}/{ref}")
        for sra in sras:
            run = sra[2]
            f = f"{run}_{config.alignmentSoftwarePathogen}"
            files = glob.glob(f"{bamDir}/{ref}/{f}.bam")
            if len(files) == 0:
                log.error(f"Sorted BAM file with the alignment for run with ID {run} against pathogen reference {ref} not found.")
                util.stopProgram()
            log.info(f"Running SNV calling with JACUSA: {run} vs {ref}...")
            cmd = f"{config.javaPath} -jar {config.jacusaPath} call-1 -p {config.threads} -r {callingDir}/{ref}/{run}.jacusa.vcf -s -f V{config.jacusa} {bamDir}/{ref}/{f}.bam"
            util.runCommand(cmd, jobName="jacusa", jobs=jobs)
    util.waitForJobs(jobs)
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        for sra in sras:
            run = sra[2]
            _jacusaToSnpEff(f"{callingDir}/{ref}/{run}.jacusa.vcf", run)

def filterAS_StrandOddsRatio(sras):
    """Calculates the AS_StrandOddsRatio value for each read using
    bcftools mpileup.

    It runs mpileup for each run BAM file. The resulting files'
    names end with .mpileup.vcf, which are saved in the
    4-snvCalling/depths directory. Then, these files are processed
    and filtered, so that only the SNVs positions that pass the filter 
    remain. These results are saved in files whose name end with 
    _filtered.csv, and are saved at 4-snvCalling/depths.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    depthsDir = baseDir + "/depths"
    util.makeDirectory(depthsDir)
    bamDir = config.workPath + "/2-alignment/pathogen/bam"
    jobs = []
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        util.makeDirectory(f"{depthsDir}/{ref}")
        for sra in sras:
            run = sra[2]
            f = f"{run}_{config.alignmentSoftwarePathogen}"
            files = glob.glob(f"{bamDir}/{ref}/{f}.bam")
            if len(files) == 0:
                log.error(f"Sorted BAM file with the alignment for run with ID {run} against pathogen reference {ref} not found.")
                util.stopProgram()
            cmd = f"{config.bcftoolsPath} mpileup{config.bcftools} -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP -O v -f {fullRef} -o {depthsDir}/{ref}/{run}.mpileup.vcf {bamDir}/{ref}/{f}.bam"
            util.runCommand(cmd, jobName="mpileup", jobs=jobs)
    util.waitForJobs(jobs)
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        for sra in sras:
            run = sra[2]
            log.info(f"Analyzing {run}...")
            mpileup = _readMpileupVcf(f"{depthsDir}/{ref}/{run}.mpileup.vcf")
            mpileup = dd.from_pandas(mpileup, chunksize=1000)
            mpileup["FORMAT2"] = mpileup["FORMAT2"].str.partition(":")[2].str.partition(":")[2].str.partition(":")[2]
            mpileup["ADF"] = mpileup["FORMAT2"].str.partition(":")[0]
            mpileup["fwdRefDepth"] = dd.to_numeric(mpileup["ADF"].str.partition(",")[0])
            mpileup["fwdAltDepth"] = dd.to_numeric(mpileup["ADF"].str.partition(",")[2].str.partition(",")[0])
            mpileup["FORMAT2"] = mpileup["FORMAT2"].str.partition(":")[2]
            mpileup["ADR"] = mpileup["FORMAT2"].str.partition(":")[0]
            mpileup["revRefDepth"] = dd.to_numeric(mpileup["ADR"].str.partition(",")[0])
            mpileup["revAltDepth"] = dd.to_numeric(mpileup["ADR"].str.partition(",")[2].str.partition(",")[0])
            mpileup["fwdRefDepth"] = mpileup["fwdRefDepth"] + 1
            mpileup["fwdAltDepth"] = mpileup["fwdAltDepth"] + 1
            mpileup["revRefDepth"] = mpileup["revRefDepth"] + 1
            mpileup["revAltDepth"] = mpileup["revAltDepth"] + 1
            mpileup["R"] = (mpileup["fwdRefDepth"] * mpileup["revAltDepth"]) / (mpileup["fwdAltDepth"] * mpileup["revRefDepth"])
            mpileup["sym"] =  mpileup["R"] + (1 / mpileup["R"])
            mpileup["refRatio"] = mpileup[["fwdRefDepth", "revRefDepth"]].min(axis=1) / mpileup[["fwdRefDepth", "revRefDepth"]].max(axis=1)
            mpileup["altRatio"] = mpileup[["fwdAltDepth", "revAltDepth"]].min(axis=1) / mpileup[["fwdAltDepth", "revAltDepth"]].max(axis=1)
            mpileup["SOR"] =  np.log(mpileup["sym"]) + np.log(mpileup["refRatio"]) - np.log(mpileup["altRatio"])
            mpileup = mpileup[mpileup["SOR"] <= config.maxAS_StrandOddsRatio]
            mpileup = mpileup[["CHROM", "Position"]].drop_duplicates()
            if len(glob.glob(f"{depthsDir}/{ref}/{run}_filtered.hdf")) >= 1:
                util.execCmd(f"rm {depthsDir}/{ref}/{run}_filtered.hdf")
            mpileup.to_hdf(f"{depthsDir}/{ref}/{run}_filtered.hdf", key=run, index = False)

def _calculateFrequencyReditools2(row):
    """Calculates the frequency of the variant in
    each row of the REDItools2 VCF files.

    It divides the support count for the alternate
    allele and divides that by the support count
    for the reference allele + alternate allele.

    Args:
        row (Series): Row of the VCF file.

    Returns:
        frequency (float): Frequency of the variant 
    """
    alt = row["ALT"]
    ct = float(row[alt])
    return ct / (float(row[row["REF"]]) + ct)

def _getVariantReadsReditools2(row):
    """Finds the total of reads for the variant in
    each row of the REDItools2 VCF files.

    Args:
        row (Series): Row of the VCF file.

    Returns:
        variantReads (int): Total of reads for the variant
    """
    alt = row["ALT"]
    return int(row[alt])

def _reditoolsToVcf(path, run):
    """Converts the given REDItools2 VCF file to a VCF format that is 
    accepted by SnpEff.

    The columns #CHROM, POSITION, ID, REF, ALT, QUAL, FILTER and INFO
    remain, while the others are removed. The resulting VCF files'
    names end with .reditools.presnpeff.vcf, and are saved in the 
    4-snvCalling/calling directory.

    Args:
        path (str): Path of the VCF file.
        run (str): Run ID.
    """
    final=[]
    snp =[]
    with open(path, 'r') as fp:
        fp.readline()
        for line in fp:
            data = line.strip().split("\t")
            chrom = data[0]
            pos = data[1]
            ref = data[2]
            if data[3] == 2:
                strand = '+'
            else:
                strand= '-'
            cov = data[4]
            mq = data[5]
            alt = ",".join([s[1] for s in data[7].split(" ")])
            freq = 0.0
            bases = data[6].split("[")[1].split("]")[0].split(",")
            Ab = bases[0]
            Cb = bases[1]
            Gb = bases[2]
            Tb = bases[3]
            if alt == 'A':
                vr = Ab
            if alt == 'C':
                vr = Cb
            if alt == 'G':
                vr = Gb
            if alt == 'T':
                vr = Tb

            readId = '.'
            readFilter = '.'
            info = ' '

            final.append([chrom, pos, readId, ref, alt, mq, readFilter, vr, cov, freq, Ab, Cb, Gb, Tb, info])

        df = pd.DataFrame(final, columns=['#CHROM', 'POSITION','ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'VARIANT_READ', 'TOTAL_READ', 'FREQUENCY', 'A', 'C', 'G', 'T', 'INFO' ])
        df["ALT"] = df["ALT"].str.split(",")
        df = df.explode(["ALT"])
        df["FREQUENCY"] = df.apply(_calculateFrequencyReditools2, axis=1)
        df = df.round({"FREQUENCY": 6})
        df["FREQUENCY"] = df["FREQUENCY"] * 100
        df["VARIANT_READ"] = df.apply(_getVariantReadsReditools2, axis=1)
        snpeff = df[['#CHROM', 'POSITION','ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]
    
    template =  "##fileformat=VCFv4.2\n"\
                "{}"
    p = pathlib.Path(path)
    with open(f"{p.parent}/{run}.reditools.vcf", "w") as f:
        f.write(template.format(df.to_csv(index=False, sep="\t")))
    with open(f"{p.parent}/{run}.reditools.presnpeff.vcf", "w") as f:
        f.write(template.format(snpeff.to_csv(index=False, sep="\t")))

def _jacusaToSnpEff(path, run):
    """Converts the given JACUSA VCF file to a VCF format that is 
    accepted by SnpEff.

    The columns #CHROM, POSITION, ID, REF, ALT, QUAL, FILTER and INFO
    remain, while the others are removed. The resulting VCF files'
    names end with .jacusa.presnpeff.vcf, and are saved in the 
    4-snvCalling/calling directory.

    Args:
        path (str): Path of the VCF file.
        run (str): Run ID.
    """
    util.execCmd(f"mv {path} {path}.pre")
    lines = []
    with open(path + ".pre", "r") as r:
        for line in r:
            if line.startswith("#CHROM"):
                lines.append(line)
                break
        lines.extend([l for l in r])
    template =  "##fileformat=VCFv4.2\n"\
                "{}"
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    )
    df["ALT"] = df["ALT"].str.split(",")
    df = df.explode(["ALT"])
    snpeff = df[['#CHROM', 'POS','ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]
    p = pathlib.Path(path)
    with open(f"{p.parent}/{run}.jacusa.vcf", "w") as f:
        f.write(template.format(df.to_csv(index=False, sep="\t")))
    with open(f"{p.parent}/{run}.jacusa.presnpeff.vcf", "w") as f:
        f.write(template.format(snpeff.to_csv(index=False, sep="\t")))
    util.execCmd(f"rm {path}.pre")

def _readMpileupVcf(path):
    """It reads the given output VCF file of mpileup as a Pandas 
    dataframe.

    Args:
        path (str): Path of the VCF file.

    Returns:
        DataFrame: Pandas dataframe with the data read from the VCF file.
    """
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