from operator import countOf
import util
import config
import logger
import glob
from Bio import SeqIO

log = logger.logger

def runHisat2(sras, host):
    fastqDir = config.workPath + "/fastq"
    indexDir = config.workPath + "/indices"
    util.makeDirectory(indexDir)
    samDir = config.workPath + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.virusReferencePath
    for fullRef in refs:
        ref = fullRef[fullRef.rfind("/") + 1:]
        log.info(f"Building Hisat2 index for {ref}...")
        ref = ref[:ref.rfind(".")]
        util.makeDirectory(f"{samDir}/{ref}")
        files = glob.glob(f"{indexDir}/{ref}_hisat2.done")
        if len(files) > 0:
            log.info(f"Hisat2 index for {ref} already built.")
        else:
            cmd = f"{config.hisat2Path}hisat2-build -p {config.threads} {fullRef} {indexDir}/{ref}"
            util.execCmd(cmd)
            util.execCmd("echo ''", f"{indexDir}/{ref}_hisat2.done")
        for sra in sras:    
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {sra[0][0]} and {sra[0][1]} to the reference file...")
                files = f"-1 {fastqDir}/{sra[0][0]} -2 {fastqDir}/{sra[0][1]}"
            else:
                log.info(f"Mapping {sra[0][0]} to the reference file...")
                files = f"-U {fastqDir}/{sra[0][0]}"
            suffix = f"_{'host' if host else 'virus'}_hisat2.sam"
            sam = sra[0][0].replace("_1.fastq", suffix).replace("_2.fastq", suffix).replace(".fastq", suffix)
            params = config.hisatH if host else config.hisatV
            cmd = f"{config.hisat2Path}hisat2 -x {indexDir}/{ref} {files} -S {samDir}/{ref}/{sam} -p {config.threads}{params}"
            util.execCmd(cmd)

def runBWA(sras, host):
    fastqDir = config.workPath + "/fastq"
    indexDir = config.workPath + "/indices"
    util.makeDirectory(indexDir)
    samDir = config.workPath + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.virusReferencePath
    params = config.bwaIndexH if host else config.bwaIndexV
    for fullRef in refs:
        ref = fullRef[fullRef.rfind("/") + 1:]
        log.info(f"Building BWA index file for {ref}...")
        shortRefName = ref[:ref.rfind(".")]
        util.makeDirectory(f"{samDir}/{shortRefName}")
        refName = ref.split("/")[-1]
        cmd = f"cp {fullRef} {indexDir}/{refName}"
        util.execCmd(cmd)
        cmd = f"{config.bwaPath} index{params} {indexDir}/{refName}"
        util.execCmd(cmd)
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {sra[0][0]} and {sra[0][1]} to the reference file...")
                files = f"{fastqDir}/{sra[0][0]} {fastqDir}/{sra[0][1]}"
            else:
                log.info(f"Mapping {sra[0][0]} to the reference file...")
                files = f"{fastqDir}/{sra[0][0]}"
            suffix = f"_{'host' if host else 'virus'}_bwa.sam"
            sam = sra[0][0].replace("_1.fastq", suffix).replace("_2.fastq", suffix).replace(".fastq", suffix)
            params = config.bwaMappingH if host else config.bwaMappingV
            cmd = f"{config.bwaPath} mem -t {config.threads}{params} {indexDir}/{refName} {files}"
            util.execCmd(cmd, f"{samDir}/{shortRefName}/{sam}")

def runSTAR(sras, host):
    fastqDir = config.workPath + "/fastq"
    indexDir = config.workPath + "/indices"
    util.makeDirectory(indexDir)
    starDir = config.workPath + "/star"
    util.makeDirectory(starDir)
    hostVirus = 'host' if host else 'virus'
    util.makeDirectory(f"{starDir}/{hostVirus}")
    samDir = config.workPath + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.virusReferencePath
    params = config.starIndexH if host else config.starIndexV
    for fullRef in refs:
        ref = fullRef[fullRef.rfind("/") + 1:]
        log.info(f"Building STAR index file for {ref}...")
        ref = ref[:ref.rfind(".")]
        util.makeDirectory(f"{indexDir}/{ref}")
        util.makeDirectory(f"{samDir}/{ref}")
        cmd = f"{config.starPath} --runThreadN {config.threads} --runMode genomeGenerate --genomeDir {indexDir}/{ref} --genomeFastaFiles {fullRef}{params}"
        util.execCmd(cmd)
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {sra[0][0]} and {sra[0][1]} to the reference file...")
                files = f"{fastqDir}/{sra[0][0]} {fastqDir}/{sra[0][1]}"
            else:
                log.info(f"Mapping {sra[0][0]} to the reference file...")
                files = f"{fastqDir}/{sra[0][0]}"
            prefix = sra[0][0].replace("_1.fastq", "").replace("_2.fastq", "").replace(".fastq", "")
            params = " --outReadsUnmapped Fastx" + config.starMappingH if host else config.starMappingV
            cmd = f"{config.starPath} --runThreadN {config.threads} --runMode alignReads --readFilesIn {files} --genomeDir {indexDir}/{ref} --outFileNamePrefix {starDir}/{hostVirus}/{prefix}{params}"
            util.execCmd(cmd)
            cmd = f"mv {fastqDir}/{sra[0][0]} {fastqDir}/{sra[0][0]}.prealign"
            util.execCmd(cmd)
            cmd = f"mv {starDir}/{hostVirus}/{prefix}Unmapped.out.mate1 {fastqDir}/{sra[0][0]}"
            util.execCmd(cmd)
            if sra[1] == "paired":
                cmd = f"mv {fastqDir}/{sra[0][1]} {fastqDir}/{sra[0][1]}.prealign"
                util.execCmd(cmd)
                cmd = f"mv {starDir}/{hostVirus}/{prefix}Unmapped.out.mate2 {fastqDir}/{sra[0][1]}"
                util.execCmd(cmd)

def runMagicBlast(sras, host):
    fastqDir = config.workPath + "/fastq"
    indexDir = config.workPath + "/db"
    util.makeDirectory(indexDir)
    samDir = config.workPath + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.virusReferencePath
    for fullRef in refs:
        ref = fullRef[fullRef.rfind("/") + 1:]
        log.info(f"Building Magic-BLAST index file for {ref}...")
        ref = ref[:ref.rfind(".")]
        util.makeDirectory(f"{samDir}/{ref}")
        cmd = f"{config.magicblastPath}makeblastdb -in {fullRef} -out {indexDir}/{ref} -parse_seqids -dbtype nucl"
        util.execCmd(cmd)
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {sra[0][0]} and {sra[0][1]} to the reference file...")
                files = f"-paired -query {fastqDir}/{sra[0][0]} -query_mate {fastqDir}/{sra[0][1]}"
            else:
                log.info(f"Mapping {sra[0][0]} to the reference file...")
                files = f"-query {fastqDir}/{sra[0][0]}"
            suffix = f"_{'host' if host else 'virus'}_magicblast.sam"
            sam = sra[0][0].replace("_1.fastq", suffix).replace("_2.fastq", suffix).replace(".fastq", suffix)
            params = config.magicblastH if host else config.magicblastV
            cmd = f"{config.magicblastPath}magicblast -db {indexDir}/{ref} {files} -out {samDir}/{ref}/{sam} -infmt fastq -num_threads {config.threads}{params}"
            util.execCmd(cmd)

def runMinimap2(sras, host):
    fastqDir = config.workPath + "/fastq"
    indexDir = config.workPath + "/indices"
    util.makeDirectory(indexDir)
    samDir = config.workPath + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.virusReferencePath
    preset = config.minimapPresetH if host else config.minimapPresetV
    params = config.minimapIndexH if host else config.minimapIndexV
    for fullRef in refs:
        ref = fullRef[fullRef.rfind("/") + 1:]
        log.info(f"Building Minimap2 index files for {ref}...")
        ref = ref[:ref.rfind(".")]
        util.makeDirectory(f"{samDir}/{ref}")
        cmd = f"{config.minimapPath} -x {preset} -t {config.threads}{params} -d {indexDir}/{ref}.mmi {fullRef}"
        util.execCmd(cmd)
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {sra[0][0]} and {sra[0][1]} to the reference file...")
                files = f"{fastqDir}/{sra[0][0]} {fastqDir}/{sra[0][1]}"
            else:
                log.info(f"Mapping {sra[0][0]} to the reference file...")
                files = f"{fastqDir}/{sra[0][0]}"
            suffix = f"_{'host' if host else 'virus'}_minimap2.sam"
            sam = sra[0][0].replace("_1.fastq", suffix).replace("_2.fastq", suffix).replace(".fastq", suffix)
            params = config.minimapMappingH if host else config.minimapMappingV
            cmd = f"{config.minimapPath} -a -x {preset}{params} {indexDir}/{ref}.mmi {files}"
            util.execCmd(cmd, f"{samDir}/{ref}/{sam}")

def runGMAP(sras, host):
    fastqDir = config.workPath + "/fastq"
    indexDir = config.workPath + "/indices"
    util.makeDirectory(indexDir)
    samDir = config.workPath + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.virusReferencePath
    params = config.gmapIndexH if host else config.gmapIndexV
    for fullRef in refs:
        ref = fullRef[fullRef.rfind("/") + 1:]
        log.info(f"Building GMAP index files for {ref}...")
        ref = ref[:ref.rfind(".")]
        util.makeDirectory(f"{samDir}/{ref}")
        cmd = f"{config.gmapPath}gmap_build{params} -D {indexDir} -d {ref} {fullRef}"
        util.execCmd(cmd)
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {sra[0][0]} and {sra[0][1]} to the reference file...")
                files = f"sampe {fastqDir}/{sra[0][0]} {fastqDir}/{sra[0][1]}"
            else:
                log.info(f"Mapping {sra[0][0]} to the reference file...")
                files = f"samse {fastqDir}/{sra[0][0]}"
            suffix = f"_{'host' if host else 'virus'}_gmap.sam"
            sam = sra[0][0].replace("_1.fastq", suffix).replace("_2.fastq", suffix).replace(".fastq", suffix)
            params = config.gmapMappingH if host else config.gmapMappingV
            cmd = f"{config.gmapPath}gmap {params} -D {indexDir} -d {ref} -f {files}"
            util.execCmd(cmd, f"{samDir}/{ref}/{sam}")

def extractUnaligned(sras):
    if config.alignmentSoftwareHost == "star":
        log.info(f"Alignment done with STAR. Skipping extraction...")
        return sras
    fastqDir = config.workPath + "/fastq"
    bamDir = config.workPath + "/bam"
    util.makeDirectory(bamDir)
    newsras = []
    for sra in sras:
        seqs = set()
        run = sra[0][0]
        run = run[:run.rfind(".")]
        runName = run.replace("_1", "")
        ref = config.hostReferencePath[0]
        ref = ref[ref.rfind("/") + 1:]
        ref = ref[:ref.rfind(".")]
        util.makeDirectory(bamDir + "/" + ref)
        files = glob.glob(f"{config.workPath}/sam/{ref}/{runName}_host_{config.alignmentSoftwareHost}.sam")
        if len(files) == 0:
            log.error(f"Alignment against host for run {runName} not found.")
            util.stopProgram()
        includeFlag = 4
        if sra[1] == "paired":
            includeFlag = 13
        cmd = f"{config.samtoolsPath} view -f {includeFlag} -F 256 -bo {config.workPath}/bam/{ref}/{runName}_host_{config.alignmentSoftwareHost}.bam {config.workPath}/sam/{ref}/{runName}_host_{config.alignmentSoftwareHost}.sam"
        util.execCmd(cmd)
        for file in sra[0]:
            cmd = f"mv {config.workPath}/fastq/{file} {config.workPath}/fastq/{file}.prealign"
            util.execCmd(cmd)
        cmd = f"{config.samtoolsPath} bam2fq {config.workPath}/bam/{ref}/{runName}_host_{config.alignmentSoftwareHost}.bam"
        if sra[1] == "paired":
            util.execCmd(cmd, f"{fastqDir}/{runName}.unmapped.fastq")
            file1 = []
            file2 = []
            count = 0
            for rec in SeqIO.parse(f"{fastqDir}/{runName}.unmapped.fastq", "fastq"):
                count += 1
                if rec.id.endswith("/1"):
                    file1.append(rec)
                elif rec.id.endswith("/2"):
                    file2.append(rec)
            with open(f"{fastqDir}/{runName}_1.fastq", "w") as f:
                SeqIO.write(file1, f, "fastq")
            with open(f"{fastqDir}/{runName}_2.fastq", "w") as f:
                SeqIO.write(file2, f, "fastq")
            if count > 0:
                newsras.append(sra)
            else:
                log.info(f"All sequences from {runName} aligned to the host. Discarding file...")
        else:
            util.execCmd(cmd, f"{fastqDir}/{runName}.fastq")
            count = 0
            for rec in SeqIO.parse(f"{fastqDir}/{runName}.fastq", "fastq"):
                count += 1
            if count > 0:
                newsras.append(sra)
            else:
                log.info(f"All sequences from {runName} aligned to the host. Discarding file...")
    return newsras
