"""This module handles everything about the alignment .
"""

from operator import countOf
import util
import config
import logger
import glob
import pathlib
import math
from Bio import SeqIO

log = logger.logger
alignmentDir = config.workPath + "/2-alignment"

def runHisat2(sras, host):
    """Runs Hisat2 aligner.
    
    It creates the database for each reference genome, then aligns each
    run against each genome. Each run-genome pair generates a SAM file,
    and it is saved in the 2-alignment/<host/pathogen>/sam directory.
    It can take the arguments listed in the hisat2.config file located
    within the config directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID      
        host (bool): Switch that tells if the alignment must be done against the host or the viral genomes.

    """
    alignmentSubDir = alignmentDir + "/" + ("host" if host else "pathogen")
    indexDir = alignmentSubDir + "/indices/hisat2"
    util.makeDirectory(indexDir)
    samDir = alignmentSubDir + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.pathogenReferenceGenomePaths
    jobsRef = {}
    jobs = []
    for fullRef in refs:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        log.info(f"Building Hisat2 index for {ref}...")
        util.makeDirectory(f"{samDir}/{ref}")
        files = glob.glob(f"{indexDir}/{ref}_hisat2.done")
        if len(files) > 0:
            log.info(f"Hisat2 index for {ref} already built.")
            jobsRef[ref] = ""
        else:
            params = config.hisatIndexH if host else config.hisatIndexV
            cmd = f"{config.hisat2Path}hisat2-build -p {config.threads}{params} {fullRef} {indexDir}/{ref}"
            jobsRef[ref] = util.runCommand(cmd, jobName="hisat2build", outFile=f"{indexDir}/{ref}_hisat2.done")
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} and {pathlib.Path(sra[0][1]).name} to the reference file...")
                files = f"-1 {sra[0][0]} -2 {sra[0][1]}"
            else:
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} to the reference file...")
                files = f"-U {sra[0][0]}"
            sam = sra[2] + "_hisat2.sam"
            params = config.hisatMappingH if host else config.hisatMappingV
            cmd = f"{config.hisat2Path}hisat2 -x {indexDir}/{ref} {files} -S {samDir}/{ref}/{sam} -p {config.threads}{params}"
            util.runCommand(cmd, jobName="hisat2", jobs=jobs, dep=jobsRef[ref])
    util.waitForJobs(jobs)

def runBWA(sras, host):
    """Runs BWA aligner.
    
    It creates the index for each reference genome, then aligns each
    run against each genome. Each run-genome pair generates a SAM file,
    and it is saved in the 2-alignment/<host/pathogen>/sam directory.
    It can take the arguments listed in the bwa.config file located
    within the config directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
        host (bool): Switch that tells if the alignment must be done against the host or the viral genomes.
    """
    alignmentSubDir = alignmentDir + "/" + ("host" if host else "pathogen")
    indexDir = alignmentSubDir + "/indices/bwa"
    util.makeDirectory(indexDir)
    samDir = alignmentSubDir + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.pathogenReferenceGenomePaths
    jobsRef = {}
    jobs = []
    for fullRef in refs:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        log.info(f"Building BWA index file for {ref}...")
        util.makeDirectory(f"{samDir}/{ref}")
        files = glob.glob(f"{indexDir}/{ref}_bwa.done")
        if len(files) > 0:
            log.info(f"BWA index for {ref} already built.")
            jobsRef[ref] = ""
        else:
            params = config.bwaIndexH if host else config.bwaIndexV
            cmd = f"cp {fullRef} {indexDir}/{ref}"
            util.execCmd(cmd)
            cmd = f"{config.bwaPath} index {params} {indexDir}/{ref}"
            jobsRef[ref] = util.runCommand(cmd, jobName="bwaIndex", outFile=f"{indexDir}/{ref}_bwa.done")
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} and {pathlib.Path(sra[0][1]).name} to the reference file...")
                files = f"{sra[0][0]} {sra[0][1]}"
            else:
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} to the reference file...")
                files = f"{sra[0][0]}"
            sam = sra[2] + "_bwa.sam"
            params = config.bwaMappingH if host else config.bwaMappingV
            cmd = f"{config.bwaPath} mem -t {config.threads} {params} {indexDir}/{ref} {files}"
            util.runCommand(cmd, jobName="bwa", jobs=jobs, dep=jobsRef[ref], outFile=f"{samDir}/{ref}/{sam}")
    util.waitForJobs(jobs)

def runSTAR(sras, host):
    """Runs STAR aligner.
    
    It creates the database for each reference genome, then aligns each
    run against each genome. Each run-genome pair generates a SAM file
    with the sequences that may have a relevant alignment, which is
    saved in the 2-alignment/<host/pathogen>/sam directory, and one or two 
    FASTQ files (depending on the run type) with the unaligned 
    sequences, which are saved in the 2-alignment/<host/pathogen>/star 
    directory with .mate1 and .mate2 extensions. It can take the 
    arguments listed in the star.config file located within the config 
    directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
        host (bool): Switch that tells if the alignment must be done against the host or the viral genomes.
    """
    alignmentSubDir = alignmentDir + "/" + ("host" if host else "pathogen")
    indexDir = alignmentSubDir + "/indices"
    util.makeDirectory(indexDir)
    starDir = alignmentSubDir + "/star"
    util.makeDirectory(starDir)
    samDir = alignmentSubDir + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.pathogenReferenceGenomePaths
    jobsRef = {}
    jobs = [] 
    for fullRef in refs:
        bases = 14
        for rec in SeqIO.parse(fullRef, "fasta"):
            bases = len(rec.seq)
        bases = min(14, int(math.log(bases, 2)/2 - 1))
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        log.info(f"Building STAR index file for {ref}...")
        util.makeDirectory(f"{indexDir}/{ref}_star")
        util.makeDirectory(f"{samDir}/{ref}")
        files = glob.glob(f"{indexDir}/{ref}_star.done")
        if len(files) > 0:
            log.info(f"STAR index for {ref} already built.")
            jobsRef[ref] = ""
        else:
            params = config.starIndexH if host else config.starIndexV
            cmd = f"{config.starPath} --runThreadN {config.threads} --runMode genomeGenerate --genomeSAindexNbases {bases} --genomeDir {indexDir}/{ref}_star --genomeFastaFiles {fullRef} {params}"
            jobsRef[ref] = util.runCommand(cmd, jobName="starBuild", outFile=f"{indexDir}/{ref}_star.done")
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} and {pathlib.Path(sra[0][1]).name} to the reference file...")
                files = f"{sra[0][0]} {sra[0][1]}"
            else:
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} to the reference file...")
                files = f"{sra[0][0]}"
            prefix = sra[2]
            params = config.starMappingH if host else config.starMappingV
            cmd = f"{config.starPath} --runThreadN {config.threads} --runMode alignReads --readFilesIn {files} --genomeDir {indexDir}/{ref}_star --outFileNamePrefix {starDir}/{ref}/{prefix} --outSAMunmapped Within {params}"
            cmd += f" && mv {starDir}/{ref}/{prefix}Aligned.out.sam {samDir}/{ref}/{prefix}_{config.alignmentSoftwareHost}.sam"
            util.runCommand(cmd, jobName="star", jobs=jobs, dep=jobsRef[ref])
    util.waitForJobs(jobs)

def runMagicBlast(sras, host):
    """Runs Magic-BLAST aligner.
    
    It creates the database for each reference genome, then aligns each
    run against each genome. Each run-genome pair generates a SAM file,
    and it is saved in the 2-alignment/<host/pathogen>/sam directory.
    It can take the arguments listed in the magicblast.config file 
    located within the config directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
        host (bool): Switch that tells if the alignment must be done against the host or the viral genomes.
    """
    alignmentSubDir = alignmentDir + "/" + ("host" if host else "pathogen")
    indexDir = alignmentSubDir + "/indices/magicblast"
    util.makeDirectory(indexDir)
    samDir = alignmentSubDir + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.pathogenReferenceGenomePaths
    jobsRef = {}
    jobs = []
    for fullRef in refs:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        log.info(f"Building Magic-BLAST index file for {ref}...")
        util.makeDirectory(f"{samDir}/{ref}")
        files = glob.glob(f"{indexDir}/{ref}_magicblast.done")
        if len(files) > 0:
            log.info(f"Magic-BLAST index for {ref} already built.")
            jobsRef[ref] = ""
        else:
            params = config.magicblastIndexH if host else config.magicblastIndexV
            cmd = f"{config.magicblastPath}/makeblastdb -in {fullRef} -out {indexDir}/{ref} -parse_seqids -dbtype nucl {params}"
            jobsRef[ref] = util.runCommand(cmd, jobName="magicblastdb", outFile=f"{indexDir}/{ref}_magicblast.done")
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} and {pathlib.Path(sra[0][1]).name} to the reference file...")
                files = f"-paired -query {sra[0][0]} -query_mate {sra[0][1]}"
            else:
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} to the reference file...")
                files = f"-query {sra[0][0]}"
            sam = sra[2] + "_magicblast.sam"
            params = config.magicblastMappingH if host else config.magicblastMappingV
            cmd = f"{config.magicblastPath}/magicblast -db {indexDir}/{ref} {files} -out {samDir}/{ref}/{sam} -infmt fastq -num_threads {config.threads} {params}"
            util.runCommand(cmd, jobName="magicblast", jobs=jobs, dep=jobsRef[ref])
    util.waitForJobs(jobs)

def runMinimap2(sras, host):
    """Runs Minimap2 aligner.
    
    It creates the database for each reference genome, then aligns each
    run against each genome. Each run-genome pair generates a SAM file,
    and it is saved in the 2-alignment/<host/pathogen>/sam directory.
    It can take the arguments listed in the minimap2.config file 
    located within the config directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
        host (bool): Switch that tells if the alignment must be done against the host or the viral genomes.
    """
    alignmentSubDir = alignmentDir + "/" + ("host" if host else "pathogen")
    indexDir = alignmentSubDir + "/indices/minimap2"
    util.makeDirectory(indexDir)
    samDir = alignmentSubDir + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.pathogenReferenceGenomePaths
    preset = config.minimapPresetH if host else config.minimapPresetV
    jobsRef = {}
    jobs = []
    for fullRef in refs:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        log.info(f"Building Minimap2 index files for {ref}...")
        util.makeDirectory(f"{samDir}/{ref}")
        files = glob.glob(f"{indexDir}/{ref}_minimap2.done")
        if len(files) > 0:
            log.info(f"Minimap2 index for {ref} already built.")
            jobsRef[ref] = ""
        else:
            params = config.minimapIndexH if host else config.minimapIndexV
            cmd = f"{config.minimapPath} -x {preset} -t {config.threads} {params} -d {indexDir}/{ref}.mmi {fullRef}"
            jobsRef[ref] = util.runCommand(cmd, jobName="minimap2build", outFile=f"{indexDir}/{ref}_minimap2.done")
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} and {pathlib.Path(sra[0][1]).name} to the reference file...")
                files = f"{sra[0][0]} {sra[0][1]}"
            else:
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} to the reference file...")
                files = f"{sra[0][0]}"
            sam = sra[2] + "_minimap2.sam"
            params = config.minimapMappingH if host else config.minimapMappingV
            cmd = f"{config.minimapPath} -a -x {preset} {params} {indexDir}/{ref}.mmi {files}"
            util.runCommand(cmd, jobName="minimap2", jobs=jobs, dep=jobsRef[ref], outFile=f"{samDir}/{ref}/{sam}")
    util.waitForJobs(jobs)

def runGMAP(sras, host):
    """Runs GMAP aligner.
    
    It creates the database for each reference genome, then aligns each
    run against each genome. Each run-genome pair generates a SAM file,
    and it is saved in the 2-alignment/<host/pathogen>/sam directory.
    It can take the arguments listed in the gmap.config file located
    within the config directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
        host (bool): Switch that tells if the alignment must be done against the host or the viral genomes.
    """
    alignmentSubDir = alignmentDir + "/" + ("host" if host else "pathogen")
    indexDir = alignmentSubDir + "/indices/gmap" 
    util.makeDirectory(indexDir)
    samDir = alignmentSubDir + "/sam"
    util.makeDirectory(samDir)
    refs = config.hostReferencePath if host else config.pathogenReferenceGenomePaths
    jobsRef = {}
    jobs = []
    for fullRef in refs:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        log.info(f"Building GMAP index files for {ref}...")
        util.makeDirectory(f"{samDir}/{ref}")
        files = glob.glob(f"{indexDir}/{ref}_gmap.done")
        if len(files) > 0:
            log.info(f"GMAP index for {ref} already built.")
            jobsRef[ref] = ""
        else:
            params = config.gmapIndexH if host else config.gmapIndexV
            cmd = f"{config.gmapPath}gmap_build {params} -D {indexDir} -d {ref} {fullRef}"
            jobsRef[ref] = util.runCommand(cmd, jobName="gmapBuild", outFile=f"{indexDir}/{ref}_gmap.done")
        for sra in sras:
            files = ""
            if sra[1] == "paired":
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} and {pathlib.Path(sra[0][1]).name} to the reference file...")
                files = f"sampe {sra[0][0]} {sra[0][1]}"
            else:
                log.info(f"Mapping {pathlib.Path(sra[0][0]).name} to the reference file...")
                files = f"samse {sra[0][0]}"
            sam = sra[2] + f"_gmap.sam"
            params = config.gmapMappingH if host else config.gmapMappingV
            cmd = f"{config.gmapPath}gmap {params} -t {config.threads} -D {indexDir} -d {ref} -f {files}"
            util.runCommand(cmd, jobName="gmap", jobs=jobs, dep=jobsRef[ref], outFile=f"{samDir}/{ref}/{sam}")
    util.waitForJobs(jobs)

def extractUnaligned(sras):
    """Extracts the unaligned sequences from the alignment against the 
    host and discards the runs where all sequences aligned.

    It filters the SAM files found in the 2-alignment/host/sam
    directory with samtools view using the filter 256 and these
    results are transformed into BAM files, found in the
    2-alignment/host/bam directory. Then, using samtools bam2fq, it 
    transforms the BAM files into FASTQ files with the final unaligned 
    sequences. 
    
    In the case of paired-end runs, the resulting FASTQ has both 3' and
    5' reads in the same file, so they are split into two files.

    In the case of the unaligned FASTQ files generated by STAR, these 
    are appended to the split FASTQ files.

    The resulting FASTQ files are saved at 2-alignment/host/fastq.

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
    alignmentSubDir = alignmentDir + "/host"
    bamDir = alignmentSubDir + "/bam"
    util.makeDirectory(bamDir)
    samDir = alignmentSubDir + "/sam"
    fastqDir = alignmentSubDir + "/fastq"
    util.makeDirectory(fastqDir)
    newsras = []
    jobs = []
    for sra in sras:
        run = sra[2]
        path = pathlib.Path(config.hostReferencePath[0])
        ref = path.parent.name
        util.makeDirectory(bamDir + "/" + ref)
        print(f"{samDir}/{ref}/{run}_{config.alignmentSoftwareHost}.sam")
        files = glob.glob(f"{samDir}/{ref}/{run}_{config.alignmentSoftwareHost}.sam")
        if len(files) == 0:
            log.error(f"Alignment against host for run with ID {run} not found.")
            util.stopProgram()
        includeFlag = 4
        excludeNotPaired = ""
        if sra[1] == "paired":
            includeFlag = 13
            # excludeNotPaired = " -F 12"

        cmd = f"{config.samtoolsPath} view -f {includeFlag} -F 256 {excludeNotPaired} -bo {bamDir}/{ref}/{run}_{config.alignmentSoftwareHost}.bam {samDir}/{ref}/{run}_{config.alignmentSoftwareHost}.sam"
        cmd += f" && {config.samtoolsPath} bam2fq {bamDir}/{ref}/{run}_{config.alignmentSoftwareHost}.bam"
        if sra[1] == "paired":
            util.runCommand(cmd, jobName="getUnaligned", jobs=jobs, outFile=f"{fastqDir}/{sra[2]}.unmapped.fastq")
        else:
            util.runCommand(cmd, jobName="getUnaligned", jobs=jobs, outFile=f"{fastqDir}/{sra[2]}.fastq")
    util.waitForJobs(jobs)

    empty = {}
    ref = pathlib.Path(config.hostReferencePath[0]).parent.name
    for sra in sras:
        run = sra[2]
        empty[run] = True

    for sra in sras:
        run = sra[2]
        if sra[1] == "paired":
            file1 = []
            file2 = []
            count = 0
            for rec in SeqIO.parse(f"{fastqDir}/{sra[2]}.unmapped.fastq", "fastq"):
                count += 1
                if rec.id.endswith("/1"):
                    file1.append(rec)
                elif rec.id.endswith("/2"):
                    file2.append(rec)
            mode = "w"
            if config.alignmentSoftwareHost == "star":
                mode = "a"
            with open(f"{fastqDir}/{run}_1.fastq", mode) as f:
                SeqIO.write(file1, f, "fastq")
            with open(f"{fastqDir}/{run}_2.fastq", mode) as f:
                SeqIO.write(file2, f, "fastq")
            if count > 0:
                empty[run] = False
            cmd = f"rm {fastqDir}/{sra[2]}.unmapped.fastq"
            util.execCmd(cmd)
        else:
            count = 0
            for rec in SeqIO.parse(f"{fastqDir}/{sra[2]}.fastq", "fastq"):
                empty[run] = False
                break
        if not empty[run]:
            newsras.append(sra)
        else:
            log.info(f"All sequences from run with ID {run} aligned to the host. Discarding file...")
    
    for i in range(len(newsras)):
        files = []
        if newsras[i][1] == "single":
            files.append(f"{fastqDir}/{newsras[i][2]}.fastq")
        else:
            files.append(f"{fastqDir}/{newsras[i][2]}_1.fastq")
            files.append(f"{fastqDir}/{newsras[i][2]}_2.fastq")
        newsras[i] = (files, newsras[i][1], newsras[i][2])
    return newsras

def sortAlignments(sras):
    """Sorts the SAM files of each alignment against the viral genomes.

    It sorts and transforms these SAM files into BAM files using the 
    samtools sort tool, and then adds some extra headers needed for
    following steps using samtools calmd. These files are saved in the 
    2-alignment/pathogen/bam directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    alignmentSubDir = config.workPath + "/2-alignment/pathogen"
    samDir = alignmentSubDir + "/sam"
    bamDir = alignmentSubDir + "/bam"
    util.makeDirectory(bamDir)
    jobs = []
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name 
        util.makeDirectory(f"{bamDir}/{ref}")
        for sra in sras:
            run = sra[2]
            files = glob.glob(f"{samDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.sam")
            if len(files) == 0:
                print(f"{samDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.sam")
                log.error(f"Alignment for run with ID {run} against pathogen reference {ref} not found.")
                util.stopProgram()
            log.info(f"Sorting run with ID {run} vs {ref} BAM file...")
            cmd = f"{config.samtoolsPath} sort -O BAM -o {bamDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.sorted.bam {samDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.sam"
            jobId = util.runCommand(cmd, jobName="sort", jobs=jobs)
            cmd = f"samtools calmd -b {bamDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.sorted.bam {fullRef}"
            util.runCommand(cmd, jobName="calmd", jobs=jobs, outFile=f"{bamDir}/{ref}/{run}_{config.alignmentSoftwarePathogen}.bam", dep=jobId)
    util.waitForJobs(jobs)
