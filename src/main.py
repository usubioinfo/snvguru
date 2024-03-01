"""This is the main module.
"""

import glob
import config
import logger
import arguments
import software
import downloads
import quality
import cropping
import alignment
import calling
import results
import report
import util
import os
import joblib
import pathlib

log = logger.logger
path = pathlib.Path(config.originalHostPath[0])
ref = path.stem


config.hostReferencePath = [f"{config.workPath}/2-alignment/host/genomes/{ref}/genome.fa"]

step = config.step

if step not in ["mergeCalling", "generateGraphs", "generateReport"] and config.resumeFrom not in ["mergeCalling", "generateGraphs", "generateReport"]:
    software.checkTools()

sras = []
useFastqInputDir = False

steps = ['hostAlignment', 'removeAligned', 'pathogenAlignment', 'sort', 'runQualimap', 
            'removeHighError', 'callReditools', 'callJacusa', 'filterSOR', 'runSnpEff', 
            'mergeCalling', 'generateGraphs']
if (
        config.cropLowQualityRuns.lower() == "on" and 
        config.runFastQC.lower() == "on" and 
        (step in steps or config.resumeFrom in steps)
):
    files = glob.glob(f"{config.workPath}/1-quality/srasCrop.pkl")
    if len(files) == 0:
        log.error("The step 'cropLowQualityRuns' has not been executed. Please, run it first, or turn it off in the main.config file.")
        util.stopProgram()
    sras = joblib.load(f"{config.workPath}/1-quality/srasCrop.pkl")

steps = ['pathogenAlignment', 'sort', 'runQualimap', 'removeHighError', 'callReditools', 
            'callJacusa', 'filterSOR', 'runSnpEff', 'mergeCalling', 'generateGraphs', 'generateReport'] 
if (
        config.removeAlignedWithHost.lower() == "on" and 
        config.runHostAlignment.lower() == "on" and 
        (step in steps or config.resumeFrom in steps)
):
    files = glob.glob(f"{config.workPath}/2-alignment/srasAlignment.pkl")
    if len(files) == 0:
        log.error("The step 'removeAligned' has not been executed. Please, run it first, or turn it off in the main.config file.")
        util.stopProgram()
    sras = joblib.load(f"{config.workPath}/2-alignment/srasAlignment.pkl")

steps = ['callReditools', 'callJacusa', 'filterSOR', 'runSnpEff', 'mergeCalling', 
        'generateGraphs', 'generateReport']
if (
    config.removeHighErrorRuns.lower() == "on" and 
    config.runQualimap.lower() == "on" and 
    (step in steps or config.resumeFrom in steps)
):
    files = glob.glob(f"{config.workPath}/3-qualimap/srasQualimap.pkl")
    if len(files) == 0:
        log.error("The step 'discardHighError' has not been executed. Please, run it first, or turn it off in the main.config file.")
        util.stopProgram()
    sras = joblib.load(f"{config.workPath}/3-qualimap/srasQualimap.pkl")

if len(sras) == 0:
    if config.source == "project":
        sras = downloads.getExperimentsList()
    elif config.source == "sra":
        with open("../sras.txt", "r") as f:
            # Skip the first line
            line = f.readline()
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip()
                sra = line.split()[0]
                sras.append(sra)

    if config.source != "file":
        sras = downloads.downloadSRAs(sras)
    else:
        useFastqInputDir = True
        fastqDir = config.inputFastqDir
        with open(f"../{config.inputType}Input.txt", "r") as f:
            # Skip the first line
            line = f.readline()
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip()
                splits = line.split()
                sra = splits[0]
                runId = splits[1]
                runType = splits[2]
                files = splits[3:]
                if runType == "single":
                    if not os.path.exists(f"{fastqDir}/{files[0]}"):
                        log.error(f"FASTQ file {files[0]} does not exist.")
                        util.stopProgram()
                    sras.append(([f"{fastqDir}/{files[0]}"], "single", runId))
                else:
                    if len(files) < 2:
                        log.error(f"{sra} is paired-end, but it does not have two files.")
                    if not os.path.exists(f"{fastqDir}/{files[0]}"):
                        log.error(f"FASTQ file {files[0]} does not exist.")
                        util.stopProgram()
                    if not os.path.exists(f"{fastqDir}/{files[1]}"):
                        log.error(f"FASTQ file {files[1]} does not exist.")
                        util.stopProgram()
                    sras.append(([f"{fastqDir}/{files[0]}", f"{fastqDir}/{files[1]}"], "paired", runId))

if config.resumeFrom in ["None", "runFastQC"] and config.runFastQC.lower() == "on" and step in ("", "runFastQC"):
    log.info("Running step runFastQC...")
    quality.runFastQC(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality"] and config.runFastQC.lower() == "on" and config.cropLowQualityRuns.lower() == "on" and step in ("", "cropLowQuality"):
    log.info("Running step cropLowQuality...")
    toCrop = cropping.detectFilesToCrop(sras)
    if len(toCrop) == 0:
        log.info("No files to be cropped.")
    if config.cropSoftware == "trimmomatic":
        cropping.runTrimmomatic(toCrop)
    elif config.cropSoftware == "trimgalore":
        cropping.runTrimGalore(toCrop)
    else:
        log.error("Invalid cropping software in configuration file. Accepted values are 'trimmomatic' and 'trimgalore'.")
        util.stopProgram()
    qualityDir = config.workPath + "/1-quality"
    fastqDir = qualityDir + "/fastq"
    for cropped in toCrop:
        for i in range(len(sras)):
            if all(f in sras[i][0] for f in cropped[0]):
                files = []
                if sras[i][1] == "single":
                    files.append(f"{fastqDir}/{sras[i][2]}.fastq")
                else:
                    files.append(f"{fastqDir}/{sras[i][2]}_1.fastq")
                    files.append(f"{fastqDir}/{sras[i][2]}_2.fastq")
                sras[i] = (files, sras[i][1], sras[i][2])
    joblib.dump(sras, f"{qualityDir}/srasCrop.pkl")

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment"] and config.runHostAlignment.lower() == "on" and step in ("", "hostAlignment"):
    log.info("Running step hostAlignment...")
    if config.originalHostPath == "":
        log.error("The path to the host reference FASTA file is required.")
        util.stopProgram()

    path = pathlib.Path(config.originalHostPath[0])
    util.makeDirectory(f"{config.workPath}/2-alignment/host/genomes/{path.stem}")
    util.execCmd(f"cp {path} {config.workPath}/2-alignment/host/genomes/{path.stem}/genome.fa")

    if config.alignmentSoftwareHost == "hisat2":
        alignment.runHisat2(sras, True)
    elif config.alignmentSoftwareHost == "bwa":
        alignment.runBWA(sras, True)
    elif config.alignmentSoftwareHost == "star":
        alignment.runSTAR(sras, True)
    elif config.alignmentSoftwareHost == "magicblast":
        alignment.runMagicBlast(sras, True)
    elif config.alignmentSoftwareHost == "minimap2":
        alignment.runMinimap2(sras, True)
    elif config.alignmentSoftwareHost == "gmap":
        alignment.runGMAP(sras, True)
    else:
        log.error("Invalid alignment software for the host in configuration file. Accepted values are 'hisat2', 'star', 'magicblast', 'minimap2' and 'gmap'.")
        util.stopProgram()

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned"] and config.runHostAlignment.lower() == "on" and config.removeAlignedWithHost.lower() == "on" and step in ("", "removeAligned"):
    log.info("Running step removeAligned...")
    sras = alignment.extractUnaligned(sras)
    if len(sras) == 0:
        log.error("All runs were discarded.")
        util.stopProgram()
    else:
        joblib.dump(sras, f"{config.workPath}/2-alignment/srasAlignment.pkl")

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment"] and step in ("", "pathogenAlignment"):
    log.info("Running step pathogenAlignment...")
    if config.pathogenReferenceGenomePaths == []:
        log.error("The paths to the pathogen reference FASTA file(s) are required.")
        util.stopProgram()

    if config.alignmentSoftwarePathogen == "hisat2":
        alignment.runHisat2(sras, False)
    elif config.alignmentSoftwarePathogen == "bwa":
        alignment.runBWA(sras, False)
    elif config.alignmentSoftwarePathogen == "star":
        alignment.runSTAR(sras, False)
    elif config.alignmentSoftwarePathogen == "magicblast":
        alignment.runMagicBlast(sras, False)
    elif config.alignmentSoftwarePathogen == "minimap2":
        alignment.runMinimap2(sras, False)
    elif config.alignmentSoftwarePathogen == "gmap":
        alignment.runGMAP(sras, False)
    else:
        log.error("Invalid alignment software for the pathogen in configuration file. Accepted values are 'hisat2', 'star', 'magicblast', 'minimap2' and 'gmap'.")
        util.stopProgram()

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort"] and step in ("", "sort"):
    log.info("Running step sort...")
    alignment.sortAlignments(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap"]  and config.runQualimap.lower() == "on" and step in ("", "runQualimap"):
    log.info("Running step runQualimap...")
    quality.runQualimap(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap", "discardHighError"] and config.runQualimap.lower() == "on" and config.removeHighErrorRuns.lower() == "on"  and step in ("", "discardHighError"):
    log.info("Running step discardHighError...")
    sras = quality.discardHighError(sras)
    if len(sras) == 0:
        log.error("All runs were discarded.")
        util.stopProgram()
    else:
        joblib.dump(sras, f"{config.workPath}/3-qualimap/srasQualimap.pkl")

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap", "discardHighError", "callReditools"] and step in ("", "callReditools") and config.callingSoftware in ["reditools", "both"]:
    log.info("Running step callReditools...")
    calling.runReditools2(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa"] and step in ("", "callJacusa") and config.callingSoftware in ["jacusa", "both"]:
    log.info("Running step callJacusa...")
    calling.runJacusa(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa", "filterSOR"] and step in ("", "filterSOR"):
    log.info("Running step filterSOR...")
    calling.filterAS_StrandOddsRatio(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa", "filterSOR", "runSnpEff"] and step in ("", "runSnpEff"):
    log.info("Running step runSnpEff...")
    results.runSnpEff(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa", "filterSOR", "runSnpEff", "mergeCalling"] and step in ("", "mergeCalling"):
    log.info("Running step mergeCalling...")
    results.mergeCalling(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa", "filterSOR", "runSnpEff", "mergeCalling", 'generateGraphs'] and step in ("", "generateGraphs"):
    log.info("Generating graphs...")
    results.generateGraphs(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "pathogenAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa", "filterSOR", "runSnpEff", "mergeCalling", 'generateGraphs', 'generateReport'] and step in ("", "generateReport"):
    log.info("Generating analysis report...")
    report.generateHtmlReport()

util.stopProgram()
