import sys
import glob
import software
import downloads
import quality
import cropping
import alignment
import sorting
import calling
import results
import util
import logger
import config
import re
import joblib
from argparse import ArgumentParser

log = logger.logger

parser = ArgumentParser()
parser.add_argument("-s", "--step", dest="step", nargs="?", type=str, required=False, help="only execute a specific step", metavar="step")
parser.add_argument("-w", "--workingDir", dest="workingDir" help="path to working directory", required=True)
parser.add_argument("-h", "--referenceHostPath", dest="referenceHostPath" help="path to reference host genome file", required=True)
parser.add_argument("-r", "--referenceVirusPath", dest="referenceVirusPath" help="path to directory with reference virus genome files", required=True)
parser.add_argument("-PENDING", "--referenceVirusPath", dest="referenceVirusPath" help="path to directory with reference virus genome files", required=True)
parser.add_argument("-d", "--download", help="download all required tools", action="store_true")

args = parser.parse_args()

if args.download == True:
    software.downloadTools()

steps = ["None", 'runFastQC', 'cropLowQuality', 'hostAlignment', 'removeAligned', 'virusAlignment', 'sort', 'runQualimap', 'discardHighError', 'callReditools', 'callJacusa', 'filterAS_StrandOddsRatio', 'mergeCalling']

step = ""
if args.step is not None:
    step = args.step
    if step not in steps:
        log.error("Invalid execution step. Accepted values are 'runFastQC', 'cropLowQuality', 'hostAlignment', 'removeAligned', 'virusAlignment', 'sort', 'runQualimap', 'removeHighError', 'callReditools', 'callJacusa', 'filterAS_StrandOddsRatio', 'mergeCalling'.")
        util.stopProgram()
    log.info(f"Executing step {step} only...")

if config.source not in ["project", "file", "sra"]:
    log.error("Invalid source in configuration file. Accepted values are 'project', 'file' and 'sra'.")
    util.stopProgram()

if config.resumeFrom not in ["None", 'runFastQC', 'cropLowQuality', 'hostAlignment', 'removeAligned', 'virusAlignment', 'sort', 'runQualimap', 'discardHighError', 'callReditools', 'callJacusa', 'filterAS_StrandOddsRatio', 'mergeCalling']:
    log.error("Invalid resuming step in configuration file. Accepted values are 'None', 'runFastQC', 'cropLowQuality', 'hostAlignment', 'removeAligned', 'virusAlignment', 'sort', 'runQualimap', 'removeHighError', 'callReditools', 'callJacusa', 'filterAS_StrandOddsRatio', 'mergeCalling'.")
    util.stopProgram()
elif config.resumeFrom != "None":
    log.info(f"Resuming from step {config.resumeFrom}...")

sras = []

if config.source == "project":
    sras = downloads.getExperimentsList()
elif config.source == "sra":
    with open(f"../config/sras.config", "r") as f:
        # Skip the first line
        line = f.readline()
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            sra = line.split()[0]
            group = line.split()[1]
            sras.append((sra, group))

if config.source != "file":
    sras = downloads.downloadSRAs(sras)
else:
    fastqDir = config.workPath + "/fastq"
    with open(f"../config/files.config", "r") as f:
        # Skip the first line
        line = f.readline()
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            sra = line.split()[0]
            type = line.split()[1]
            if type == "single":
                files = glob.glob(f"{fastqDir}/{sra}*.fastq")
                files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
                for file in files:
                    cmd = f"cat {file}"
                    util.execCmd(cmd, f"{fastqDir}/{sra}.fastq")
                sras.append(([f"{sra}.fastq"], "single"))
            else:
                files = glob.glob(f"{fastqDir}/{sra}*.fastq")
                files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
                for i in range(len(files)):
                    files[i] = files[i].split("/")[-1]
                sras.append((files, "paired"))

if config.resumeFrom in ["None", "runFastQC"] and config.runFastQC.lower() == "on" and step in ("", "runFastQC"):
    quality.runFastQC(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality"] and config.runFastQC.lower() == "on" and config.cropLowQualityRuns.lower() == "on" and step in ("", "cropLowQuality"):
    toCrop = cropping.detectFilesToCrop(sras)
    if config.cropSoftware == "trimmomatic":
        cropping.runTrimmomatic(toCrop)
    elif config.cropSoftware == "trimgalore":
        cropping.runTrimGalore(toCrop)
    else:
        log.error("Invalid cropping software in configuration file. Accepted values are 'trimmomatic' and 'trimgalore'.")
        util.stopProgram()

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment"] and config.runHostAlignment.lower() == "on" and step in ("", "hostAlignment"):
    if config.hostReferencePath == "":
        log.error("The path to the host reference FASTA file is required.")
        util.stopProgram()
        
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
    sras = alignment.extractUnaligned(sras)
    if len(sras) == 0:
        log.error("All runs were discarded.")
        util.stopProgram()
    else:
        joblib.dump(sras, f"{config.workPath}/srasAlignment.pkl")

if config.runHostAlignment.lower() == "on" and config.removeAlignedWithHost.lower() == "on" and step in ("", "virusAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa", "filterAS_StrandOddsRatio", "mergeCalling"):
    files = glob.glob(f"{config.workPath}/srasAlignment.pkl")
    if len(files) == 0:
        log.error("The step 'removeAlignedWithHost' has not been executed. Please, run it first, or turn it off in the main.config file.")
        util.stopProgram()
    sras = joblib.load(f"{config.workPath}/srasAlignment.pkl")

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "virusAlignment"] and step in ("", "virusAlignment"):
    if config.virusReferencePath == "":
        log.error("The path to the virus reference FASTA file(s) is required.")
        util.stopProgram()

    if config.alignmentSoftwareVirus == "hisat2":
        alignment.runHisat2(sras, False)
    elif config.alignmentSoftwareVirus == "bwa":
        alignment.runBWA(sras, False)
    elif config.alignmentSoftwareVirus == "star":
        alignment.runSTAR(sras, False)
    elif config.alignmentSoftwareVirus == "magicblast":
        alignment.runMagicBlast(sras, False)
    elif config.alignmentSoftwareVirus == "minimap2":
        alignment.runMinimap2(sras, False)
    elif config.alignmentSoftwareVirus == "gmap":
        alignment.runGMAP(sras, False)
    else:
        log.error("Invalid alignment software for the virus in configuration file. Accepted values are 'hisat2', 'star', 'magicblast', 'minimap2' and 'gmap'.")
        util.stopProgram()

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "virusAlignment", "sort"] and step in ("", "sort"):
    sorting.runSamtools(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "virusAlignment", "sort", "runQualimap"]  and config.runQualimap.lower() == "on" and step in ("", "runQualimap"):
    quality.runQualimap(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "virusAlignment", "sort", "runQualimap", "discardHighError"] and config.runQualimap.lower() == "on" and config.removeHighErrorRuns.lower() == "on"  and step in ("", "discardHighError"):
    sras = quality.discardHighError(sras)
    if len(sras) == 0:
        log.error("All runs were discarded.")
        util.stopProgram()
    else:
        joblib.dump(sras, f"{config.workPath}/srasQualimap.pkl")

if config.runQualimap.lower() == "on" and config.removeHighErrorRuns.lower() == "on" and step in ("", "callReditools", "callJacusa", "filterAS_StrandOddsRatio", "mergeCalling"):
    files = glob.glob(f"{config.workPath}/srasQualimap.pkl")
    if len(files) == 0:
        log.error("The step 'discardHighError' has not been executed. Please, run it first, or turn it off in the main.config file.")
        util.stopProgram()
    sras = joblib.load(f"{config.workPath}/srasQualimap.pkl")
    
if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "virusAlignment", "sort", "runQualimap", "discardHighError", "callReditools"] and step in ("", "callReditools"):
    calling.runReditools2(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "virusAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa"] and step in ("", "callJacusa"):
    calling.runJacusa(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "virusAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa", "filterAS_StrandOddsRatio"] and step in ("", "filterAS_StrandOddsRatio"):
    calling.filterAS_StrandOddsRatio(sras)

if config.resumeFrom in ["None", "runFastQC", "cropLowQuality", "hostAlignment", "removeAligned", "virusAlignment", "sort", "runQualimap", "discardHighError", "callReditools", "callJacusa", "filterAS_StrandOddsRatio", "mergeCalling"] and step in ("", "mergeCalling"):
    results.mergeCalling(sras)

util.stopProgram()
