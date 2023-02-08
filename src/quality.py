import util
import config
import logger
import glob

log = logger.logger


def runFastQC(sras):
    fastqDir = config.workPath + "/fastq"
    fastqcDir = config.workPath + "/fastqc"
    util.makeDirectory(fastqcDir)
    fastqcDirTemp = config.workPath + "/fastqctmp"
    util.makeDirectory(fastqcDirTemp)
    for sra in sras:
        files = []
        files = f"{fastqDir}/" + f" {fastqDir}/".join(sra[0])
        filesLog = " ".join(sra[0])
        log.info(f"Analyzing the quality of {filesLog}...")
        cmd = f"{config.fastqcPath} --extract -t {config.threads} -o {fastqcDir}/ {files}{config.fastqcContaminantsFile}{config.fastqcAdapterFile}{config.fastqcMinLength}{config.fastqcQuiet}"
        util.execCmd(cmd)


def runQualimap(sras):
    qualDir = config.workPath + "/qualimap"
    util.makeDirectory(qualDir)
    bamDir = config.workPath + "/bam"
    for fullRef in config.virusReferencePath:
        ref = fullRef[fullRef.rfind("/") + 1:]
        ref = ref[:ref.rfind(".")]
        util.makeDirectory(f"{qualDir}/{ref}")
        for sra in sras:
            run = sra[0][0]
            run = run[:run.rfind(".")]
            type = sra[1]
            if type == "paired":
                run = run.replace("_1", "")
            util.makeDirectory(f"{qualDir}/{ref}/{run}")
            files = glob.glob(f"{bamDir}/{ref}/{run}_virus_{config.alignmentSoftwareVirus}.bam")
            if len(files) == 0:
                log.error(f"Sorted BAM file with the alignment for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            log.info(f"Analyzing the alignment of {run} vs {ref}...")
            cmd = f"{config.qualimapPath} bamqc -bam {bamDir}/{ref}/{run}_virus_{config.alignmentSoftwareVirus}.bam -outdir {qualDir}/{ref}/{run} -nt {config.threads}{config.qualimap}"
            util.execCmd(cmd)

def discardHighError(sras):
    log.info(f"Discarding high error runs...")
    qualDir = config.workPath + "/qualimap"
    runs = []
    newSras = []
    for fullRef in config.virusReferencePath:
        ref = fullRef[fullRef.rfind("/") + 1:]
        ref = ref[:ref.rfind(".")]
        for sra in sras:
            run = sra[0][0]
            run = run[:run.rfind(".")]
            type = sra[1]
            if type == "paired":
                run = run.replace("_1", "")
            files = glob.glob(f"{qualDir}/{ref}/{run}/genome_results.txt")
            if len(files) == 0:
                log.error(f"Qualimap results for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            with open(f"{qualDir}/{ref}/{run}/genome_results.txt") as f:
                for line in f:
                    if line.lstrip().startswith("general error rate"):
                        value = float(line.strip().split("=")[1].strip())
                        if value < config.qualimapMaxError:
                            if run not in runs:
                                runs.append(run)
                                newSras.append(sra)
    for sra in sras:
        run = sra[0][0]
        run = run[:run.rfind(".")]
        type = sra[1]
        if type == "paired":
            run = run.replace("_1", "")
        if run not in runs:
            log.info(f"Discarded run {run}.")
    return newSras
