import util
import config
import logger
import glob

log = logger.logger


def runSamtools(sras):
    samDir = config.workPath + "/sam"
    bamDir = config.workPath + "/bam"
    util.makeDirectory(bamDir)
    for fullRef in config.virusReferencePath:
        ref = fullRef[fullRef.rfind("/") + 1:]
        ref = ref[:ref.rfind(".")]        
        util.makeDirectory(f"{bamDir}/{ref}")
        for sra in sras:
            run = sra[0][0]
            run = run[:run.rfind(".")]
            type = sra[1]
            if type == "paired":
                run = run.replace("_1", "")
            files = glob.glob(f"{config.workPath}/sam/{ref}/{run}_virus_{config.alignmentSoftwareVirus}.sam")
            if len(files) == 0:
                log.error(f"Alignment for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            log.info(f"Sorting {run} vs {ref} BAM file...")
            cmd = f"{config.samtoolsPath} sort -O BAM -o {bamDir}/{ref}/{run}_virus_{config.alignmentSoftwareVirus}.bam {samDir}/{ref}/{run}_virus_{config.alignmentSoftwareVirus}.sam"
            util.execCmd(cmd)