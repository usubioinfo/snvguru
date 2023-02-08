import util
import config
import logger

log = logger.logger

def filterCommon(sras):
    callingDir = config.workPath + "/calling"
    for fullRef in config.virusReferencePath:
        ref = fullRef[fullRef.rfind("/") + 1:]
        ref = ref[:ref.rfind(".")]
        dir = f"{callingDir}/{ref}" 
        for sra in sras:
            run = sra[0][0]
            run = run[:run.rfind(".")]
            type = sra[1]
            if type == "paired":
                run = run.replace("_1", "")

            with open(f"{run}.calling.txt") as f:
                first = True
                for line in f:
                    if first:
                        first = False
                        continue
                    data = line.split()
                    pos = data[1]
                    refBase = data[2]
                    coverage = data[3]
                    qual = data[4]
                    baseCt = data[5][1:-1].split(",")
                    cts = {"A": baseCt[0], "C": baseCt[1], "G": baseCt[2], "T": baseCt[3]}
                    freq = data[7]


def calculateAllelicDepths(sras):
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
        file = f"{run}_virus_{config.alignmentSoftware}"
        cmd = f"{config.bcftoolsPath} mpileup -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP -O v -A -C -I -d 1000000{config.bcftools} -o {depthsDir}/{run}.vcf {bamDir}/{file}.bam"
        util.execCmd(cmd)
