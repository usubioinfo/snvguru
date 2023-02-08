import util
import config
import logger
import glob

log = logger.logger

def detectFilesToCrop(sras):
    fastqcDir = config.workPath + "/fastqc"
    toCrop = []
    for sra in sras:
        type = sra[1]
        for file in sra[0]:
            added = False
            dir = file.replace(".fastq", "_fastqc")
            files = glob.glob(f"{fastqcDir}/{dir}/fastqc_data.txt")
            if len(files) == 0:
                log.error(f"FastQC analysis for file {file} not found.")
                util.stopProgram()
            with open(f"{fastqcDir}/{dir}/fastqc_data.txt") as f:
                values = []
                for _ in range(13):
                    f.readline()
                for line in f:
                    if line.startswith(">>END_MODULE"):
                        break
                    line = line.strip()
                    val = float(line.split()[1])
                    values.append(val)
                maxValue = -1
                for val in values:
                    if val < config.cropMinMeanQuality:
                        toCrop.append(sra)
                        added = True
                        break
                    elif val > maxValue:
                        maxValue = val
                    elif maxValue - config.cropMaxDecay > val:
                        toCrop.append(sra)
                        added = True
                        break
                if added:
                    break
            
    return toCrop

def runTrimmomatic(files):
    for file in files:
        fastqDir = config.workPath + "/fastq"
        filename = file[0]
        type = file[1]
        filesString = " and ".join(filename)
        log.info(f"Cropping {filesString}...")
        if len(glob.glob(f"{fastqDir}/{filename[0]}.uncropped")) > 0 and len(glob.glob(f"{fastqDir}/{filename[0]}")) > 0:
            log.info(f"{filesString} already cropped. Skipping...")
        else:
            cmd = f"mv {fastqDir}/{filename[0]} {fastqDir}/{filename[0]}.uncropped"
            util.execCmd(cmd)
            if type == "single":
                cmd = f"java -jar {config.trimmomaticPath} SE {fastqDir}/{filename[0]}.uncropped {fastqDir}/{filename[0]} CROP:{config.cropSize}"
            else:
                cmd = f"mv {fastqDir}/{filename[1]} {fastqDir}/{filename[1]}.uncropped"
                util.execCmd(cmd)
                cmd = f"java -jar {config.trimmomaticPath} PE {fastqDir}/{filename[0]}.uncropped {fastqDir}/{filename[1]}.uncropped {fastqDir}/{filename[0]} {fastqDir}/{filename[0]}.cropped.unpaired {fastqDir}/{filename[1]} {fastqDir}/{filename[1]}.cropped.unpaired CROP:{config.cropSize}"
            util.execCmd(cmd)

def runTrimGalore(files):
    for file in files:
        fastqDir = config.workPath + "/fastq"
        filename = file[0]
        type = file[1]
        filesString = " and ".join(filename)
        log.info(f"Cropping {filesString}...")
        if len(glob.glob(f"{fastqDir}/{filename[0]}.uncropped")) > 0 and len(glob.glob(f"{fastqDir}/{filename[0]}")) > 0:
            log.info(f"{filesString} already cropped. Skipping...")
        else:
            cmd = f"mv {fastqDir}/{filename[0]} {fastqDir}/{filename[0]}.uncropped"
            util.execCmd(cmd)
            if type == "single":
                cmd = f"{config.trimGalorePath} -o {fastqDir} --hardtrim5 {config.cropSize} {fastqDir}/{filename[0]}.uncropped"
                util.execCmd(cmd)
                cmd = f"mv {fastqDir}/{filename[0]}.uncropped.{config.cropSize}bp_5prime.fq {fastqDir}/{filename[0]}"
                util.execCmd(cmd)
            else:
                cmd = f"mv {fastqDir}/{filename[1]} {fastqDir}/{filename[1]}.uncropped"
                util.execCmd(cmd)
                cmd = f"{config.trimGalorePath} -o {fastqDir} --paired --hardtrim5 {config.cropSize} {fastqDir}/{filename[0]}.uncropped {fastqDir}/{filename[1]}.uncropped"
                util.execCmd(cmd)
                cmd = f"mv {fastqDir}/{filename[0]}.uncropped.{config.cropSize}bp_5prime.fq {fastqDir}/{filename[0]}"
                util.execCmd(cmd)
                cmd = f"mv {fastqDir}/{filename[1]}.uncropped.{config.cropSize}bp_5prime.fq {fastqDir}/{filename[1]}"
                util.execCmd(cmd)
