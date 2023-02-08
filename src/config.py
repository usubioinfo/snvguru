import os
import glob
import re

os.chdir(os.path.dirname(__file__))

hisatH = ""
hisatV = ""
bwaIndexH = ""
bwaMappingH = ""
bwaIndexV = ""
bwaMappingV = ""
starIndexH = ""
starMappingH = ""
starIndexV = ""
starMappingV = ""
magicblastH = ""
magicblastV = ""
minimapIndexH = ""
minimapPresetH = ""
minimapMappingH = ""
minimapIndexV = ""
minimapPresetV = ""
minimapMappingV = ""
gmapIndexH = ""
gmapMappingH = ""
gmapIndexV = ""
gmapMappingV = ""
reditools = ""
jacusa = ""
bcftools = ""
qualimap = ""

with open("../config/main.config", "r") as f:
    for line in f:
        line = line.strip()
        # Source
        if line.startswith("source"):
            source = line.split()[1]
        # Threads
        elif line.startswith("threads"):
            val = line.split()[1]
            if val == "None":
                threads = ""
            else:
                threads = f"{int(val)}"
        # Work path
        elif line.startswith("workPath"):
            workPath = line.split()[1]
            if not os.path.exists(workPath):
                os.makedirs(workPath)
        # Log
        elif line.startswith("logFile"):
            logFile = True if line.split()[1] == "True" else False
        elif line.startswith("logConsole"):
            logConsole = True if line.split()[1] == "True" else False
        # Optional steps
        elif line.startswith("runFastQC"):
            runFastQC = line.split()[1]
        elif line.startswith("cropLowQualityRuns"):
            cropLowQualityRuns = line.split()[1]
        elif line.startswith("runHostAlignment"):
            runHostAlignment = line.split()[1]
        elif line.startswith("removeAlignedWithHost"):
            removeAlignedWithHost = line.split()[1]
        elif line.startswith("runQualimap"):
            runQualimap = line.split()[1]
        elif line.startswith("removeHighErrorRuns"):
            removeHighErrorRuns = line.split()[1]
        # Resume from a specific step
        elif line.startswith("resumeFrom"):
            resumeFrom = line.split()[1]
        # SRA toolkit path
        elif line.startswith("sratoolkitPath"):
            val = line.split()[1]
            if val == "None":
                sratoolkitPath = ""
            else:
                sratoolkitPath = f"{val}/"
        # SAMtools path
        elif line.startswith("samtoolsPath"):
            val = line.split()[1]
            if val == "None":
                samtoolsPath = "samtools"
            else:
                samtoolsPath = val
        # FastQC
        elif line.startswith("fastqcPath"):
            val = line.split()[1]
            if val == "None":
                fastqcPath = "fastqc"
            else:
                fastqcPath = val
        elif line.startswith("fastqcContaminantsFile"):
            val = line.split()[1]
            if val == "None":
                fastqcContaminantsFile = ""
            else:
                fastqcContaminantsFile = f" -c {val}"
        elif line.startswith("fastqcAdapterFile"):
            val = line.split()[1]
            if val == "None":
                fastqcAdapterFile = ""
            else:
                fastqcAdapterFile = f" -a {val}"
        elif line.startswith("fastqcMinLength"):
            val = line.split()[1]
            if val == "None":
                fastqcMinLength = ""
            else:
                fastqcMinLength = f" --min_length {int(val)}"
        elif line.startswith("fastqcQuiet"):
            val = line.split()[1]
            if val == "None" or val == "False":
                fastqcQuiet = ""
            else:
                fastqcQuiet = " -q"
        # Qualimap
        elif line.startswith("qualimapPath"):
            val = line.split()[1]
            if val == "None":
                qualimapPath = "qualimap"
            else:
                qualimapPath = val
        elif line.startswith("qualimapSeqProtocol"):
            val = line.split()[1]
            if val == "None":
                qualimapSeqProtocol = ""
            else:
                qualimapSeqProtocol = f" -p {val}"
        elif line.startswith("qualimapMinHomopol"):
            val = line.split()[1]
            if val == "None":
                qualimap += ""
            else:
                qualimap += f" -hm {val}"
        elif line.startswith("qualimapReadsPerChunk"):
            val = line.split()[1]
            if val == "None":
                qualimap += ""
            else:
                qualimap += f" -nr {val}"
        elif line.startswith("qualimapWindows"):
            val = line.split()[1]
            if val == "None":
                qualimap += ""
            else:
                qualimap += f" -nw {val}"
        elif line.startswith("qualimapMaxError"):
            val = line.split()[1]
            qualimapMaxError = float(val) / 100
        
        # Cropping
        elif line.startswith("cropSoftware"):
            val = line.split()[1]
            if val == "None":
                cropSoftware = ""
            else:
                cropSoftware = val
        elif line.startswith("trimmomaticPath"):
            val = line.split()[1]
            if val == "None":
                trimmomaticPath = "trimmomatic.jar"
            else:
                trimmomaticPath = val
        elif line.startswith("trimGalorePath"):
            val = line.split()[1]
            if val == "None":
                trimGalorePath = "trim_galore"
            else:
                trimGalorePath = val
        elif line.startswith("cropMinMeanQuality"):
            val = line.split()[1]
            if val == "None":
                cropMinMeanQuality = -1
            else:
                cropMinMeanQuality = float(val)
        elif line.startswith("cropMaxDecay"):
            val = line.split()[1]
            if val == "None":
                cropMaxDecay = 999999
            else:
                cropMaxDecay = float(val)
        elif line.startswith("cropSize"):
            val = line.split()[1]
            if val == "None":
                cropSize = ""
            else:
                cropSize = int(val)
        # Alignment
        elif line.startswith("alignmentSoftware"):
            vals = line.split()
            if vals[1] == "None":
                alignmentSoftwareHost = ""
            else:
                alignmentSoftwareHost = vals[1]
            if vals[2] == "None":
                alignmentSoftwareVirus = ""
            else:
                alignmentSoftwareVirus = vals[2]
        elif line.startswith("hostReferencePath"):
            val = line.split()[1]
            if val == "None":
                hostReferencePath = ""
            else:
                hostReferencePath = [val]
        elif line.startswith("virusReferencePath"):
            val = line.split()[1]
            if val == "None":
                virusReferencePath = [""]
            elif val.endswith("/"):
                files = glob.glob(f"{val}*.fasta")
                files.sort(key=lambda x: [int(c) if c.isdigit(
                ) else c for c in re.split(r'(\d+)', x)])
                virusReferencePath = files
            else:
                virusReferencePath = [val]
        # Calling - REDitools2
        elif line.startswith("reditoolsCommand"):
            val = line.split()[1]
            if val == "None":
                reditoolsPath = "python2 reditools.py"
            else:
                reditoolsPath = val
        elif line.startswith("reditoolsStrand"):
            val = line.split()[1]
            if val == "None":
                reditoolsStrand = 0
            else:
                reditoolsStrand = int(val)
        elif line.startswith("reditoolsHomopolSpan"):
            val = line.split()[1]
            if val == "None":
                reditoolsHomopolSpan = 4
            else:
                reditoolsHomopolSpan = int(val)
        elif line.startswith("reditoolsMinBasePos"):
            val = line.split()[1]
            if val == "None":
                reditools += ""
            else:
                reditools += f" -mbp {int(val)}"
        elif line.startswith("reditoolsMaxBasePos"):
            val = line.split()[1]
            if val == "None":
                reditools += ""
            else:
                reditools += f" -Mbp {int(val)}"
        elif line.startswith("reditoolsMinReadLen"):
            val = line.split()[1]
            if val == "None":
                reditools += ""
            else:
                reditools += f" -mrl {int(val)}"
        # Calling - JACUSA
        elif line.startswith("jacusaPath"):
            val = line.split()[1]
            if val == "None":
                jacusaPath = "jacusa.jar"
            else:
                jacusaPath = val
        elif line.startswith("jacusaFeatureFilter"):
            val = line.split()[1]
            if val == "None":
                jacusaFeatureFilter = "B,I,Y"
            else:
                jacusaFeatureFilter = val
        elif line.startswith("jacusa_flag_"):
            split = line.split()
            param = split[0].replace("jacusa_flag_", "")
            param = f" -{param} "
            if split[1] == "true":
                jacusa += param
        elif line.startswith("jacusa_"):
            split = line.split()
            param = split[0].replace("jacusa_", "")
            param = f" -{param} "
            jacusa += param + split[1]
        # CALLING - bcftools
        elif line.startswith("bcftoolsPath"):
            val = line.split()[1]
            if val == "None":
                bcftoolsPath = "bcftools"
            else:
                bcftoolsPath = val
        elif line.startswith("bcftools_flag_"):
            split = line.split()
            param = split[0].replace("bcftools_flag_", "")
            param = f" -{param} "
            if split[1] == "true":
                bcftools += param
        elif line.startswith("bcftools_"):
            split = line.split()
            param = split[0].replace("bcftools_", "")
            param = f" -{param} "
            bcftools += param + split[1]
        # Calling - General
        elif line.startswith("callingReadMinQuality"):
            val = line.split()[1]
            reditools += f" -q {int(val)}"
            jacusa += f" -m {int(val)}"
            bcftools += f" -q {int(val)}"
        elif line.startswith("callingBaseMinQuality"):
            val = line.split()[1]
            reditools += f" -bq {int(val)}"
            jacusa += f" -q {int(val)}"
            bcftools += f" -Q {int(val)}"
        # Results
        elif line.startswith("minSNVCoverage"):
            val = line.split()[1]
            minSNVCoverage = int(val)
        elif line.startswith("minMainReadSupport"):
            val = line.split()[1]
            minMainReadSupport = int(val)
        elif line.startswith("minRecurringReadSupport"):
            val = line.split()[1]
            minRecurringReadSupport = int(val)
        elif line.startswith("minFrequency"):
            val = line.split()[1]
            minFrequency = float(val) / 100
        elif line.startswith("maxAS_StrandOddsRatio"):
            val = line.split()[1]
            maxAS_StrandOddsRatio = float(val)
        # Figures
        elif line.startswith("figDistributionWidth"):
            val = line.split()[1]
            figDistributionWidth = float(val)
        elif line.startswith("figDistributionTicksX"):
            val = line.split()[1]
            figDistributionTicksX = int(val)
        elif line.startswith("figDistributionBinSize"):
            val = line.split()[1]
            figDistributionBinSize = int(val)
        elif line.startswith("figDPI"):
            val = line.split()[1]
            figDPI = int(val)
# Hisat2
with open("../config/alignHisat2.config", "r") as f:
    for line in f:
        if line.startswith("hisat2Path"):
            val = line.split()[1]
            if val == "None":
                hisat2Path = ""
            else:
                hisat2Path = f"{val}/"
        elif line.startswith("hisat2Mismatch"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                hisatH += ""
            else:
                hisatH += f" --mp {val1}"
            if val2 == "None":
                hisatV += ""
            else:
                hisatV += f" --mp {val2}"
        elif line.startswith("hisat2SoftClipping"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                hisatH += ""
            else:
                hisatH += f" --sp {val1}"
            if val2 == "None":
                hisatV += ""
            else:
                hisatV += f" --sp {val2}"
        elif line.startswith("hisat2NoSoftClipping"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1.lower() == "true":
                hisatH += " --no-softclip"
            else:
                hisatH += ""
            if val2.lower() == "true":
                hisatV += " --no-softclip"
            else:
                hisatV += ""
        elif line.startswith("hisat2NonACGT"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                hisatH += ""
            else:
                hisatH += f" --np {val1}"
            if val2 == "None":
                hisatV += ""
            else:
                hisatV += f" --np {val2}"
        elif line.startswith("hisat2ReadGapOpenExtend"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                hisatH += ""
            else:
                hisatH += f" --rdg {val1}"
            if val2 == "None":
                hisatV += ""
            else:
                hisatV += f" --rdg {val2}"
        elif line.startswith("hisat2RefGapOpenExtend"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                hisatH += ""
            else:
                hisatH += f" --rfg {val1}"
            if val2 == "None":
                hisatV += ""
            else:
                hisatV += f" --rfg {val2}"
        elif line.startswith("hisat2MinScore"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                hisatH += ""
            else:
                hisatH += f" --score-min {val1}"
            if val2 == "None":
                hisatV += ""
            else:
                hisatV += f" --score-min {val2}"

# BWA
with open("../config/alignBWA.config", "r") as f:
    for line in f:
        if line.startswith("bwaPath"):
            val = line.split()[1]
            if val == "None":
                bwaPath = "bwa"
            else:
                bwaPath = val
        elif line.startswith("bwaIndexAlgorithm"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaIndexH += ""
            else:
                bwaIndexH += f" -a {val1}"
            if val2 == "None":
                bwaIndexV += ""
            else:
                bwaIndexV += f" -a {val2}"
        elif line.startswith("bwaMinSeedLength"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -k {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -k {int(val2)}"
        elif line.startswith("bwaBandwidth"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -w {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -w {int(val2)}"
        elif line.startswith("bwaOffDiagonalXDropoff"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -d {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -d {int(val2)}"
        elif line.startswith("bwaReseeding"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -r {float(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -r {float(val2)}"
        elif line.startswith("bwaMinOccurence"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -c {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -c {int(val2)}"
        elif line.startswith("bwaMatch"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -A {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -A {int(val2)}"
        elif line.startswith("bwaMismatch"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -B {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -B {int(val2)}"
        elif line.startswith("bwaGapOpen"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -O {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -O {int(val2)}"
        elif line.startswith("bwaGapExt"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -E {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -E {int(val2)}"
        elif line.startswith("bwaClipping"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -L {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -L {int(val2)}"
        elif line.startswith("bwaUnpaired"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -U {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -U {int(val2)}"
        elif line.startswith("bwaMinScore"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                bwaMappingH += ""
            else:
                bwaMappingH += f" -T {int(val1)}"
            if val2 == "None":
                bwaMappingV += ""
            else:
                bwaMappingV += f" -T {int(val2)}"
        elif line.startswith("bwaOutputSinglePaired"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "true":
                bwaMappingH += " -a"
            if val2 == "true":
                bwaMappingV += " -a"

# STAR
with open("../config/alignSTAR.config", "r") as f:
    for line in f:
        if line.startswith("starPath"):
            val = line.split()[1]
            if val == "None":
                starPath = "STAR"
            else:
                starPath = val
        elif line.startswith("starGTF"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starIndexH += ""
            else:
                starIndexH += f" --sjdbGTFfile {val1}"
            if val2 == "None":
                starIndexV += ""
            else:
                starIndexV += f" --sjdbGTFfile {val2}"
        elif line.startswith("starOverhang"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starIndexH += ""
            else:
                starIndexH += f" --sjdbOverhang {int(val1)}"
            if val2 == "None":
                starIndexV += ""
            else:
                starIndexV += f" --sjdbOverhang {int(val2)}"
        elif line.startswith("starGenomeSAindexNBases"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --genomeSAindexNbases {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --genomeSAindexNbases {int(val2)}"
        elif line.startswith("starGap") and line.split()[0] == "starGap":
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --scoreGap {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --scoreGap {int(val2)}"
        elif line.startswith("starGapNonCanonical"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --scoreGapNoncan {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --scoreGapNoncan {int(val2)}"
        elif line.startswith("starGapGCAG"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --scoreGapGCAG {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --scoreGapGCAG {int(val2)}"
        elif line.startswith("starGapATAC"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --scoreGapATAC {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --scoreGapATAC {int(val2)}"
        elif line.startswith("starDeleteOpen"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --scoreDelOpen {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --scoreDelOpen {int(val2)}"
        elif line.startswith("starDeleteExtend"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --scoreDelBase {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --scoreDelBase {int(val2)}"
        elif line.startswith("starInsertOpen"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --scoreInsOpen {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --scoreInsOpen {int(val2)}"
        elif line.startswith("starInsertExtend"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                starMappingH += ""
            else:
                starMappingH += f" --scoreInsBase {int(val1)}"
            if val2 == "None":
                starMappingV += ""
            else:
                starMappingV += f" --scoreInsBase {int(val2)}"
        elif line.startswith("starIndex_flag_"):
            params = line.split()
            param = line.split()[0].replace("starIndex_flag_", "")
            if params[1] == "true":
                starIndexH += f" --{param}"
            if params[2] == "true":
                starIndexV += f" --{param}"
        elif line.startswith("starIndex_"):
            params = line.split()
            param = line.split()[0].replace("starIndex_", "")
            starIndexH += f" --{param} {params[1]}"
            starIndexV += f" --{param} {params[2]}"
        elif line.startswith("starMapping_flag_"):
            params = line.split()
            param = line.split()[0].replace("starMapping_flag_", "")
            if params[1].lower() == "true":
                starMappingH += f" --{param}"
            if params[2].lower() == "true":
                starMappingV += f" --{param}"
        elif line.startswith("starMapping_"):
            params = line.split()
            param = line.split()[0].replace("starMapping_", "")
            starMappingH += f" --{param} {params[1]}"
            starMappingV += f" --{param} {params[2]}"

# Magic-BLAST
with open("../config/alignMagicBLAST.config", "r") as f:
    for line in f:
        if line.startswith("magicblastPath"):
            val = line.split()[1]
            if val == "None":
                magicblastPath = ""
            else:
                magicblastPath = f"{val}/"
        elif line.startswith("magicblastWordSize"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -word_size {int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -word_size {int(val2)}"
        elif line.startswith("magicblastGapOpen"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -gapopen {int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -gapopen {int(val2)}"
        elif line.startswith("magicblastGapExtend"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -gapextend {int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -gapextend {int(val2)}"
        elif line.startswith("magicBlastMismatch"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -penalty -{int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -penalty -{int(val2)}"
        elif line.startswith("magicblastMaxIntronLen"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -max_intron_length {int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -max_intron_length {int(val2)}"
        elif line.startswith("magicblastPercIdentity"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -perc_identity {float(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -perc_identity {float(val2)}"
        elif line.startswith("magicblastLcaseMask"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None" or val1.lower() == "false":
                magicblastH += ""
            else:
                magicblastH += " -lcase_masking"
            if val2 == "None" or val2.lower() == "false":
                magicblastV += ""
            else:
                magicblastV += " -lcase_masking"
        elif line.startswith("magicblastTranscriptome"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "false":
                magicblastH += ""
            else:
                magicblastH += f" -reftype transcriptome"
            if val2 == "false":
                magicblastV += ""
            else:
                magicblastV += f" -reftype transcriptome"
        elif line.startswith("magicblastMaxDbWordCount"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -max_db_word_count {int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -max_db_word_count {int(val2)}"
        elif line.startswith("magicblastLimitLookup"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1.lower() == "false":
                magicblastH += " -limit_lookup false"
            else:
                magicblastH += ""
            if val2.lower() == "false":
                magicblastV += " -limit_lookup false"
            else:
                magicblastV += ""
        elif line.startswith("magicblastLookupStride"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -lookup_stride {int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -lookup_stride {int(val2)}"
        elif line.startswith("magicblastValidateSeqs"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1.lower() == "false":
                magicblastH += " -validate_seqs false"
            else:
                magicblastH += ""
            if val2.lower() == "false":
                magicblastV += " -validate_seqs false"
            else:
                magicblastV += ""
        elif line.startswith("magicblastScore"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -score {int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -score {int(val2)}"
        elif line.startswith("magicblastMaxEditDist"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                magicblastH += ""
            else:
                magicblastH += f" -max_edit_dist {int(val1)}"
            if val2 == "None":
                magicblastV += ""
            else:
                magicblastV += f" -max_edit_dist {int(val2)}"
        elif line.startswith("magicblastSplice"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1.lower() == "false":
                magicblastH += " -splice false"
            else:
                magicblastH += ""
            if val2.lower() == "false":
                magicblastV += " -splice false"
            else:
                magicblastV += ""

# Minimap2
with open("../config/alignMinimap2.config", "r") as f:
    for line in f:
        if line.startswith("minimapPath"):
            val = line.split()[1]
            if val == "None":
                minimapPath = "minimap2"
            else:
                minimapPath = val
        elif line.startswith("minimapKmerSize"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapIndexH += ""
            else:
                minimapIndexH += f" -k {int(val1)}"
            if val1 == "None":
                minimapIndexV += ""
            else:
                minimapIndexV += f" -k {int(val2)}"
        elif line.startswith("minimapWindowSize"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapIndexH += ""
            else:
                minimapIndexH += f" -w {int(val1)}"
            if val2 == "None":
                minimapIndexV += ""
            else:
                minimapIndexV += f" -w {int(val2)}"
        elif line.startswith("minimapSplitIndex"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapIndexH += ""
            else:
                minimapIndexH += f" -I {val1}"
            if val2 == "None":
                minimapIndexV += ""
            else:
                minimapIndexV += f" -I {val2}"
        elif line.startswith("minimapPreset"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapPresetH += "splice"
            else:
                minimapPresetH += val1
            if val2 == "None":
                minimapPresetV += "splice"
            else:
                minimapPresetV += val2
        elif line.startswith("minimapGTAG"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -u {val1}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -u {val2}"
        elif line.startswith("minimapFraction"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -f {val1}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -f {val2}"
        elif line.startswith("minimapStopChainElongation"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -g {val1}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -g {val2}"
        elif line.startswith("minimapMaxIntron"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -G {val1}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -G {val2}"
        elif line.startswith("minimapChainAlign"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -r {val1}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -r {val2}"
        elif line.startswith("minimapMinMinimizers"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -n {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -n {int(val2)}"
        elif line.startswith("minimapChainingScore"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -m {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -m {int(val2)}"
        elif line.startswith("minimapSecondary"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" --secondary={val1}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" --secondary={val2}"
        elif line.startswith("minimapMin2to1ScoreRatio"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -p {val1}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -p {val2}"
        elif line.startswith("minimapRetainSecondary"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -N {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -N {int(val2)}"
        elif line.startswith("minimapMatch"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -A {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -A {int(val2)}"
        elif line.startswith("minimapMismatch"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -B {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -B {int(val2)}"
        elif line.startswith("minimapOpen"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -O {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -O {int(val2)}"
        elif line.startswith("minimapExtension"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -E {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -E {int(val2)}"
        elif line.startswith("minimapZDrop"):
            val1 = line.split()[1]
            vals1 = val1.split(",")
            val2 = line.split()[2]
            vals2 = val2.split(",")
            if val1 == "None":
                minimapMappingH += ""
            else:
                if len(vals1) == 1:
                    minimapMappingH += f" -z {int(vals1[0])}"
                else:
                    minimapMappingH += f" -z {int(vals1[0])},{int(vals1[1])}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                if len(vals2) == 1:
                    minimapMappingV += f" -z {int(vals2[0])}"
                else:
                    minimapMappingV += f" -z {int(vals2[0])},{int(vals2[1])}"
        elif line.startswith("minimapMinPeakDP"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -s {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -s {int(val2)}"
        elif line.startswith("minimapMinibatch"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                minimapMappingH += ""
            else:
                minimapMappingH += f" -K {int(val1)}"
            if val2 == "None":
                minimapMappingV += ""
            else:
                minimapMappingV += f" -K {int(val2)}"
        elif line.startswith("minimapIndex_flag_"):
            split = line.split()
            param = split[0].replace("minimapIndex_flag_", "")
            if len(param) == 1:
                param = f" -{param} "
            else:
                param = f" --{param} "
            if split[1] == "true":
                minimapIndexH += param
            if split[2] == "true":
                minimapIndexV += param
        elif line.startswith("minimapIndex_"):
            split = line.split()
            param = split[0].replace("minimapIndex_", "")
            if len(param) == 1:
                param = f" -{param} "
            else:
                param = f" --{param} "
            minimapIndexH += param + split[1]
            minimapIndexV += param + split[2]
        elif line.startswith("minimapMapping_flag_"):
            split = line.split()
            param = split[0].replace("minimapMapping_flag_", "")
            if len(param) == 1:
                param = f" -{param} "
            else:
                param = f" --{param} "
            if split[1] == "true":
                minimapMappingH += param
            if split[2] == "true":
                minimapMappingV += param
        elif line.startswith("minimapMapping_"):
            split = line.split()
            param = split[0].replace("minimapMapping_", "")
            if len(param) == 1:
                param = f" -{param} "
            else:
                param = f" --{param} "
            minimapMappingH += param + split[1]
            minimapMappingV += param + split[2]

# GMAP
with open("../config/alignGMAP.config", "r") as f:
    for line in f:
        if line.startswith("gmapPath"):
            val = line.split()[1]
            if val == "None":
                gmapPath = ""
            else:
                gmapPath = f"{val}/"
        elif line.startswith("gmapIndexKmer"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapIndexH += ""
            else:
                gmapIndexH += f" -k {int(val1)}"
            if val2 == "None":
                gmapIndexV += ""
            else:
                gmapIndexV += f" -k {int(val2)}"
        elif line.startswith("gmapInterval"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapIndexH += ""
            else:
                gmapIndexH += f" -q {int(val1)}"
            if val2 == "None":
                gmapIndexV += ""
            else:
                gmapIndexV += f" -q {int(val2)}"
        elif line.startswith("gmapMapKmer"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" -k {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" -k {int(val2)}"
        elif line.startswith("gmapSampling"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --sampling {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --sampling {int(val2)}"
        elif line.startswith("gmapPart"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" -q {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" -q {int(val2)}"
        elif line.startswith("gmapBuffer"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --input-buffer-size {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --input-buffer-size {int(val2)}"
        elif line.startswith("gmapSplicing"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "false":
                gmapMappingH += " --nosplicing"
            else:
                gmapMappingH += ""
            if val2 == "false":
                gmapMappingV += " --nosplicing"
            else:
                gmapMappingV += ""
        elif line.startswith("gmapMinIntron"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --min-intronlength {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --min-intronlength {int(val2)}"
        elif line.startswith("gmapMaxIntronMiddle"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --max-intronlength-middle {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --max-intronlength-middle {int(val2)}"
        elif line.startswith("gmapMaxIntronEnds"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --max-intronlength-ends {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --max-intronlength-ends {int(val2)}"
        elif line.startswith("gmapSplitLargeIntrons"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "true":
                gmapMappingH += " --split-large-introns"
            else:
                gmapMappingH += f""
            if val2 == "true":
                gmapMappingV += " --split-large-introns"
            else:
                gmapMappingV += f""
        elif line.startswith("gmapTrimEndExons"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --trim-end-exons {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --trim-end-exons {int(val2)}"
        elif line.startswith("gmapMaxKnownSplice"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" -w {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" -w {int(val2)}"
        elif line.startswith("gmapMaxIntron"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" -L {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" -L {int(val2)}"
        elif line.startswith("gmapChimeraMargin"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" -x {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" -x {int(val2)}"
        elif line.startswith("gmapDirection"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" -z {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" -z {int(val1)}"
        elif line.startswith("gmapCanonical"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --canonical-mode {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --canonical-mode {int(val2)}"
        elif line.startswith("gmapCrossSpecies"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "true":
                gmapMappingH += "--cross-species"
            else:
                gmapMappingH += ""
            if val2 == "true":
                gmapMappingV += "--cross-species"
            else:
                gmapMappingV += ""
        elif line.startswith("gmapAllowCloseIndels"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --allow-close-indels {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --allow-close-indels {int(val2)}"
        elif line.startswith("gmapMicroexonSpliceProb"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" --microexon-spliceprob {float(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" --microexon-spliceprob {float(val2)}"
        elif line.startswith("gmapPruning"):
            val1 = line.split()[1]
            val2 = line.split()[2]
            if val1 == "None":
                gmapMappingH += ""
            else:
                gmapMappingH += f" -p {int(val1)}"
            if val2 == "None":
                gmapMappingV += ""
            else:
                gmapMappingV += f" -p {int(val2)}"
        elif line.startswith("gmapIndex_flag_"):
            split = line.split()
            param = split[0].replace("gmapIndex_flag_", "")
            if len(param) == 1:
                param = f" -{param}"
            else:
                param = f" --{param}"
            if split[1] == "true":
                gmapIndexH += param
            if split[2] == "true":
                gmapIndexV += param
        elif line.startswith("gmapIndex_"):
            split = line.split()
            param = split[0].replace("gmapIndex_", "")
            if len(param) == 1:
                param = f" -{param}"
            else:
                param = f" --{param}"
            gmapIndexH += param + split[1]
            gmapIndexV += param + split[2]
        elif line.startswith("gmapMapping_flag_"):
            split = line.split()
            param = split[0].replace("gmapMapping_flag_", "")
            if len(param) == 1:
                param = f" -{param}"
            else:
                param = f" --{param}"
            gmapMappingH += param + split[1]
            gmapMappingV += param + split[2]
        elif line.startswith("gmapMapping_"):
            split = line.split()
            param = split[0].replace("gmapMapping_", "")
            if len(param) == 1:
                param = f" -{param}"
            else:
                param = f" --{param}"
            if split[1] == "true":
                gmapMappingH += param
            if split[2] == "true":
                gmapMappingV += param
