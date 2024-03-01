"""This module reads the configuration files and sets up the global 
configuration variables.
"""

import os
import glob
import re
import shlex

os.chdir(os.path.dirname(__file__))

step = ""
hisatIndexH = hisatIndexV = hisatMappingH = hisatMappingV = ""
bwaIndexH = bwaMappingH = bwaIndexV = bwaMappingV = ""
starIndexH = starMappingH = starIndexV = starMappingV = ""
magicblastIndexH = magicblastIndexV = magicblastMappingH = magicblastMappingV = ""
minimapIndexH = minimapPresetH = minimapMappingH = minimapIndexV = minimapPresetV = minimapMappingV = ""
gmapIndexH = gmapMappingH = gmapIndexV = gmapMappingV = ""
reditools = jacusa = bcftools = qualimap = fastqc = ""
pathogenReferenceGenomePaths = []
pathogenReferenceGenesPaths = []
pathogenReferenceProteinPaths = []
pathogenReferenceGenesFormats = []
workPath = ""

def setWorkPath(workPathP):
    global workPath 
    workPath = workPathP

def loadConfig(configDir="../config/", stepP=None, resumeFromP=None):
    global reditools, jacusa, bcftools, workPath
    with open(f"{configDir}/main.config", "r") as f:
        for line in f:
            line = line.strip()
            # Source
            if line.startswith("source"):
                global source
                source = line.split()[1]
            # Threads
            elif line.startswith("threads"):
                global threads
                val = line.split()[1]
                if val == "None":
                    threads = ""
                else:
                    threads = f"{int(val)}"
            # SLURM
            # Run on SLURM
            elif line.startswith("slurmTime"):
                global slurmTime
                val = line.split()[1]
                slurmTime = f"{int(val)}"
            elif line.startswith("slurmMem"):
                global slurmMem
                val = line.split()[1]
                slurmMem = val
            elif line.startswith("slurmCpus"):
                global slurmCpus
                val = line.split()[1]
                slurmCpus = f"{int(val)}"
            elif line.startswith("slurm"):
                global slurm
                val = line.split()[1]
                if val.lower() == "true":
                    slurm = True
                else:
                    slurm = False
            # Work path
            elif line.startswith("workPath") and workPath == "":
                workPath = line.split()[1]
                if not os.path.exists(workPath):
                    os.makedirs(workPath)
                elif not stepP and not resumeFromP:
                    i = 1
                    while os.path.exists(f'{workPath}_{i}'):
                        i += 1
                    os.mkdir(f'{workPath}_{i}')
                    workPath = f'{workPath}_{i}'
            # Tools path
            elif line.startswith("toolsPath"):
                global toolsPath
                toolsPath = line.split()[1]
                if not os.path.exists(toolsPath):
                    os.makedirs(toolsPath)
            # Input FASTQ type -- Single-end, paired-end or mixed
            elif line.startswith("inputType"):
                global inputType
                inputType = line.split()[1]
            # Input FASTQ directory -- ONLY FOR SOURCE = 'file'
            elif line.startswith("inputFastqDir"):
                global inputFastqDir
                inputFastqDir = line.split()[1]
            # Log
            elif line.startswith("logFile"):
                global logFile
                logFile = True if line.split()[1] == "True" else False
            elif line.startswith("logConsole"):
                global logConsole
                logConsole = True if line.split()[1] == "True" else False
            # Optional steps
            elif line.startswith("runFastQC"):
                global runFastQC
                runFastQC = line.split()[1]
            elif line.startswith("cropLowQualityRuns"):
                global cropLowQualityRuns
                cropLowQualityRuns = line.split()[1]
            elif line.startswith("runHostAlignment"):
                global runHostAlignment
                runHostAlignment = line.split()[1]
            elif line.startswith("removeAlignedWithHost"):
                global removeAlignedWithHost
                removeAlignedWithHost = line.split()[1]
            elif line.startswith("runQualimap"):
                global runQualimap
                runQualimap = line.split()[1]
            elif line.startswith("removeHighErrorRuns"):
                global removeHighErrorRuns
                removeHighErrorRuns = line.split()[1]
            # Resume from a specific step
            elif line.startswith("resumeFrom"):
                global resumeFrom
                resumeFrom = line.split()[1]
                if resumeFromP:
                    resumeFrom = resumeFromP
            # Java path
            elif line.startswith("javaPath"):
                global javaPath
                val = line.split()[1]
                if val == "None":
                    javaPath = f"{toolsPath}/jdk-13.0.2/bin/java"
                else:
                    javaPath = val
            # SRA toolkit path
            elif line.startswith("sratoolkitPath"):
                global sratoolkitPath
                val = line.split()[1]
                if val == "None":
                    sratoolkitPath = ""
                else:
                    sratoolkitPath = f"{val}/"
            # SAMtools path
            elif line.startswith("samtoolsPath"):
                global samtoolsPath
                val = line.split()[1]
                if val == "None":
                    samtoolsPath = "samtools"
                else:
                    samtoolsPath = val
            # FastQC
            elif line.startswith("fastqcPath"):
                global fastqcPath
                val = line.split()[1]
                if val == "None":
                    fastqcPath = f"{toolsPath}/FastQC/fastqc"
                else:
                    fastqcPath = val
            # Qualimap
            elif line.startswith("qualimapPath"):
                global qualimapPath
                val = line.split()[1]
                if val == "None":
                    qualimapPath = f"{toolsPath}/qualimap_v2.2.1/qualimap"
                else:
                    qualimapPath = val
            elif line.startswith("qualimapMaxError"):
                global qualimapMaxError
                val = line.split()[1]
                qualimapMaxError = float(val) / 100
            # Cropping
            elif line.startswith("cropSoftware"):
                global cropSoftware
                val = line.split()[1]
                if val == "None":
                    cropSoftware = ""
                else:
                    cropSoftware = val
            elif line.startswith("trimmomaticPath"):
                global trimmomaticPath
                val = line.split()[1]
                if val == "None":
                    trimmomaticPath = f"{toolsPath}/Trimmomatic-0.39/trimmomatic-0.39.jar"
                else:
                    trimmomaticPath = val
            elif line.startswith("trimGalorePath"):
                global trimGalorePath
                val = line.split()[1]
                if val == "None":
                    trimGalorePath = "trim_galore"
                else:
                    trimGalorePath = val
            elif line.startswith("cropMinMeanQuality"):
                global cropMinMeanQuality
                val = line.split()[1]
                if val == "None":
                    cropMinMeanQuality = -1
                else:
                    cropMinMeanQuality = float(val)
            elif line.startswith("cropMaxDecay"):
                global cropMaxDecay
                val = line.split()[1]
                if val == "None":
                    cropMaxDecay = 999999
                else:
                    cropMaxDecay = float(val)
            elif line.startswith("cropSize"):
                global cropSize
                val = line.split()[1]
                if val == "None":
                    cropSize = ""
                else:
                    cropSize = int(val)
            # Alignment
            elif line.startswith("alignmentSoftwareHost"):
                global alignmentSoftwareHost
                vals = line.split()
                if vals[1] == "None":
                    alignmentSoftwareHost = ""
                else:
                    alignmentSoftwareHost = vals[1]
            elif line.startswith("alignmentSoftwarePathogen"):
                global alignmentSoftwarePathogen
                vals = line.split()
                if vals[1] == "None":
                    alignmentSoftwarePathogen = ""
                else:
                    alignmentSoftwarePathogen = vals[1]
            elif line.startswith("hostReferencePath"):
                global originalHostPath
                val = line.split()[1]
                if val == "None":
                    originalHostPath = ""
                else:
                    originalHostPath = [val]
            elif line.startswith("originalPathogenGenomePaths"):
                global originalPathogenGenomePaths
                val = line.split()[1]
                if val == "None":
                    originalPathogenGenomePaths = [""]
                else:
                    originalPathogenGenomePaths = val
            elif line.startswith("originalPathogenGenesPaths"):
                global originalPathogenGenesPaths
                val = line.split()[1]
                if val == "None":
                    originalPathogenGenesPaths = [""]
                else:
                    originalPathogenGenesPaths = val
            elif line.startswith("originalPathogenProteinPaths"):
                global originalPathogenProteinPaths
                val = line.split()[1]
                if val == "None":
                    originalPathogenProteinPaths = [""]
                else:
                    originalPathogenProteinPaths = val
            # Calling - General
            elif line.startswith("bcftoolsPath"):
                global bcftoolsPath
                val = line.split()[1]
                if val == "None":
                    bcftoolsPath = "bcftools"
                else:
                    bcftoolsPath = val
            elif line.startswith("jacusaPath"):
                global jacusaPath
                val = line.split()[1]
                if val == "None":
                    jacusaPath = f"{toolsPath}/jacusa/JACUSA_v2.0.2-RC.jar"
                else:
                    jacusaPath = val
            elif line.startswith("reditoolsCommand"):
                global reditoolsCommand
                val = shlex.split(line)[1]
                if val == "None":
                    reditoolsCommand = f"{toolsPath}/REDItools2/env/bin/python2 {toolsPath}/REDItools2/src/cineca/reditools.py"
                else:
                    reditoolsCommand = val
            elif line.startswith("callingReadMinQuality"):
                global callingReadMinQuality
                val = line.split()[1]
                callingReadMinQuality = int(val)
                reditools += f" -q {int(val)}"
                jacusa += f" -m {int(val)}"
                bcftools += f" -q {int(val)}"
            elif line.startswith("callingBaseMinQuality"):
                global callingBaseMinQuality
                val = line.split()[1]
                callingBaseMinQuality = int(val)
                reditools += f" -bq {int(val)}"
                jacusa += f" -q {int(val)}"
                bcftools += f" -Q {int(val)}"
            elif line.startswith("callingSoftware"):
                global callingSoftware
                val = line.split()[1]
                if val == "None":
                    callingSoftware = ""
                else:
                    callingSoftware = val
            # Results
            elif line.startswith("snpEffPath"):
                global snpEffPath
                val = line.split()[1]
                if val == "None":
                    snpEffPath = f"{toolsPath}/snpEff/snpEff.jar"
                else:
                    snpEffPath = val
            elif line.startswith("minSNVCoverage"):
                global minSNVCoverage
                val = line.split()[1]
                minSNVCoverage = int(val)
            elif line.startswith("minMainReadSupport"):
                global minMainReadSupport
                val = line.split()[1]
                minMainReadSupport = int(val)
            elif line.startswith("minRecurringReadSupport"):
                global minRecurringReadSupport
                val = line.split()[1]
                minRecurringReadSupport = int(val)
            elif line.startswith("minFrequency"):
                global minFrequency
                val = line.split()[1]
                minFrequency = float(val) / 100.0
            elif line.startswith("maxAS_StrandOddsRatio"):
                global maxAS_StrandOddsRatio
                val = line.split()[1]
                maxAS_StrandOddsRatio = float(val)
            # Figures
            elif line.startswith("figDPI"):
                global figDPI
                val = line.split()[1]
                figDPI = int(val)
            elif line.startswith("figPositionGraphMinIntensity"):
                global figPositionGraphMinIntensity
                val = line.split()[1]
                figPositionGraphMinIntensity = float(val)
            elif line.startswith("figPositionGraph"):
                global figPositionGraph
                val = line.split()[1]
                figPositionGraph = val
            elif line.startswith("figGenerateMutationCountBarPlotsPerRun"):
                global figGenerateMutationCountBarPlotsPerRun
                val = line.split()[1]
                figGenerateMutationCountBarPlotsPerRun = val
            elif line.startswith("figGenerateMutationCountPerGeneBarPlotPerRun"):
                global figGenerateMutationCountPerGeneBarPlotPerRun
                val = line.split()[1]
                figGenerateMutationCountPerGeneBarPlotPerRun = val
            elif line.startswith("figGenerateMutationCountBoxPlotsPerRun"):
                global figGenerateMutationCountBoxPlotsPerRun
                val = line.split()[1]
                figGenerateMutationCountBoxPlotsPerRun = val
            elif line.startswith("figGenerateMutationCountPerGeneBoxPlotPerRun"):
                global figGenerateMutationCountPerGeneBoxPlotPerRun
                val = line.split()[1]
                figGenerateMutationCountPerGeneBoxPlotPerRun = val                
            elif line.startswith("figGenerateFrequencyPerMutationPerPositionGraphsPerRun"):
                global figGenerateFrequencyPerMutationPerPositionGraphsPerRun
                val = line.split()[1]
                figGenerateFrequencyPerMutationPerPositionGraphsPerRun = val                
            elif line.startswith("figGenerateGlobalMutationCountBarPlots"):
                global figGenerateGlobalMutationCountBarPlots
                val = line.split()[1]
                figGenerateGlobalMutationCountBarPlots = val                
            elif line.startswith("figGenerateGlobalMutationCountPerGeneBarPlots"):
                global figGenerateGlobalMutationCountPerGeneBarPlots
                val = line.split()[1]
                figGenerateGlobalMutationCountPerGeneBarPlots = val         
            elif line.startswith("figGenerateGlobalMutationCountPerRunBarPlots"):
                global figGenerateGlobalMutationCountPerRunBarPlots
                val = line.split()[1]
                figGenerateGlobalMutationCountPerRunBarPlots = val                    
            elif line.startswith("figGenerateGlobalMutationCountBoxPlots"):
                global figGenerateGlobalMutationCountBoxPlots
                val = line.split()[1]
                figGenerateGlobalMutationCountBoxPlots = val                
            elif line.startswith("figGenerateGlobalMutationCountPerGeneBoxPlots"):
                global figGenerateGlobalMutationCountPerGeneBoxPlots
                val = line.split()[1]
                figGenerateGlobalMutationCountPerGeneBoxPlots = val          
            elif line.startswith("figGenerateGlobalMutationCountPerRunBoxPlots"):
                global figGenerateGlobalMutationCountPerRunBoxPlots
                val = line.split()[1]
                figGenerateGlobalMutationCountPerRunBoxPlots = val              
            elif line.startswith("figGenerateGlobalFrequencyPerMutationStripPlots"):
                global figGenerateGlobalFrequencyPerMutationStripPlots
                val = line.split()[1]
                figGenerateGlobalFrequencyPerMutationStripPlots = val                
            elif line.startswith("figGenerateGlobalFrequencyPerGeneStripPlots"):
                global figGenerateGlobalFrequencyPerGeneStripPlots
                val = line.split()[1]
                figGenerateGlobalFrequencyPerGeneStripPlots = val                
            elif line.startswith("figGenerateGlobalFrequencyPerRunStripPlots"):
                global figGenerateGlobalFrequencyPerRunStripPlots
                val = line.split()[1]
                figGenerateGlobalFrequencyPerRunStripPlots = val                
            elif line.startswith("figGenerateGlobalDistributionHistogramsPlots"):
                global figGenerateGlobalDistributionHistogramsPlots
                val = line.split()[1]
                figGenerateGlobalDistributionHistogramsPlots = val                
            elif line.startswith("figGenerateGlobalRegressionPlots"):
                global figGenerateGlobalRegressionPlots
                val = line.split()[1]
                figGenerateGlobalRegressionPlots = val                
            elif line.startswith("figGenerateGlobalPresencePerRunPerPositionGraphs"):
                global figGenerateGlobalPresencePerRunPerPositionGraphs
                val = line.split()[1]
                figGenerateGlobalPresencePerRunPerPositionGraphs = val                
            # Bar plots
            elif line.startswith("figBarPlotPerGeneWidth"):
                global figBarPlotPerGeneWidth
                val = line.split()[1]
                figBarPlotPerGeneWidth = int(val)
            elif line.startswith("figBarPlotPerGeneHeight"):
                global figBarPlotPerGeneHeight
                val = line.split()[1]
                figBarPlotPerGeneHeight = int(val)
            elif line.startswith("figBarPlotPerGeneTickSize"):
                global figBarPlotPerGeneTickSize
                val = line.split()[1]
                figBarPlotPerGeneTickSize = int(val)
            elif line.startswith("figBarPlotPerGeneAxisLabelSize"):
                global figBarPlotPerGeneAxisLabelSize
                val = line.split()[1]
                figBarPlotPerGeneAxisLabelSize = int(val)
            elif line.startswith("figBarPlotPerGeneTitleSize"):
                global figBarPlotPerGeneTitleSize
                val = line.split()[1]
                figBarPlotPerGeneTitleSize = int(val)
            elif line.startswith("figBarPlotPerGeneNumber"):
                global figBarPlotPerGeneNumber
                val = line.split()[1]
                figBarPlotPerGeneNumber = int(val)
            elif line.startswith("figBarPlotPerRunWidth"):
                global figBarPlotPerRunWidth
                val = line.split()[1]
                figBarPlotPerRunWidth = int(val)
            elif line.startswith("figBarPlotPerRunBarHeight"):
                global figBarPlotPerRunBarHeight
                val = line.split()[1]
                figBarPlotPerRunBarHeight = float(val)
            elif line.startswith("figBarPlotCountWidth"):
                global figBarPlotCountWidth
                val = line.split()[1]
                figBarPlotCountWidth = int(val)
            elif line.startswith("figBarPlotCountHeight"):
                global figBarPlotCountHeight
                val = line.split()[1]
                figBarPlotCountHeight = int(val)
            # Box plots
            elif line.startswith("figBoxPlotPerGeneWidth"):
                global figBoxPlotPerGeneWidth
                val = line.split()[1]
                figBoxPlotPerGeneWidth = int(val)
            elif line.startswith("figBoxPlotPerGeneHeight"):
                global figBoxPlotPerGeneHeight
                val = line.split()[1]
                figBoxPlotPerGeneHeight = int(val)
            elif line.startswith("figBoxPlotPerGeneTickSize"):
                global figBoxPlotPerGeneTickSize
                val = line.split()[1]
                figBoxPlotPerGeneTickSize = int(val)
            elif line.startswith("figBoxPlotPerGeneAxisLabelSize"):
                global figBoxPlotPerGeneAxisLabelSize
                val = line.split()[1]
                figBoxPlotPerGeneAxisLabelSize = int(val)
            elif line.startswith("figBoxPlotPerGeneTitleSize"):
                global figBoxPlotPerGeneTitleSize
                val = line.split()[1]
                figBoxPlotPerGeneTitleSize = int(val)
            elif line.startswith("figBoxPlotPerGeneNumber"):
                global figBoxPlotPerGeneNumber
                val = line.split()[1]
                figBoxPlotPerGeneNumber = int(val)
            elif line.startswith("figBoxPlotPerRunWidth"):
                global figBoxPlotPerRunWidth
                val = line.split()[1]
                figBoxPlotPerRunWidth = int(val)
            elif line.startswith("figBoxPlotPerRunBarHeight"):
                global figBoxPlotPerRunBarHeight
                val = line.split()[1]
                figBoxPlotPerRunBarHeight = float(val)
            elif line.startswith("figBoxPlotCountWidth"):
                global figBoxPlotCountWidth
                val = line.split()[1]
                figBoxPlotCountWidth = int(val)
            elif line.startswith("figBoxPlotCountHeight"):
                global figBoxPlotCountHeight
                val = line.split()[1]
                figBoxPlotCountHeight = int(val)
            # Strip plots
            elif line.startswith("figStripPlotDotSize"):
                global figStripPlotDotSize
                val = line.split()[1]
                figStripPlotDotSize = int(val)
            elif line.startswith("figStripPlotPerGeneNumber"):
                global figStripPlotPerGeneNumber
                val = line.split()[1]
                figStripPlotPerGeneNumber = int(val)
            elif line.startswith("figStripPlotPerGeneTitleSize"):
                global figStripPlotPerGeneTitleSize
                val = line.split()[1]
                figStripPlotPerGeneTitleSize = int(val)
            elif line.startswith("figStripPlotPerMutationWidth"):
                global figStripPlotPerMutationWidth
                val = line.split()[1]
                figStripPlotPerMutationWidth = int(val)
            elif line.startswith("figStripPlotPerMutationHeight"):
                global figStripPlotPerMutationHeight
                val = line.split()[1]
                figStripPlotPerMutationHeight = int(val)
            elif line.startswith("figStripPlotPerRunWidth"):
                global figStripPlotPerRunWidth
                val = line.split()[1]
                figStripPlotPerRunWidth = int(val)
            elif line.startswith("figStripPlotPerRunBarHeight"):
                global figStripPlotPerRunBarHeight
                val = line.split()[1]
                figStripPlotPerRunBarHeight = float(val)
            elif line.startswith("figStripPlotPerGeneWidth"):
                global figStripPlotPerGeneWidth
                val = line.split()[1]
                figStripPlotPerGeneWidth = int(val)
            elif line.startswith("figStripPlotPerGeneHeight"):
                global figStripPlotPerGeneHeight
                val = line.split()[1]
                figStripPlotPerGeneHeight = int(val)
            # Histograms
            elif line.startswith("figDistributionTicksY"):
                global figDistributionTicksY
                val = line.split()[1]
                figDistributionTicksY = int(val)
            elif line.startswith("figDistributionBinSize"):
                global figDistributionBinSize
                val = line.split()[1]
                figDistributionBinSize = int(val)
            elif line.startswith("figDistributionWidth"):
                global figDistributionWidth
                val = line.split()[1]
                figDistributionWidth = int(val)
            elif line.startswith("figDistributionHeight"):
                global figDistributionHeight
                val = line.split()[1]
                figDistributionHeight = int(val)
            # Regression                
            elif line.startswith("figRegressionHeight"):
                global figRegressionHeight
                val = line.split()[1]
                figRegressionHeight = int(val)
            elif line.startswith("figRegressionWidth"):
                global figRegressionWidth
                val = line.split()[1]
                figRegressionWidth = int(val)
            # Circos                
            elif line.startswith("figPositionsPerCircos"):
                global figPositionsPerCircos
                val = line.split()[1]
                figPositionsPerCircos = int(val)
            elif line.startswith("figCircosSize"):
                global figCircosSize
                val = line.split()[1]
                figCircosSize = int(val)
            elif line.startswith("figCircosCenterSize"):
                global figCircosCenterSize
                val = line.split()[1]
                figCircosCenterSize = int(val)
            elif line.startswith("figCircosTitleSize"):
                global figCircosTitleSize
                val = line.split()[1]
                figCircosTitleSize = int(val)
            elif line.startswith("figCircosSampleLabelSize"):
                global figCircosSampleLabelSize
                val = line.split()[1]
                figCircosSampleLabelSize = int(val)
            elif line.startswith("figCircosMutationLabelSize"):
                global figCircosMutationLabelSize
                val = line.split()[1]
                figCircosMutationLabelSize = int(val)
            elif line.startswith("figCircosPositionLabelSize"):
                global figCircosPositionLabelSize
                val = line.split()[1]
                figCircosPositionLabelSize = int(val)
            elif line.startswith("figCircosColorBarLabelSize"):
                global figCircosColorBarLabelSize
                val = line.split()[1]
                figCircosColorBarLabelSize = int(val)
            elif line.startswith("figCircosColorBarTickSize"):
                global figCircosColorBarTickSize
                val = line.split()[1]
                figCircosColorBarTickSize = int(val)
            elif line.startswith("figCircosMutationsBlankDegrees"):
                global figCircosMutationsBlankDegrees
                val = line.split()[1]
                figCircosMutationsBlankDegrees = int(val)
            elif line.startswith("figCircosSamplesBlankDegrees"):
                global figCircosSamplesBlankDegrees
                val = line.split()[1]
                figCircosSamplesBlankDegrees = int(val)
            elif line.startswith("figCircosGeneLabelSize"):
                global figCircosGeneLabelSize
                val = line.split()[1]
                figCircosGeneLabelSize = int(val)
            elif line.startswith("figCircosSamplesColor"):
                global figCircosSamplesColor
                val = line.split()[1]
                figCircosSamplesColor = val
            elif line.startswith("figCircosMutationsColor"):
                global figCircosMutationsColor
                val = line.split()[1]
                figCircosMutationsColor = val
            # Heatmap
            elif line.startswith("figPositionsPerHeatmap"):
                global figPositionsPerHeatmap
                val = line.split()[1]
                figPositionsPerHeatmap = int(val)
            elif line.startswith("figHeatmapSize"):
                global figHeatmapSize
                val = line.split()[1]
                figHeatmapSize = int(val)
            elif line.startswith("figHeatmapTitlePadding"):
                global figHeatmapTitlePadding
                val = line.split()[1]
                figHeatmapTitlePadding = int(val)
            elif line.startswith("figHeatmapTitleSize"):
                global figHeatmapTitleSize
                val = line.split()[1]
                figHeatmapTitleSize = int(val)
            elif line.startswith("figHeatmapColorBarTitlePadding"):
                global figHeatmapColorBarTitlePadding
                val = line.split()[1]
                figHeatmapColorBarTitlePadding = int(val)
            elif line.startswith("figHeatmapColorBarLabelSize"):
                global figHeatmapColorBarLabelSize
                val = line.split()[1]
                figHeatmapColorBarLabelSize = int(val)
            elif line.startswith("figHeatmapTickSize"):
                global figHeatmapTickSize
                val = line.split()[1]
                figHeatmapTickSize = int(val)
            elif line.startswith("figHeatmapGeneLabelSize"):
                global figHeatmapGeneLabelSize
                val = line.split()[1]
                figHeatmapGeneLabelSize = int(val)
            elif line.startswith("figHeatmapSamplesColor"):
                global figHeatmapSamplesColor
                val = line.split()[1]
                figHeatmapSamplesColor = val
            elif line.startswith("figHeatmapMutationsColor"):
                global figHeatmapMutationsColor
                val = line.split()[1]
                figHeatmapMutationsColor = val
                
            # Aligner paths
            if line.startswith("bwaPath"):
                global bwaPath
                val = line.split()[1]
                if val == "None":
                    bwaPath = "bwa"
                else:
                    bwaPath = val
            if line.startswith("hisat2Path"):
                global hisat2Path
                val = line.split()[1]
                if val == "None":
                    hisat2Path = f"{toolsPath}/hisat2-2.1.0/"
                else:
                    hisat2Path = f"{val}/"
            if line.startswith("starPath"):
                global starPath
                val = line.split()[1]
                if val == "None":
                    starPath = "STAR"
                else:
                    starPath = val
            if line.startswith("magicblastPath"):
                global magicblastPath
                val = line.split()[1]
                if val == "None":
                    magicblastPath = f"{toolsPath}/ncbi-magicblast-1.6.0/bin/"
                else:
                    magicblastPath = f"{val}/"
            if line.startswith("minimapPath"):
                global minimapPath
                val = line.split()[1]
                if val == "None":
                    minimapPath = f"{toolsPath}/minimap2-2.20_x64-linux/minimap2"
                else:
                    minimapPath = val
            if line.startswith("gmapPath"):
                global gmapPath
                val = line.split()[1]
                if val == "None":
                    gmapPath = ""
                else:
                    gmapPath = f"{val}/"

    # READ QUALITY - FastQC
    with open(f"{configDir}fastqc.config", "r") as f:
        global fastqc
        line = f.readline()
        while line:
            if not line.startswith("#"):
                fastqc += line.split("#")[0].strip() + " "
            line = f.readline()

    # ALIGNMENT - Hisat2
    with open(f"{configDir}hisat2.config", "r") as f:
        global hisatIndexH, hisatMappingH, hisatIndexV, hisatMappingV
        line = f.readline()
        while line:
            if line.startswith("#"):
                line = f.readline()
            elif line.startswith("hostIndex:"):
                line = f.readline()
                while line.startswith("\t"):
                    if not line.startswith("#"):
                        hisatIndexH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("hostAlignment:"):
                line = f.readline()
                while line.startswith("\t"):
                    splits = line.split()
                    if not splits[0].startswith("#"):
                        hisatMappingH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenIndex:"):
                line = f.readline()
                while line.startswith("\t"):
                    if not line.startswith("#"):
                        hisatIndexV += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenAlignment:"):
                line = f.readline()
                while line.startswith("\t"):
                    splits = line.split()
                    if not line.startswith("#"):
                        hisatMappingV += line.split("#")[0].strip() + " "
                    line = f.readline()
            else:
                line = f.readline()

    # ALIGNMENT - BWA
    with open(f"{configDir}bwa.config", "r") as f:
        global bwaIndexH, bwaMappingH, bwaIndexV, bwaMappingV
        line = f.readline()
        while line:
            if line.startswith("#"):
                line = f.readline()
            elif line.startswith("hostIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        bwaIndexH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("hostAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if not splits[0].startswith("#"):
                        bwaMappingH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        bwaIndexV += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if not line.startswith("#"):
                        bwaMappingV += line.split("#")[0].strip() + " "
                    line = f.readline()
            else:
                line = f.readline()

    # ALIGNMENT - STAR
    with open(f"{configDir}star.config", "r") as f:
        global starIndexH, starMappingH, starIndexV, starMappingV
        line = f.readline()
        while line:
            if line.startswith("#"):
                line = f.readline()
            elif line.startswith("hostIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        starIndexH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("hostAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if not splits[0].startswith("#"):
                        starMappingH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        starIndexV += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if not line.startswith("#"):
                        starMappingV += line.split("#")[0].strip() + " "
                    line = f.readline()
            else:
                line = f.readline()

    # ALIGNMENT - Magic-BLAST
    with open(f"{configDir}magicblast.config", "r") as f:
        global magicblastIndexH, magicblastMappingH, magicblastIndexV, magicblastMappingV
        line = f.readline()
        while line:
            if line.startswith("#"):
                line = f.readline()
            elif line.startswith("hostIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        magicblastIndexH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("hostAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if not splits[0].startswith("#"):
                        magicblastMappingH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        magicblastIndexV += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if not line.startswith("#"):
                        magicblastMappingV += line.split("#")[0].strip() + " "
                    line = f.readline()
            else:
                line = f.readline()

    # ALIGNMENT - Minimap2
    with open(f"{configDir}minimap2.config", "r") as f:
        global minimapIndexH, minimapMappingH, minimapIndexV, minimapMappingV
        line = f.readline()
        while line:
            if line.startswith("#"):
                line = f.readline()
            elif line.startswith("hostIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        minimapIndexH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("hostAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if splits[0] == "-x":
                        minimapPresetH = splits[1]
                    elif not splits[0].startswith("#"):
                        minimapMappingH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        minimapIndexV += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if splits[0] == "-x":
                        minimapPresetV = splits[1]
                    elif not line.startswith("#"):
                        minimapMappingV += line.split("#")[0].strip() + " "
                    line = f.readline()
            else:
                line = f.readline()

    # ALIGNMENT - GMAP
    with open(f"{configDir}gmap.config", "r") as f:
        global gmapIndexH, gmapMappingH, gmapIndexV, gmapMappingV
        line = f.readline()
        while line:
            if line.startswith("#"):
                line = f.readline()
            elif line.startswith("hostIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        gmapIndexH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("hostAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if not splits[0].startswith("#"):
                        gmapMappingH += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenIndex:"):
                line = f.readline()
                while line.startswith(" "):
                    if not line.startswith("#"):
                        gmapIndexV += line.split("#")[0].strip() + " "
                    line = f.readline()
            elif line.startswith("pathogenAlignment:"):
                line = f.readline()
                while line.startswith(" "):
                    splits = line.split()
                    if not line.startswith("#"):
                        gmapMappingV += line.split("#")[0].strip() + " "
                    line = f.readline()
            else:
                line = f.readline()

    # CALLING - REDItools2
    with open(f"{configDir}reditools.config", "r") as f:
        line = f.readline()
        while line:
            if not line.startswith("#"):
                reditools += line.split("#")[0].strip() + " "
            line = f.readline()

    # CALLING - JACUSA
    with open(f"{configDir}jacusa.config", "r") as f:
        line = f.readline()
        while line:
            if not line.startswith("#"):
                jacusa += line.split("#")[0].strip() + " "
            line = f.readline()

    # CALLING - bcftools
    with open(f"{configDir}bcftools.config", "r") as f:
        line = f.readline()
        while line:
            if not line.startswith("#"):
                bcftools += line.split("#")[0].strip() + " "
            line = f.readline()

    # ALIGNMENT QUALITY - Qualimap
    with open(f"{configDir}qualimap.config", "r") as f:
        global qualimap
        line = f.readline()
        while line:
            if not line.startswith("#"):
                qualimap += line.split("#")[0].strip() + " "
            line = f.readline()