from jinja2 import Environment, FileSystemLoader
import config
from Bio import Entrez
import xml.etree.ElementTree as ET
import pathlib
import glob
import re
import os
import subprocess
import pandas as pd

def generateHtmlReport():
    sourcePath = os.getcwd()
    workPath = config.workPath
    os.chdir(workPath)
    workPath = "../" + os.path.basename(os.getcwd()) + "/"

    readsSource = config.source
    if readsSource == "project":
        readsSource = "List of BioProject accessions"
    if readsSource == "sra":
        readsSource = "List of SRA experiments"
    if readsSource == "file":
        readsSource = "List of files"

    totalSamples = 0
    inputs = []
    if config.source == "file":
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
                    file1 = files[0]
                    if not os.path.isabs(file1):
                        file1 = os.path.join(sourcePath, file1)
                        file1 = os.path.relpath(file1, os.getcwd())
                    files = [file1, ""]
                else:
                    file1 = files[0]
                    file2 = files[1]
                    if not os.path.isabs(file1):
                        file1 = os.path.join(sourcePath, file1)
                        file1 = os.path.relpath(file1, os.getcwd())
                        file2 = os.path.join(sourcePath, file2)
                        file2 = os.path.relpath(file2, os.getcwd())
                    files = [file1, file2]
                inputs.append((runId, runType, files[0], files[1]))
                totalSamples += 1

    elif config.source == "sra":
        with open("../sras.txt", "r") as f:
            # Skip the first line
            line = f.readline()
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip()
                sra = line.split()[0]
                files = sorted(glob.glob(f"{workPath}/data/fastq/{sra}*.fastq"))
                files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
                runType = ""
                if len(files) > 1: 
                    runType = "paired"
                    file1 = files[0]
                    file2 = files[1]
                    if not os.path.isabs(file1):
                        file1 = os.path.join(sourcePath, file1)
                        file1 = os.path.relpath(file1, os.getcwd())
                        file2 = os.path.join(sourcePath, file2)
                        file2 = os.path.relpath(file2, os.getcwd())
                    files = [file1, file2]
                elif len(files) == 1:
                    runType = "single"
                    file1 = files[0]
                    if not os.path.isabs(file1):
                        file1 = os.path.join(sourcePath, file1)
                        file1 = os.path.relpath(file1, os.getcwd())
                    files = [file1, ""]
                inputs.append((sra, runType, files[0], files[1]))
                totalSamples += 1

    elif config.source == "project":
        with open("../projects.txt", "r") as f:
            # Skip the first line
            line = f.readline()
            line = f.readline().strip()
            types = []
            ids = []
            while line != None and line != "":
                if line.startswith("#"):
                    line = f.readline().strip()
                    continue
                project = line.split()[0]
                # Search each project and get the list of experiment IDs
                search = Entrez.esearch(db="sra", term=f"{project}[BioProject]", rettype="gb", retmode="text", retmax="10000")
                read = Entrez.read(search)
                search.close()
                ids += read["IdList"]
                count = int(read["Count"])
                line = f.readline().strip()

            # Get sras for the experiments
            search = Entrez.epost(db='sra', id=",".join(ids), retmode='text', rettype='gb')
            read = Entrez.read(search)
            search.close()
            webenv = read["WebEnv"]
            key = read["QueryKey"]
            batchSize = 200
            sras = []
            for start in range(0, count, batchSize):
                search = Entrez.efetch(
                    db="sra",
                    retmode="xml",
                    retstart=start,
                    retmax=batchSize,
                    webenv=webenv,
                    query_key=key,
                    idtype="acc",
                )
                data = search.read()
                search.close()
                root = ET.fromstring(data)
                for el in root.iter("RUN_SET"):
                    for el2 in el.iter("RUN"):
                        totalSamples += 1
                        sra = el2.attrib["accession"]
                        files = glob.glob(f"{workPath}data/fastq/{sra}*.fastq")
                        files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
                        runType = ""
                        if len(files) > 1: 
                            runType = "paired"
                            file1 = files[0]
                            file2 = files[1]
                            if not os.path.isabs(file1):
                                file1 = os.path.join(sourcePath, file1)
                                file1 = os.path.relpath(file1, os.getcwd())
                                file2 = os.path.join(sourcePath, file2)
                                file2 = os.path.relpath(file2, os.getcwd())
                            files = [file1, file2]
                        elif len(files) == 1:
                            runType = "single"
                            file1 = files[0]
                            if not os.path.isabs(file1):
                                file1 = os.path.join(sourcePath, file1)
                                file1 = os.path.relpath(file1, os.getcwd())
                            files = [file1, ""]
                        inputs.append((sra, runType, files[0], files[1]))

    hostReference = config.hostReferencePath[0]
    if not os.path.isabs(hostReference):
        hostReference = os.path.join(sourcePath, hostReference)
        hostReference = os.path.relpath(hostReference, os.getcwd())

    pathogenReference = []
    for genome, geneFormat in zip(config.pathogenReferenceGenomePaths, config.pathogenReferenceGenesFormats):
        path = pathlib.Path(genome)
        parentFolder = path.parent
        if not os.path.isabs(parentFolder):
            parentFolder = os.path.join(sourcePath, parentFolder)
            parentFolder = os.path.relpath(parentFolder, os.getcwd())
        if geneFormat == "gtf22":
            geneFormat = "gtf"
        elif geneFormat == "gff3":
            geneFormat = "gff"
        elif geneFormat == "genbank":
            geneFormat = "gbk"
        elif geneFormat == "refSeq":
            geneFormat = "refseq"
        pathogenReference.append((f"{parentFolder}/genome.fa", f"{parentFolder}/protein.fa", f"{parentFolder}/genes.{geneFormat}"))

    totalSequences = []
    samFiles = sorted(glob.glob(f"2-alignment/host/sam/*/*_{config.alignmentSoftwareHost}.sam"))
    for sam in samFiles:
        name = os.path.splitext(os.path.basename(sam))[0]
        nameWithoutTool = name.replace("_" + config.alignmentSoftwareHost, "")

        preCount = 0
        postCount = 0
        
        command = f'{config.samtoolsPath} view {sam} | cut -f1 | sort -u | wc -l'
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if result.returncode == 0:
            preCount = int(result.stdout.strip())
        else:
            print(f"Error while running the command: {result.stderr}")
            preCount = -1


        bam = sam.replace("/sam/", "/bam/").replace(".sam", ".bam")
        command = f'{config.samtoolsPath} view {bam} | cut -f1 | sort -u | wc -l'
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if result.returncode == 0:
            postCount = int(result.stdout.strip())
        else:
            print(f"Error while running the command: {result.stderr}")
            postCount = -1

        # os.remove(samPost)

        totalSequences.append((nameWithoutTool, preCount, postCount))

    alignmentSoftware = {
        "star": ("STAR", "https://github.com/alexdobin/STAR", "This tool is a widely used RNA-seq read aligner for short and long reads, particularly well-suited for mapping reads to genomes with complex structures, such as those with many introns and alternative splicing events."),
        "hisat2": ("HISAT2", "http://daehwankimlab.github.io/hisat2/", "This tool is a widely used RNA-seq read aligner for short reads, particularly well-suited for ekaryotic transcriptomes with complex splicing patterns."),
        "bwa": ("BWA", "https://bio-bwa.sourceforge.net/", "This tool is a widely used read aligner made for short DNA reads. Since it is not splice-aware, it is not suggested for RNA-seq data."),
        "minimap2": ("Minimap2", "https://github.com/lh3/minimap2", "It is a versatile tool that can handle both DNA and RNA long reads. While it can handle short reads, it is optimized for long reads."),
        "magicblast": ("Magic-BLAST", "https://ncbi.github.io/magicblast/", "It is a versatile tool that can handle both DNA and RNA short and long reads. It is built on the BLAST algorithm."),
        "gmap": ("GMAP", "https://github.com/juliangehring/GMAP-GSNAP", "It is a tool that was especially developed for aligning cDNA reads, performing even better with long reads."),
    }

    cropSoftware = config.cropSoftware
    cropWeb = ""
    if cropSoftware == "trimmomatic":
        cropSoftware = "Trimmomatic"
        cropWeb = "http://www.usadellab.org/cms/?page=trimmomatic"
    elif cropSoftware == "trimgalore":
        cropSoftware = "Trim Galore"
        cropWeb = "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"

    errors = []
    for genome, _, _ in pathogenReference:
        path = pathlib.Path(genome)
        ref = path.parent.name
        for run, _, _, _ in inputs:
            files = sorted(glob.glob(f"{workPath}/3-qualimap/{ref}/{run}/genome_results.txt"))
            with open(f"{workPath}3-qualimap/{ref}/{run}/genome_results.txt") as f:
                for line in f:
                    if line.lstrip().startswith("general error rate"):
                        value = float('{:.2f}'.format(float(line.strip().split("=")[1].strip()) * 100))
                        errors.append((ref, run, value))

    for genome, _, _ in pathogenReference:
        path = pathlib.Path(genome)
        ref = path.parent.name
        df = pd.read_csv(f"{workPath}6-visualization/{ref}/csv/globalCommon.csv", nrows=100)
        mutationCountBarPlot = f"{workPath}6-visualization/{ref}/graphs/common.mutationCountBarPlot.png"
        mutationCountBoxPlot = f"{workPath}6-visualization/{ref}/graphs/common.mutationCountBoxPlot.png"
        mutationCountPerRunBarPlot = f"{workPath}6-visualization/{ref}/graphs/common.mutationsPerRunCountBarPlot.png"
        mutationCountPerRunBoxPlot = f"{workPath}6-visualization/{ref}/graphs/common.mutationsPerRunCountBoxPlot.png"
        frequencyPerMutationPlot = f"{workPath}6-visualization/{ref}/graphs/common.frequencyPerMutation.png"
        frequencyPerRunPlot = f"{workPath}6-visualization/{ref}/graphs/common.frequencyPerRun.png"

        geneMutationCountBarPlot = ""
        geneMutationCountBoxPlot = ""
        circosPlot = ""
        heatmapPlot = ""
        circosPerRunPlot = ""
        heatmapPerRunPlot = ""
        frequencyPerGenePlot = ""
        regressionPlot = ""
        histogramPlot = ""

        if config.figGenerateGlobalPresencePerRunPerPositionGraphs.lower() == "on":
            if config.figPositionGraph in ["both", "circos"]:
                circosPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/graphs/circos/*.png"))
                circosPlot = circosPlots[0]
            if config.figPositionGraph in ["both", "heatmap"]:
                heatmapPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/graphs/heatmap/*.png"))
                heatmapPlot = heatmapPlots[0]

        if config.figGenerateGlobalMutationCountPerGeneBarPlots.lower() == "on":
            geneMutationCountBarPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/graphs/geneBarPlot/*.png"))
            geneMutationCountBarPlot = geneMutationCountBarPlots[0]
        elif config.figGenerateMutationCountPerGeneBarPlotPerRun.lower() == "on":
            geneMutationCountBarPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/*/graphs/geneBarPlot/*.png"))
            geneMutationCountBarPlot = geneMutationCountBarPlots[0]
        
        if config.figGenerateGlobalMutationCountPerGeneBoxPlots.lower() == "on":
            geneMutationCountBoxPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/graphs/geneBoxPlot/*.png"))
            geneMutationCountBoxPlot = geneMutationCountBoxPlots[0]
        elif config.figGenerateMutationCountPerGeneBoxPlotPerRun.lower() == "on":
            geneMutationCountBoxPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/*/graphs/geneBoxPlot/*.png"))
            geneMutationCountBoxPlot = geneMutationCountBoxPlots[0]
            
        if config.figGenerateFrequencyPerMutationPerPositionGraphsPerRun.lower() == "on":
            if config.figPositionGraph in ["both", "circos"]:
                circosPerRunPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/*/graphs/circos/*.png"))
                circosPerRunPlot = circosPerRunPlots[0]
            if config.figPositionGraph in ["both", "heatmap"]:
                heatmapPerRunPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/*/graphs/heatmap/*.png"))
                heatmapPerRunPlot = heatmapPlots[0]
        
        if config.figGenerateGlobalFrequencyPerGeneStripPlots.lower() == "on":
            frequencyPerGenePlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/graphs/frequencyPerGene/*.png"))
            frequencyPerGenePlot = frequencyPerGenePlots[0]

        if config.figGenerateGlobalRegressionPlots.lower() == "on":
            regressionPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/graphs/*.regression.png"))
            regressionPlot = regressionPlots[0]

        if config.figGenerateGlobalDistributionHistogramsPlots.lower() == "on":
            histogramPlots = sorted(glob.glob(f"{workPath}6-visualization/{ref}/graphs/*.histogram.png"))
            histogramPlot = histogramPlots[0]

        break

    env = Environment(loader=FileSystemLoader(sourcePath))
    template = env.get_template(f'report_template.html')
    content = template.render(
        readsSource = readsSource,
        totalSamples =  totalSamples,
        resultsDirectory = workPath,
        referenceHostGenome = hostReference,
        hostAlignmentSoftware = alignmentSoftware[config.alignmentSoftwareHost],
        pathogenReference = pathogenReference,
        pathogenAlignmentSoftware = alignmentSoftware[config.alignmentSoftwarePathogen],
        inputs = inputs,
        cropSoftware = cropSoftware,
        cropWeb = cropWeb,
        cropQuality = config.cropMinMeanQuality,
        cropDecay = config.cropMaxDecay,
        cropSize = config.cropSize,
        hostAlignment = config.runHostAlignment.lower(),
        removeAligned = config.removeAlignedWithHost.lower(),
        totalSequences = totalSequences,
        qualimap = config.runQualimap.lower(),
        removeHighError = config.removeHighErrorRuns.lower(),
        minSNVCoverage = config.minSNVCoverage,
        minMainReadSupport = config.minMainReadSupport,
        minRecurringReadSupport = config.minRecurringReadSupport,
        minFrequency = config.minFrequency,
        maxAS_StrandOddsRatio = config.maxAS_StrandOddsRatio,
        callingReadMinQuality = config.callingReadMinQuality,
        callingBaseMinQuality = config.callingBaseMinQuality,
        maxError = float(config.qualimapMaxError * 100),
        errors = errors,
        callingSoftware = config.callingSoftware,
        df = df,
        figGenerateMutationCountBarPlotsPerRun = config.figGenerateMutationCountBarPlotsPerRun.lower(),
        figGenerateGlobalMutationCountBarPlots = config.figGenerateGlobalMutationCountBarPlots.lower(),
        figGenerateMutationCountBoxPlotsPerRun = config.figGenerateMutationCountBoxPlotsPerRun.lower(),
        figGenerateGlobalMutationCountBoxPlots = config.figGenerateGlobalMutationCountBoxPlots.lower(),
        figGenerateMutationCountPerGeneBarPlotPerRun = config.figGenerateMutationCountPerGeneBarPlotPerRun.lower(),
        figGenerateGlobalMutationCountPerGeneBarPlots = config.figGenerateGlobalMutationCountPerGeneBarPlots.lower(),
        figGenerateMutationCountPerGeneBoxPlotPerRun = config.figGenerateMutationCountPerGeneBoxPlotPerRun.lower(),
        figGenerateGlobalMutationCountPerGeneBoxPlots = config.figGenerateGlobalMutationCountPerGeneBoxPlots.lower(),
        figGenerateGlobalMutationCountPerRunBarPlots = config.figGenerateGlobalMutationCountPerRunBarPlots.lower(),
        figGenerateGlobalMutationCountPerRunBoxPlots = config.figGenerateGlobalMutationCountPerRunBoxPlots.lower(),
        figGenerateGlobalFrequencyPerMutationStripPlots = config.figGenerateGlobalFrequencyPerMutationStripPlots.lower(),
        figGenerateGlobalFrequencyPerGeneStripPlots = config.figGenerateGlobalFrequencyPerGeneStripPlots.lower(),
        figGenerateGlobalFrequencyPerRunStripPlots = config.figGenerateGlobalFrequencyPerRunStripPlots.lower(),
        figGenerateGlobalDistributionHistogramsPlots = config.figGenerateGlobalDistributionHistogramsPlots.lower(),
        figGenerateGlobalRegressionPlots = config.figGenerateGlobalRegressionPlots.lower(),
        figGenerateGlobalPresencePerRunPerPositionGraphs = config.figGenerateGlobalPresencePerRunPerPositionGraphs.lower(),
        figGenerateFrequencyPerMutationPerPositionGraphsPerRun = config.figGenerateFrequencyPerMutationPerPositionGraphsPerRun.lower(),
        
        mutationCountBarPlot = mutationCountBarPlot,
        mutationCountBoxPlot = mutationCountBoxPlot,
        mutationCountPerRunBarPlot = mutationCountPerRunBarPlot,
        mutationCountPerRunBoxPlot = mutationCountPerRunBoxPlot,
        frequencyPerMutationPlot = frequencyPerMutationPlot,
        frequencyPerRunPlot = frequencyPerRunPlot,
        figDistributionBinSize = config.figDistributionBinSize,
        histogramPlot = histogramPlot,
        geneMutationCountBarPlot = geneMutationCountBarPlot,
        figBarPlotPerGeneNumber = config.figBoxPlotPerGeneNumber,
        geneMutationCountBoxPlot = geneMutationCountBoxPlot,
        figPositionsPerCircos = config.figPositionsPerCircos,
        circosPlot = circosPlot,
        circosPerRunPlot = circosPerRunPlot,
        figPositionsPerHeatmap = config.figPositionsPerHeatmap,
        heatmapPlot = heatmapPlot,
        heatmapPerRunPlot = heatmapPerRunPlot,
        figStripPlotPerGeneNumber = config.figStripPlotPerGeneNumber,
        frequencyPerGenePlot = frequencyPerGenePlot,
        regressionPlot = regressionPlot,
        figPositionGraph = config.figPositionGraph, 
    )

    with open(f"{workPath}analysis_report.html", "w") as f:
        f.write(content)