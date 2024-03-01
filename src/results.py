"""This module handles everything about the final analysis and results 
output.
"""

import util
import config
import logger
import pandas as pd
import io
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib import gridspec
import glob
import pathlib
import warnings
import dask.dataframe as dd
import math
import multiprocessing

pd.options.mode.chained_assignment = None
plt.rcParams.update({'figure.max_open_warning': 0})
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

numProcesses = max(1, multiprocessing.cpu_count())

log = logger.logger

callingDir = config.workPath + "/4-snvCalling"
snpeffDir = config.workPath + "/5-snpeff"
resultsDir = config.workPath + "/6-visualization"

tab10 = ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", 
            "#17BECF"]
tab20 = ["#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", "#FF9896", "#9467BD", 
        "#C5B0D5", "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D",
        "#17BECF", "#9EDAE5"]

def mergeCalling(sras):
    """Merges the results from REDItools2, JACUSA, SnpEff and 
    bcftools mpileup, and generates the final CSV files and the 
    graphs.

    If configured to execute REDItools2, it reads the REDItools2 
    VCF and merges it with the SnpEff
    REDItools2 VCF file by the SNV position, so the SNVs found by 
    REDItools2 that were not found by SnpEff are discarded. Then, it 
    filters the SNVs by the minimum SNV coverage and the minimum read 
    support and saves those results in a different dataframe. These
    are saved into an CSV file named reditools.csv.
    
    If configured to execute JACUSA2, it reads the JACUSA VCF and 
    merges it with the SnpEff
    JACUSA VCF file by the SNV position, so the SNVs found by 
    JACUSA that were not found by SnpEff are discarded. Then, it 
    filters the SNVs by the minimum SNV coverage and the minimum read 
    support and saves those results in a different dataframe. These
    are saved into an CSV file named jacusa.csv.

    Finally, if both tools were executed, both filtered dataframes 
    are merged by the SNV position,
    so that only the SNVs common to both outputs remain. These
    are saved into an CSV file named runCommon.csv. If only one tool 
    was executed, then its respective output dataframe is also saved
    as runCommon.csv.

    These dataframes are also serialized as jacusa.h5, reditools.h5 
    and common.h5.

    The JACUSA, REDItools2 and common dataframes are 
    concatenated by dataframe type to the dataframes of the other
    runs. In the end, there will be a global JACUSA dataframe, a 
    global REDItools2 datafram and a global common dataframe.

    The JACUSA dataframe is saved into an CSV file named 
    globalJacusa.csv and jacusa.h5. Same thing happens with REDItools2
    (globalReditools.csv and reditools.h5) and common 
    (globalCommon.csv and common.h5) dataframes.

    All this process is repeated for every reference viral genome. The
    output files are found in the 6-visualization directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    reditoolsDir = callingDir + "/calling/reditools"
    jacusaDir = callingDir + "/calling/jacusa"
    depthsDir = callingDir + "/depths"
    util.makeDirectory(resultsDir)
    globalReditools = pd.DataFrame()
    globalJacusa = pd.DataFrame()
    for fullRef in config.pathogenReferenceGenomePaths:
        globalMain = None
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        util.makeDirectory(f"{resultsDir}/{ref}/csv")
        for sra in sras:
            run = sra[2]
            util.makeDirectory(f"{resultsDir}/{ref}/{run}/csv")
            log.info(f"Analyzing calling results from {run} vs {ref}...")
            
            files = glob.glob(f"{depthsDir}/{ref}/{run}_filtered.hdf")
            if len(files) == 0:
                log.error(f"AS_StrandOddsRatio filter for run with ID {run} against pathogen reference {ref} not found.")
                util.stopProgram()
            sor = dd.read_hdf(f"{depthsDir}/{ref}/{run}_filtered.hdf", key=run, chunksize=10000)
            sor["Concat"] = sor["CHROM"] + "||" + sor["Position"].astype(str)
            sor = sor.set_index("Concat")
            if config.callingSoftware in ["reditools", "both"]:
                files = glob.glob(f"{reditoolsDir}/{ref}/{run}.reditools.txt")
                if len(files) == 0:
                    log.error(f"REDitools calling for run with ID {run} against pathogen reference {ref} not found.")
                    util.stopProgram()
                reditools = _readReditoolsVcf(f"{reditoolsDir}/{ref}/{run}.reditools.vcf")
                reditools["RefReads"] = reditools.apply(_findRefReads, axis=1)
                reditools["AltReads"] = reditools.apply(_findAltReads, axis=1)
                reditools = reditools.set_index(["CHROM", "Position", "Alt"])
                reditoolsSnpeff = _readSnpEffVcf(f"{snpeffDir}/{ref}/{run}.snpeff.reditools.vcf")
                reditoolsSnpeff = reditoolsSnpeff.set_index(["CHROM", "Position", "ALT"])
                reditools = reditools.merge(reditoolsSnpeff, left_index=True, right_index=True) 
                # recReditools = reditools.copy()
                reditools = reditools.reset_index()
                reditools = reditools[reditools["Frequency"] >= config.minFrequency]
                reditools = reditools[reditools["TotalReads"] >= config.minSNVCoverage]
                reditools["FilteredSub"] = reditools.apply(lambda x: _filterMutationCt(x, config.minMainReadSupport), axis=1)
                reditools = reditools[reditools["FilteredSub"] != ""]
                reditools = reditools.drop(columns=["FilteredSub"])
                outputColumns = ["CHROM", "Position", "Reference", "Alt", "Type", "AAVar", "GeneName", "GeneID", "RefReads", "AltReads", "TotalReads", "Frequency", "A", "C", "G", "T"]
                reditools = reditools[outputColumns]
                
                # Save to file
                if len(glob.glob(f"{resultsDir}/{ref}/{run}/csv/reditools.csv")) > 0:
                    util.execCmd(f"rm {resultsDir}/{ref}/{run}/csv/reditools.csv")
                reditools.to_csv(f"{resultsDir}/{ref}/{run}/csv/reditools.csv", index=False)
                
                reditoolsRun = reditools.copy()
                reditoolsRun["Sample"] = run
                globalReditools = pd.concat([globalReditools, reditoolsRun], ignore_index=True)
                
                reditools["Concat"] = reditools["CHROM"] + "||" + reditools["Position"].astype(str)
                reditools = reditools.set_index("Concat") 
                reditools = dd.from_pandas(reditools, chunksize=10000)
                reditools = dd.merge(reditools, sor[[]], left_index=True, right_index=True, how="inner")
                reditools = reditools.compute()
                reditools = reditools.reset_index()
                reditools.drop(columns=["Concat"])

            if config.callingSoftware in ["jacusa", "both"]:
                files = glob.glob(f"{jacusaDir}/{ref}/{run}.jacusa.vcf")
                if len(files) == 0:
                    log.error(f"JACUSA calling for run with ID {run} against pathogen reference {ref} not found.")
                    util.stopProgram()
                jacusa = _readJacusaVcf(f"{jacusaDir}/{ref}/{run}.jacusa.vcf")
                jacusa = jacusa.set_index(["CHROM", "Position", "Alt"])
                jacusaSnpeff = _readSnpEffVcf(f"{snpeffDir}/{ref}/{run}.snpeff.jacusa.vcf")
                jacusaSnpeff = jacusaSnpeff.set_index(["CHROM", "Position", "ALT"])
                jacusa = jacusa.merge(jacusaSnpeff, left_index=True, right_index=True)     
                jacusa = jacusa.reset_index()
                jacusa["INFO2"] = jacusa["INFO2"].str.partition(":")[2]
                jacusa["BC"] = jacusa["INFO2"].str.partition(":")[0]
                jacusa["A"] = pd.to_numeric(jacusa["BC"].str.partition(",")[0])
                jacusa["BC"] = jacusa["BC"].str.partition(",")[2]
                jacusa["C"] = pd.to_numeric(jacusa["BC"].str.partition(",")[0])
                jacusa["BC"] = jacusa["BC"].str.partition(",")[2]
                jacusa["G"] = pd.to_numeric(jacusa["BC"].str.partition(",")[0])
                jacusa["T"] = pd.to_numeric(jacusa["BC"].str.partition(",")[2])
                jacusa["Frequency"] = jacusa.apply(_calculateFrequencyJacusa, axis=1)
                jacusa["RefReads"] = jacusa.apply(_findRefReads, axis=1)
                jacusa["AltReads"] = jacusa.apply(_findAltReads, axis=1)
                jacusa["TotalReads"] = jacusa["A"] + jacusa["C"] + jacusa["G"] + jacusa["T"]
                # recJacusa = jacusa.copy()
                jacusa = jacusa[jacusa["Frequency"] >= config.minFrequency]
                jacusa = jacusa[jacusa["TotalReads"] >= config.minSNVCoverage]
                jacusa["FilteredSub"] = jacusa.apply(lambda x: _filterMutationCt(x, config.minMainReadSupport), axis=1)
                jacusa = jacusa[jacusa["FilteredSub"] != ""]
                jacusa = jacusa.drop(columns=["FilteredSub"])
                outputColumns = ["CHROM", "Position", "Reference", "Alt", "Type", "AAVar", "GeneName", "GeneID", "RefReads", "AltReads", "TotalReads", "Frequency", "A", "C", "G", "T"]
                jacusa = jacusa[outputColumns]

                # Save to file
                if len(glob.glob(f"{resultsDir}/{ref}/{run}/csv/jacusa.csv")) > 0:
                    util.execCmd(f"rm {resultsDir}/{ref}/{run}/csv/jacusa.csv")
                jacusa.to_csv(f"{resultsDir}/{ref}/{run}/csv/jacusa.csv", index=False)
                jacusaRun = jacusa.copy()
                jacusaRun["Sample"] = run
                globalJacusa = pd.concat([globalJacusa, jacusaRun])
                
                jacusa["Concat"] = jacusa["CHROM"] + "||" + jacusa["Position"].astype(str)
                jacusa = jacusa.set_index("Concat") 
                jacusa = dd.from_pandas(jacusa, chunksize=10000)
                jacusa = dd.merge(jacusa, sor[[]], left_index=True, right_index=True, how="inner")
                jacusa = jacusa.compute()
                jacusa = jacusa.reset_index()
                jacusa.drop(columns=["Concat"])

            if config.callingSoftware == "both":
                # recJacusa = recJacusa.rename(columns={"A": "JacA", "C": "JacC", "G": "JacG", "T": "JacT", "TotalReads": "JacTotalReads", "RefReads": "JacRefReads", "AltReads": "JacAltReads", "Frequency": "JacFrequency"})
                jacusa = jacusa.rename(columns={"A": "JacA", "C": "JacC", "G": "JacG", "T": "JacT", "TotalReads": "JacTotalReads", "RefReads": "JacRefReads", "AltReads": "JacAltReads", "Frequency": "JacFrequency"})
                jacusa = jacusa.set_index(["CHROM", "Position", "Alt"])
                reditools = reditools.set_index(["CHROM", "Position", "Alt"])
                main = pd.merge(reditools, jacusa[["JacRefReads",  "JacAltReads", "JacTotalReads", "JacFrequency", "JacA", "JacC", "JacG", "JacT"]], left_index=True, right_index=True, how="inner")
                reditools = reditools.reset_index()
                jacusa = jacusa.reset_index()
                main = main.reset_index()
                # rec = dd.merge(recReditools, recJacusa[["CHROM", "Position", "Alt", "JacRefReads", "JacAltReads", "JacTotalReads", "JacFrequency", "JacA", "JacC", "JacG", "JacT"]], left_on=["CHROM", "Position", "Alt"], right_on=["CHROM", "Position", "Alt"], how="inner")
            elif config.callingSoftware == "jacusa":
                main = jacusa
                # rec = recJacusa 
            elif config.callingSoftware == "reditools":
                main = reditools
                # rec = recReditools

            main["Sample"] = run
            # rec["Sample"] = run

            if globalMain is None:
                globalMain = main.copy()
            else:
                globalMain = pd.concat([globalMain, main])

            # Save to files
            if config.callingSoftware == "both":
                outputColumns = outputColumns + ["JacRefReads", "JacAltReads", "JacTotalReads", "JacFrequency", "JacA", "JacC", "JacG", "JacT", "Sample"]
            else:
                outputColumns = outputColumns + ["Sample"]
            main = main[outputColumns]
            if len(glob.glob(f"{resultsDir}/{ref}/{run}/csv/runCommon.csv")) > 0:
                util.execCmd(f"rm {resultsDir}/{ref}/{run}/csv/runCommon.csv")
            main.to_csv(f"{resultsDir}/{ref}/{run}/csv/runCommon.csv", index=False)
            
            if config.callingSoftware == "both":
                reditools = reditools.assign(Mutation=lambda row: row.Reference + row.Alt)
                jacusa = jacusa.assign(Mutation=lambda row: row.Reference + row.Alt)
            main = main.assign(Mutation=lambda row: row.Reference + row.Alt)

            reditools.to_hdf(f"{resultsDir}/{ref}/{run}/reditools.h5", "key")
            jacusa.to_hdf(f"{resultsDir}/{ref}/{run}/jacusa.h5", "key")
            main.to_hdf(f"{resultsDir}/{ref}/{run}/common.h5", "key")

        globalMain = globalMain.drop_duplicates()
        positions = len(globalMain[["CHROM", "Position"]].drop_duplicates())

        log.info(f"{positions} positions were found.")
        log.info("Generating global files...")
        if config.callingSoftware == "both":
            if len(glob.glob(f"{resultsDir}/{ref}/csv/globalReditools.csv")) > 0:
                util.execCmd(f"rm {resultsDir}/{ref}/csv/globalReditools.csv")
            globalReditools.to_csv(f"{resultsDir}/{ref}/csv/globalReditools.csv", index=False)
            if len(glob.glob(f"{resultsDir}/{ref}/csv/globalJacusa.csv")) > 0:
                util.execCmd(f"rm {resultsDir}/{ref}/csv/globalJacusa.csv")
            globalJacusa.to_csv(f"{resultsDir}/{ref}/csv/globalJacusa.csv", index=False)
        if len(glob.glob(f"{resultsDir}/{ref}/csv/globalCommon.csv")) > 0:
            util.execCmd(f"rm {resultsDir}/{ref}/csv/globalCommon.csv")
        globalMain = globalMain.drop(columns="Concat")
        globalMain.to_csv(f"{resultsDir}/{ref}/csv/globalCommon.csv", index=False)

        globalMain.loc[:, ["Mutation"]] = globalMain["Reference"] + globalMain["Alt"]
        
        if config.callingSoftware == "both":
            globalReditools.loc[:, ["Mutation"]] = globalReditools["Reference"] + globalReditools["Alt"]
            globalJacusa.loc[:, ["Mutation"]] = globalJacusa["Reference"] + globalJacusa["Alt"]

        globalReditools.to_hdf(f"{resultsDir}/{ref}/reditools.h5", "key")
        globalJacusa.to_hdf(f"{resultsDir}/{ref}/jacusa.h5", "key")
        globalMain.to_hdf(f"{resultsDir}/{ref}/common.h5", "key")


def generateGraphs(sras):
    """Generates all the graphs that are configured to run in the main.config file. 

    It reads all the h5 files and produces the required graphs.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    for fullRef in config.pathogenReferenceGenomePaths:
        path = pathlib.Path(fullRef)
        ref = path.parent.name
        util.makeDirectory(f"{resultsDir}/{ref}/graphs")
        util.makeDirectory(f"{resultsDir}/{ref}/graphs/geneBarPlot")
        util.makeDirectory(f"{resultsDir}/{ref}/graphs/geneBoxPlot")
        util.makeDirectory(f"{resultsDir}/{ref}/graphs/frequencyPerGene")
        for sra in sras:
            run = sra[2]
            log.info(f"Generating graphs for {run} vs {ref}...")
            util.makeDirectory(f"{resultsDir}/{ref}/{run}/graphs")
            util.makeDirectory(f"{resultsDir}/{ref}/{run}/graphs/geneBarPlot")
            util.makeDirectory(f"{resultsDir}/{ref}/{run}/graphs/geneBoxPlot")

            reditools = pd.read_hdf(f"{resultsDir}/{ref}/{run}/reditools.h5", "key")
            jacusa = pd.read_hdf(f"{resultsDir}/{ref}/{run}/jacusa.h5", "key")
            main = pd.read_hdf(f"{resultsDir}/{ref}/{run}/common.h5", "key")
            if config.figGenerateMutationCountBarPlotsPerRun.lower() == "on":
                if config.callingSoftware == "both":
                    _graphMutationCountBarPlot(reditools, ref, "reditools", run)
                    _graphMutationCountBarPlot(jacusa, ref, "jacusa", run)
                elif config.callingSoftware == "jacusa":
                    _graphMutationCountBarPlot(jacusa, ref, "jacusa", run)
                elif config.callingSoftware == "reditools":
                    _graphMutationCountBarPlot(reditools, ref, "reditools", run)
                _graphMutationCountBarPlot(main, ref, "common", run)
            if config.figGenerateMutationCountPerGeneBarPlotPerRun.lower() == "on":
                _graphMutationsPerGeneCountBarPlot(main, ref, "common", run)
            if config.figGenerateMutationCountBoxPlotsPerRun.lower() == "on":
                _graphMutationCountBoxPlot(main, ref, "common", run)
            if config.figGenerateMutationCountPerGeneBoxPlotPerRun.lower() == "on":
                _graphMutationsPerGeneCountBoxPlot(main, ref, "common", run)
            if config.figGenerateFrequencyPerMutationPerPositionGraphsPerRun.lower() == "on":
                if config.figPositionGraph in ["both", "circos"]:
                    util.makeDirectory(f"{resultsDir}/{ref}/{run}/graphs/circos")
                    _graphCircosFrequencyPerMutation(main, ref, run)
                if config.figPositionGraph in ["both", "heatmap"]:
                    util.makeDirectory(f"{resultsDir}/{ref}/{run}/graphs/heatmap")
                    _graphHeatmapFrequencyPerMutation(main, ref, run)

        log.info(f"Generating global graphs for {ref}...")
        globalReditools = pd.read_hdf(f"{resultsDir}/{ref}/reditools.h5", "key")
        globalJacusa = pd.read_hdf(f"{resultsDir}/{ref}/jacusa.h5", "key")
        globalMain = pd.read_hdf(f"{resultsDir}/{ref}/common.h5", "key")
    
        if config.figGenerateGlobalMutationCountBarPlots.lower() == "on":
            if config.callingSoftware == "both":
                globalReditools = globalReditools
                globalJacusa = globalJacusa
                _graphMutationCountBarPlot(globalReditools, ref, "reditools")
                _graphMutationCountBarPlot(globalJacusa, ref, "jacusa")
            elif config.callingSoftware == "jacusa":
                globalJacusa = globalJacusa
                _graphMutationCountBarPlot(globalJacusa, ref, "jacusa")
            elif config.callingSoftware == "reditools":
                globalReditools = globalReditools
                _graphMutationCountBarPlot(globalReditools, ref, "reditools")
            _graphMutationCountBarPlot(globalMain, ref, "common")
        if config.figGenerateGlobalMutationCountPerGeneBarPlots.lower() == "on":
            _graphMutationsPerGeneCountBarPlot(globalMain, ref, "common")
        if config.figGenerateGlobalMutationCountPerRunBarPlots.lower() == "on":
            _graphMutationsPerRunCountBarPlot(globalMain, ref, "common")
        if config.figGenerateGlobalMutationCountBoxPlots.lower() == "on":
            _graphMutationCountBoxPlot(globalMain, ref, "common")
        if config.figGenerateGlobalMutationCountPerGeneBoxPlots.lower() == "on":
            _graphMutationsPerGeneCountBoxPlot(globalMain, ref, "common")
        if config.figGenerateGlobalMutationCountPerRunBoxPlots.lower() == "on":
            _graphMutationsPerRunCountBoxPlot(globalMain, ref, "common")
        if config.figGenerateGlobalFrequencyPerMutationStripPlots.lower() == "on":
            _graphFrequencyPerMutation(globalMain, ref, "common")
        if config.figGenerateGlobalFrequencyPerGeneStripPlots.lower() == "on":
            _graphFrequencyPerGene(globalMain, ref, "common")
        if config.figGenerateGlobalFrequencyPerRunStripPlots.lower() == "on":
            _graphFrequencyPerRun(globalMain, ref, "common")
        if config.figGenerateGlobalDistributionHistogramsPlots.lower() == "on":
            _graphHistograms(globalMain, ref, "common")
        if config.figGenerateGlobalRegressionPlots.lower() == "on":
            _graphRegression(globalMain, ref, "common")
        if config.figGenerateGlobalPresencePerRunPerPositionGraphs.lower() == "on":
            if config.figPositionGraph in ["both", "circos"]:
                util.makeDirectory(f"{resultsDir}/{ref}/graphs/circos")
                _graphCircosPresencePerRun(globalMain, ref)
            if config.figPositionGraph in ["both", "heatmap"]:
                util.makeDirectory(f"{resultsDir}/{ref}/graphs/heatmap")
                _graphHeatmapPresencePerRun(globalMain, ref)

def runSnpEff(sras):
    """Runs SnpEff on the VCF files resulting from JACUSA and REDItools2.

    First, it creates a configuration SnpEff file, adding the database
    names for every reference genome. Then, it creates those databases.
    Finally, each VCF file is analyzed against each database and outputs
    a VCF file per run per tool. These files' names end with
    .snpeff.jacusa.vcf and .snpeff.reditools.vcf, and are found in the 
    5-snpeff directory.

    Args:
        sras (list): List of tuples with the following data:
            * A list of the paths for the input run (one file if single-end, two if paired-end) 
            * Run type. "single" if single-end, "paired" if paired-end
            * Run ID  
    """
    util.makeDirectory(snpeffDir)
    reditoolsDir = callingDir + "/calling/reditools"
    jacusaDir = callingDir + "/calling/jacusa"
    util.execCmd(f"cp {pathlib.Path(config.snpEffPath).parent}/snpEff.config {config.workPath}/")
    with open(f"{config.workPath}/snpEff.config", "w") as f:
        for fullRef in config.pathogenReferenceGenomePaths:
            path = pathlib.Path(fullRef)
            ref = path.parent.name
            f.write(f"{ref}.genome : {ref}\n")
    jobs = []
    for fullRef, genesFormat in zip(config.pathogenReferenceGenomePaths, config.pathogenReferenceGenesFormats):
        util.makeDirectory(f"{snpeffDir}/{ref}")
        path = pathlib.Path(fullRef)
        ref = path.parent.name 
        util.execCmd(f"{config.javaPath} -jar {config.snpEffPath} build -{genesFormat} -c {config.workPath}/snpEff.config -v {ref}")
        for sra in sras:
            run = sra[2]
            if config.callingSoftware in ["reditools", "both"]:
                util.runCommand(f"{config.javaPath} -jar {config.snpEffPath} -c {config.workPath}/snpEff.config {ref} {reditoolsDir}/{ref}/{run}.reditools.presnpeff.vcf", outFile=f"{snpeffDir}/{ref}/{run}.snpeff.reditools.vcf", jobs=jobs, jobName="snpeff")
            if config.callingSoftware in ["jacusa", "both"]:
                util.runCommand(f"{config.javaPath} -jar {config.snpEffPath} -c {config.workPath}/snpEff.config {ref} {jacusaDir}/{ref}/{run}.jacusa.presnpeff.vcf", outFile=f"{snpeffDir}/{ref}/{run}.snpeff.jacusa.vcf", jobs=jobs, jobName="snpeff")
    util.waitForJobs(jobs)

def _filterMutationCt(row, min):
    """Verifies that the alternate allele has the minimum support.

    The minimum support value is provided in the main.config file,
    found in the config directory.

    Args:
        row (Series): Row of the VCF file.
        min (int): Minimum support value for the alternate allele.

    Returns:
        str: The alternate allele letter if its support value is 
        greater than the minimum. Empty string otherwise.
    """
    if row["Reference"] != "A" and row["Alt"] == "A" and row["A"] >= min:
        return "A"
    if row["Reference"] != "C" and row["Alt"] == "C" and row["C"] >= min:
        return "C"
    if row["Reference"] != "G" and row["Alt"] == "G" and row["G"] >= min:
        return "G"
    if row["Reference"] != "T" and row["Alt"] == "T" and row["T"] >= min:
        return "T"
    return ""

def _calculateFrequencyJacusa(row):
    """Calculates the frequency of the variant in
    each row of the JACUSA VCF files.

    It divides the support count for the alternate
    allele and divides that by the support count
    for the reference allele + alternate allele.

    Args:
        row (Series): Row of the VCF file.

    Returns:
        _type_: _description_
    """
    alt = row["Alt"]
    ct = row[alt]
    return round(float(ct / (row[row["Reference"]] + ct)), 6) * 100

def _readReditoolsVcf(path):
    """Transforms the REDItools2 VCF file into a Pandas DataFrame.

    Args:
        path (_type_): _description_

    Returns:
        DataFrame: The data read from the JACUSA VCF file.
    """
    lines = []
    with open(path, 'r') as f:
        start = False
        for line in f:
            if start:
                lines.append(line.strip())
            if line.startswith("#CHROM"):
                lines.append(line.strip())
                start = True
    csv = pd.read_csv(
        io.StringIO('\n'.join(lines)),
        dtype={'#CHROM': str, 'POSITION': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 
            'VARIANT_READ': int, 'TOTAL_READ': int, 'FREQUENCY': float, "A": int, "C": int, "T": int, "G": int},
        delimiter="\t"
    ).rename(columns={'#CHROM': 'CHROM', "11": "INFO2", "REF": "Reference", "POSITION": "Position", "ALT": "Alt", "FREQUENCY": "Frequency", "TOTAL_READ": "TotalReads", "VARIANT_READ": "AltReads"})
    
    return csv

def _readJacusaVcf(path):
    """Transforms the JACUSA VCF file into a Pandas DataFrame.

    Args:
        path (_type_): _description_

    Returns:
        DataFrame: The data read from the JACUSA VCF file.
    """
    lines = []
    with open(path, 'r') as f:
        start = False
        for line in f:
            if start:
                lines.append(line.strip())
            if line.startswith("#CHROM"):
                lines.append(line.strip())
                start = True
    csv = pd.read_csv(
        io.StringIO('\n'.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str, '11': str},
        delimiter="\t"
    ).rename(columns={'#CHROM': 'CHROM', "11": "INFO2", "REF": "Reference", "POS": "Position", "ALT": "Alt"})
    
    return csv

def _readSnpEffVcf(path):
    """Transforms the SnpEff VCF file into a Pandas DataFrame.

    Args:
        path (_type_): _description_

    Returns:
        DataFrame: The data read from the JACUSA VCF file.
    """
    lines = []
    counter = 1
    with open(path, 'r') as f:
        start = False
        for line in f:
            counter += 1
            if line.startswith("#CHROM"):
                start = True
                break
    csv = pd.read_csv(
        path,
        names = ["CHROM", "Position", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
        dtype={'CHROM': str, 'Position': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str},
        delimiter="\t", skiprows=counter
    )
    csv["Type"], csv["AAVar"], csv["GeneName"], csv["GeneID"] = zip(*csv.apply(_getSnpeffTypeAAVarGene, axis=1))
    csv = csv.drop(columns=['ID', 'REF', 'QUAL', 'FILTER', 'INFO'])
    return csv

def _getSnpeffTypeAAVarGene(row):
    """Gets the type of effect of the mutation on the amino 
    acid, the amino acid variant in the given row of a SnpEff 
    VCF file, the gene name and the gene ID.

    Args:
        row (Series): Row of a SnpEff VCF file.

    Returns:
        DataFrame: The type of effect and the amino acid variant.
    """
    varType = ""
    aaVar = ""
    geneName = ""
    geneID = ""
    alt = row["ALT"]
    info = row["INFO"]
    annotations = info.split("=")[1].split(",")
    for annotation in annotations:
        allele = annotation[0]
        rightSide = annotation[2:]
        if alt == allele:
            splits = rightSide.split("|")
            varType = splits[0]
            geneName = splits[2]
            geneID = splits[3]
            aaVar = splits[9]
            break
    return varType, aaVar, geneName, geneID

def _checkReditoolsSnpeffTypeExists(row):
    """Checks whether a row in the SnpEff VCF file for REDItools2
    has a type of effect.

    If the row has no type of effect, the alternate allele will become
    empty.

    Args:
        row (Series): Row of a SnpEff VCF file.

    Returns:
        Series: Modified row only with the alleles that contain
            a type of effect.
    """
    alts = row["Alt"]
    types = row["Type"]
    newAlts = []
    newTypes = []
    for alt in alts.split(","):
        for tp in types.split(","):
            if alt[1] == tp[0]:
                newAlts.append(alt[1])
                newTypes.append(tp[2:])
                break
    row["Alt"] = ",".join(newAlts)
    row["Type"] = ",".join(newTypes)
    return row

def _checkJacusaSnpeffTypeExists(row):
    """Checks whether a row in the SnpEff VCF file for JACUSA
    has a type of effect.

    If the row has no type of effect, the alternate allele will become
    empty.

    Args:
        row (Series): Row of a SnpEff VCF file.

    Returns:
        Series: Modified row only with the alleles that contain
            a type of effect.
    """
    alts = row["Alt"]
    types = row["Type"]
    newAlts = []
    newTypes = []
    for alt in alts.split(","):
        for tp in types.split(","):
            if alt == tp[0]:
                newAlts.append(alt)
                newTypes.append(tp[2:])
                break
    row["Alt"] = ",".join(newAlts)
    row["Type"] = ",".join(newTypes)
    return row

def _findAltReads(row):
    """Finds the support count of the alternate allele.

    Args:
        row (Series): Row of a VCF file.

    Returns:
        int: Support count of the alternate allele.
    """
    alt = row["Alt"]
    if alt == "A":
        return row["A"]
    if alt == "C":
        return row["C"]
    if alt == "G":
        return row["G"]
    if alt == "T":
        return row["T"]

def _findRefReads(row):
    """Finds the support count of the reference allele.

    Args:
        row (Series): Row of a VCF file.

    Returns:
        int: Support count of the reference allele.
    """
    alt = row["Reference"]
    if alt == "A":
        return row["A"]
    if alt == "C":
        return row["C"]
    if alt == "G":
        return row["G"]
    if alt == "T":
        return row["T"]

def _graphMutationCountBarPlot(df, ref, source, run=None):
    """Generates a bar plot with the counts per mutation for the given 
    DataFrame.

    It gets the count per possible mutation and finds the mean and 
    median of the reads for both Jacusa and REDItools. Then, it plots
    the counts on a bar plot, and adds the read mean and median on the
    lower side of the graph.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (_type_): _description_
        run (_type_, optional): _description_. Defaults to None.
    """
    height = config.figBarPlotCountHeight
    width = config.figBarPlotPerRunWidth
    counts = {"AC": 0, "AG": 0, "AT": 0, "CA": 0, "CG": 0, "CT": 0, "GA": 0, "GC": 0, "GT": 0, "TA": 0, "TC": 0, "TG": 0}
    for _, row in df.iterrows():
        counts[row["Mutation"]] += 1
    if source == "jacusa" and run is not None:
        mean = df["JacAltReads"].mean()
        median = df["JacAltReads"].median()
    else:
        mean = df["AltReads"].mean()
        median = df["AltReads"].median()
    mean = "{:.2f}".format(mean)
    median = "{:.2f}".format(median)
    ax = sns.barplot(x=list(counts.keys()), y=list(counts.values()))
    ax.set_title(f"Global SNV count per mutation for {ref}")
    ax.set_xlabel("Mutation")
    ax.set_ylabel("SNV count")
    ax.text(0.02, -0.1, f"Mean coverage: {mean}\nMedian coverage: {median}", ha="left", va="top", transform=ax.transAxes)
    fig = ax.get_figure()
    fig.set_size_inches(width, height)
    fileName = source + ".mutationCountBarPlot"
    if run is None:
        fig.savefig(f"{resultsDir}/{ref}/graphs/{fileName}.png", dpi=config.figDPI, bbox_inches="tight") 
    else:
        fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/{fileName}.png", dpi=config.figDPI, bbox_inches="tight") 
    fig.clf()

def _graphMutationsPerGeneCountBarPlot(df, ref, source, run=None):
    """Generates a bar plot with the counts per mutation per gene for 
    the given DataFrame.

    It gets the count per possible mutation per gene. Then, it plots
    the counts on a bar plot.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (_type_): _description_
        run (_type_, optional): _description_. Defaults to None.
    """
    mutations = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    groupSize = config.figBarPlotPerGeneNumber
    height = config.figBarPlotPerGeneHeight
    width = config.figBarPlotPerGeneWidth
    axisLabelSize = config.figBarPlotPerGeneAxisLabelSize
    tickSize = config.figBarPlotPerGeneTickSize
    titleSize = config.figBarPlotPerGeneTitleSize

    if len(df.GeneName.unique()) > groupSize:
        filePath = f"{resultsDir}/{ref}/graphs/geneBarPlot/_groups.txt"
        if run is not None:
            filePath = f"{resultsDir}/{ref}/{run}/graphs/geneBarPlot/_groups.txt"
        with open(filePath, "w") as f:
            f.write("Chromosome\tGeneName\tGroup")
        for chrom in df.CHROM.unique():
            countsPerMutation = dict.fromkeys(mutations, 0)
            counts = dict.fromkeys(df["GeneName"])
            for key in counts:
                counts[key] = countsPerMutation.copy()
            for _, row in df[df.CHROM == chrom].iterrows():
                counts[row["GeneName"]][row["Mutation"]] += 1
            
            genePosDf = df[df["CHROM"] == chrom][["GeneName", "Position"]]
            genePosDf = genePosDf.merge(genePosDf.groupby("GeneName").min().rename(columns={"Position": "GenePos"}), left_on="GeneName", right_index=True)
            genePosDf = genePosDf.set_index("GeneName")
            
            dfChrom = pd.DataFrame.from_dict(counts, orient="index")
            dfChrom = dfChrom.merge(genePosDf, left_index=True, right_index=True)
            dfChrom = dfChrom.reset_index().rename(columns={"index": "GeneName"})
            dfChrom = dfChrom.sort_values("GenePos", ascending=True)
            dfChrom = dfChrom.drop(columns=["Position", "GenePos"]).drop_duplicates()

            for group in range(math.ceil(len(dfChrom) / groupSize)):
                dfTemp = dfChrom.iloc[group * groupSize:(group+1) * groupSize]
                
                fig, ax = plt.subplots()
                fig.set_size_inches(width, height)
                barWidth = 0.5

                i = 0
                sum = 0
                for mutation in dfTemp.columns[1:]:
                    data = dfTemp.loc[:,mutation][::-1]
                    ax.barh(dfTemp.GeneName[::-1], data, barWidth, left=sum, color = tab20[i])
                    sum = np.add(sum, data)
                    i += 1

                ax.margins(0.01)
                ax.legend(dfTemp.columns[1:], loc="upper left", bbox_to_anchor=(1.05, 1.0))
                ax.set_title(f"SNV count per gene for {ref} - Group {group}", fontsize=titleSize)
                ax.set_xlabel("SNV Count", fontsize=axisLabelSize)
                ax.set_ylabel("Gene", fontsize=axisLabelSize)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_position(('outward', 8))
                ax.spines['bottom'].set_position(('outward', 5))
                plt.xticks(fontsize=tickSize) 
                plt.yticks(fontsize=tickSize) 
                fig.tight_layout(pad=1)
                fileName = f"{chrom}_{group}"
                if run is None:
                    fig.savefig(f"{resultsDir}/{ref}/graphs/geneBarPlot/{fileName}.png", bbox_inches='tight') 
                else:
                    fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/geneBarPlot/{fileName}.png", bbox_inches='tight') 
                fig.clf()
                # plt.gcf().set_size_inches(0, 0)

                with open(filePath, "a") as f:
                    for gene in dfTemp.GeneName.unique():
                        f.write("\n" + chrom + "\t" + gene + "\t" + str(int(group)))
    else:
        countsPerMutation = dict.fromkeys(mutations, 0)
        counts = dict.fromkeys(df["GeneName"])
        for key in counts:
            counts[key] = countsPerMutation.copy()
        for _, row in df.iterrows():
            counts[row["GeneName"]][row["Mutation"]] += 1
        
        genePosDf = df[["CHROM", "GeneName", "Position"]]
        genePosDf = genePosDf.merge(genePosDf.groupby(["CHROM", "GeneName"]).min().rename(columns={"Position": "GenePos"}), left_on=["CHROM", "GeneName"], right_index=True)
        genePosDf = genePosDf.set_index("GeneName")
        
        dfChrom = pd.DataFrame.from_dict(counts, orient="index")
        dfChrom = dfChrom.merge(genePosDf, left_index=True, right_index=True)
        dfChrom = dfChrom.reset_index().rename(columns={"index": "GeneName"})
        dfChrom = dfChrom.sort_values(["CHROM", "GenePos"], ascending=True)
        dfChrom = dfChrom.drop(columns=["Position", "GenePos"]).drop_duplicates()

        for group in range(math.ceil(len(dfChrom) / groupSize)):
            dfTemp = dfChrom.iloc[group * groupSize:(group+1) * groupSize]
            
            fig, ax = plt.subplots()
            barWidth = 0.5

            i = 0
            total = 0
            for mutation in dfTemp.columns[1:-1]:
                data = dfTemp.loc[:,mutation][::-1]
                ax.barh(dfTemp.GeneName[::-1], data, barWidth, left=total, color = tab20[i])
                total = np.add(total, data)
                i += 1

            ax.margins(0.01)
            ax.legend(dfTemp.columns[1:], loc="upper left", bbox_to_anchor=(1.05, 1.0))
            ax.set_title(f"SNV count per gene for {ref}", fontsize=titleSize)
            ax.set_xlabel("SNV Count", fontsize=axisLabelSize)
            ax.set_ylabel("Gene", fontsize=axisLabelSize)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_position(('outward', 8))
            ax.spines['bottom'].set_position(('outward', 5))
            plt.xticks(fontsize=tickSize) 
            plt.yticks(fontsize=tickSize) 
            fig.set_size_inches(width, height)
            fig.tight_layout(pad=1)
            if run is None:
                fileName = f"{ref}"
                fig.savefig(f"{resultsDir}/{ref}/graphs/geneBarPlot/{fileName}.mutationsPerGeneCountBarPlot.png", bbox_inches='tight') 
            else:
                fileName = f"common"
                fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/geneBarPlot/{fileName}.mutationsPerGeneCountBarPlot.png", bbox_inches='tight') 
            fig.clf()

def _graphMutationsPerRunCountBarPlot(df, ref, source, run=None):
    """Generates a bar plot with the counts per mutation per run for 
    the given DataFrame.

    It gets the count per possible mutation per run. Then, it plots
    the counts on a bar plot.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (_type_): _description_
        run (_type_, optional): _description_. Defaults to None.
    """
    barHeight = config.figBarPlotPerRunBarHeight
    width = config.figBarPlotPerRunWidth
    run = None
    source = "common"
    counts = dict.fromkeys(df["Sample"], 0)
    for _, row in df.iterrows():
        counts[row["Sample"]] += 1
    if source == "jacusa" and run is not None:
        mean = df["JacAltReads"].mean()
        median = df["JacAltReads"].median()
    else:
        mean = df["AltReads"].mean()
        median = df["AltReads"].median()
    mean = "{:.2f}".format(mean)
    median = "{:.2f}".format(median)
    keys = list(counts.keys())
    values = list(counts.values())
    fig, ax = plt.subplots(figsize=(width, len(keys) * barHeight))
    sns.barplot(y=keys, x=values)
    ax.set_title(f"Global SNV count per run for {ref}")
    ax.set_xlabel("Sample")
    ax.set_ylabel("SNV count")
    ax.tick_params(labelrotation=45)
    plt.draw()
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    fig = ax.get_figure()
    fileName = source + ".mutationsPerRunCountBarPlot"
    fig.savefig(f"{resultsDir}/{ref}/graphs/{fileName}.png", dpi=config.figDPI, bbox_inches="tight") 
    fig.clf()

def _graphMutationCountBoxPlot(df, ref, source, run=None):
    """Generates a box plot with the counts per mutation for the given 
    DataFrame.
    
    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (_type_): _description_
        run (_type_, optional): _description_. Defaults to None.
    """
    height = config.figBoxPlotCountHeight
    width = config.figBoxPlotPerRunWidth
    ax = sns.boxplot(x="Mutation", y="AltReads", data=df,
            order=["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"])
    ax.set_title(f"Mutated counts for {ref}")
    ax.set_xlabel("Mutation")
    ax.set_ylabel("Mutated reads per SNV")
    fig = ax.get_figure()
    fig.set_size_inches(width, height)
    for j in range(0, len(df["Mutation"].unique()), 2):
        plt.axvspan(j - 0.5, j+0.5, facecolor='0.2', alpha=0.1)
    fileName = source + ".mutationCountBoxPlot"
    if run is None:
        fig.savefig(f"{resultsDir}/{ref}/graphs/{fileName}.png", dpi=config.figDPI, bbox_inches="tight") 
    else:
        fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/{fileName}.png", dpi=config.figDPI, bbox_inches="tight") 
    fig.clf()

def _graphMutationsPerGeneCountBoxPlot(df, ref, source, run=None):
    """Generates a box plot with the mutation counts per gene per 
    run for the given DataFrame.

    It plots the mutation counts on a box plot.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (_type_): _description_
        run (_type_, optional): _description_. Defaults to None.
    """
    groupSize = config.figBoxPlotPerGeneNumber
    height = config.figBoxPlotPerGeneHeight
    width = config.figBoxPlotPerGeneWidth
    axisLabelSize = config.figBoxPlotPerGeneAxisLabelSize
    tickSize = config.figBoxPlotPerGeneTickSize
    titleSize = config.figBoxPlotPerGeneTitleSize
    
    genePosDf = df[["CHROM", "GeneName", "Position"]]
    genePosDf = genePosDf.merge(genePosDf.groupby(["CHROM", "GeneName"]).min().rename(columns={"Position": "GenePos"}), left_on=["CHROM", "GeneName"], right_index=True)
    genePosDf = genePosDf[["GeneName", "GenePos"]].drop_duplicates()
    genePosDf = genePosDf.set_index("GeneName")
    df = df.set_index("GeneName")
    df = df.merge(genePosDf, left_index=True, right_index=True)
    df = df.reset_index()
    df = df.sort_values("GenePos")
    df = df.set_index("GeneName")

    if len(df.index.unique()) > groupSize:
        filePath = f"{resultsDir}/{ref}/graphs/geneBoxPlot/_groups.txt"
        if run is not None:
            filePath = f"{resultsDir}/{ref}/{run}/graphs/geneBoxPlot/_groups.txt"
        with open(filePath, "w") as f:
            f.write("Chromosome\tGeneName\tGroup")
        for chrom in df.CHROM.unique():
            genes = df[df.CHROM == chrom].index.unique()
            for i in range(math.ceil(len(genes) / groupSize)):
                plt.figure(figsize=(width, height))
                data = df.loc[genes[i * groupSize:(i+1) * groupSize]]
                data = data.reset_index()
                ax = sns.boxplot(y="GeneName", x="AltReads", data=data, width=0.5)
                ax.set_title(f"SNV count per gene for {ref} - Group {i}", fontsize=titleSize)
                ax.set_xlabel("Mutation count", fontsize=axisLabelSize)
                ax.set_ylabel("Gene Name", fontsize=axisLabelSize)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_position(('outward', 8))
                ax.spines['bottom'].set_position(('outward', 5))
                ax.margins(0.01, tight=False)
                for j in range(0, len(data["GeneName"].unique()), 2):
                    plt.axhspan(j - 0.5, j+0.5, facecolor='0.2', alpha=0.1)
                fig = ax.get_figure()
                fig.tight_layout(pad=1)
                fig.set_size_inches(width, height)
                plt.xticks(fontsize=tickSize) 
                plt.yticks(fontsize=tickSize) 
                fileName = f"{chrom}_{i}"
                if run is None:
                    fig.savefig(f"{resultsDir}/{ref}/graphs/geneBoxPlot/{fileName}.png", bbox_inches='tight') 
                else:
                    fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/geneBoxPlot/{fileName}.png", bbox_inches='tight') 
                fig.clf()
                # plt.gcf().set_size_inches(0, 0)
                with open(filePath, "a") as f:
                    for gene in genes:
                        f.write("\n" + chrom + "\t" + gene + "\t" + str(int(i)))
    else:
        genes = df.index.unique()
        for i in range(math.ceil(len(genes) / groupSize)):
            plt.figure(figsize=(width, height))
            data = df.loc[genes[i * groupSize:(i+1) * groupSize]]
            data = data.reset_index()
            ax = sns.boxplot(y="GeneName", x="AltReads", data=data, width=0.5)
            ax.set_title(f"Mutation counts per SNV for {ref}", fontsize=titleSize)
            ax.set_xlabel("Mutation count", fontsize=axisLabelSize)
            ax.set_ylabel("Gene Name", fontsize=axisLabelSize)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_position(('outward', 8))
            ax.spines['bottom'].set_position(('outward', 5))
            ax.margins(0.01, tight=False)
            for j in range(0, len(data["GeneName"].unique()), 2):
                plt.axhspan(j - 0.5, j+0.5, facecolor='0.2', alpha=0.1)
            fig = ax.get_figure()
            fig.tight_layout(pad=1)
            fig.set_size_inches(width, height)
            plt.xticks(fontsize=tickSize) 
            plt.yticks(fontsize=tickSize) 
            if run is None:
                fileName = f"{ref}"
                fig.savefig(f"{resultsDir}/{ref}/graphs/geneBoxPlot/{fileName}.mutationsPerGeneCountBoxPlot.png", bbox_inches='tight') 
            else:
                fileName = f"common"
                fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/geneBoxPlot/{fileName}.mutationsPerGeneCountBoxPlot.png", bbox_inches='tight') 
            fig.clf()
            # plt.gcf().set_size_inches(0, 0)

def _graphMutationsPerRunCountBoxPlot(df, ref, source, run=None):
    """Generates a box plot with the counts per mutation per run for 
    the given DataFrame.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (_type_): _description_
        run (_type_, optional): _description_. Defaults to None.
    """
    barHeight = config.figBoxPlotPerRunBarHeight
    width = config.figBoxPlotPerRunWidth
    ax = sns.boxplot(y="Sample", x="AltReads", data=df)
    ax.set_title(f"Mutation counts for {ref}")
    ax.set_ylabel("Sample")
    ax.set_xlabel("Mutation count per SNV")
    ax.tick_params(labelrotation=45)
    fig = ax.get_figure()
    fig.set_size_inches(width, len(df.Sample.unique()) * barHeight)
    for j in range(0, len(df["Sample"].unique()), 2):
        plt.axhspan(j - 0.5, j+0.5, facecolor='0.2', alpha=0.1)
    fileName = source + ".mutationsPerRunCountBoxPlot"
    if run is None:
        fig.savefig(f"{resultsDir}/{ref}/graphs/{fileName}.png", bbox_inches='tight') 
    else:
        fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/{fileName}.png", bbox_inches='tight') 
    fig.clf()

def _graphFrequencyPerMutation(df, ref, source=""):
    """Generates a strip plot with the SNV frequency (%) of each 
    row in the VCF file, each strip being a different mutation.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (str, optional): _description_. Defaults to "".
    """
    height = config.figStripPlotPerMutationHeight
    width = config.figStripPlotPerMutationWidth
    dotSize = config.figStripPlotDotSize
    df["Mutation"] = df["Reference"] + df["Alt"]
    ax = sns.stripplot(x="Mutation", y="Frequency", data=df, hue="Mutation", size=dotSize,
        order = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"])
    ax.set_title(f"SNV frequency (%) for {ref}")
    ax.set_xlabel("Mutation")
    ax.set_ylabel("SNV frequency (%)")
    ax.legend_.remove()
    fig = ax.get_figure()
    fig.set_size_inches(width, height)
    fig.savefig(f"{resultsDir}/{ref}/graphs/{source}.frequencyPerMutation.png", dpi=config.figDPI, bbox_inches="tight") 
    fig.clf()

def _graphFrequencyPerGene(df, ref, source=""):
    """Generates a strip plot with the SNV frequency (%) of each 
    row in the VCF file, each strip being a different mutation.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (str, optional): _description_. Defaults to "".
    """
    dotSize = config.figStripPlotDotSize
    height = config.figStripPlotPerGeneHeight
    width = config.figStripPlotPerGeneWidth
    groupSize = config.figStripPlotPerGeneNumber
    titleSize = config.figStripPlotPerGeneTitleSize
    genePosDf = df[["GeneName", "Position"]]
    genePosDf = genePosDf.merge(genePosDf.groupby("GeneName").min().rename(columns={"Position": "GenePos"}), left_on="GeneName", right_index=True)
    genePosDf = genePosDf[["GeneName", "GenePos"]].drop_duplicates()
    genePosDf = genePosDf.set_index("GeneName")
    df = df.set_index("GeneName")
    df = df.merge(genePosDf, left_index=True, right_index=True)
    df = df.reset_index()
    df = df.sort_values("GenePos")
    df = df.set_index("GeneName")
    df["Mutation"] = df["Reference"] + df["Alt"]

    filePath = f"{resultsDir}/{ref}/graphs/frequencyPerGene/_groups.txt"
    with open(filePath, "w") as f:
        f.write("Chromosome\tGeneName\tGroup")

    for chrom in df.CHROM.unique():
        genes = df[df.CHROM == chrom].index.unique()
        for i in range(math.ceil(len(genes) / groupSize)):
            plt.figure(figsize=(width, height))
            data = df.loc[genes[i * groupSize:(i+1) * groupSize]]
            ax = sns.stripplot(y="GeneName", x="Frequency", data=data, hue="GeneName", size=dotSize)
            ax.set_title(f"SNV frequency (%) for {ref}  - Group {i}", fontsize=titleSize)
            ax.set_ylabel("Gene")
            ax.set_xlabel("SNV frequency (%)")
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_position(('outward', 8))
            ax.spines['bottom'].set_position(('outward', 5))
            ax.margins(0.01, tight=False)
            for j in range(0, len(data.index.unique()), 2):
                plt.axhspan(j - 0.5, j+0.5, facecolor='0.2', alpha=0.1)
            ax.legend_.remove()
            fig = ax.get_figure()
            fig.tight_layout(pad=0)
            fig.set_size_inches(width, height)
            fig.savefig(f"{resultsDir}/{ref}/graphs/frequencyPerGene/{chrom}_{i}.png", bbox_inches='tight') 
            fig.clf()

            with open(filePath, "a") as f:
                for gene in data.index.unique():
                    f.write("\n" + chrom + "\t" + gene + "\t" + str(int(i)))

def _graphFrequencyPerRun(df, ref, source=""):
    """Generates a strip plot with the SNV frequency (%) of each 
    row in the VCF file, each strip being a different run.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (str, optional): _description_. Defaults to "".
    """
    barHeight = config.figStripPlotPerRunBarHeight
    width = config.figStripPlotPerRunWidth
    dotSize = config.figStripPlotDotSize
    df["Mutation"] = df["Reference"] + df["Alt"]
    ax = sns.stripplot(y="Sample", x="Frequency", data=df, hue="Sample", size=dotSize)
    ax.set_title(f"SNV frequency (%) for {ref}")
    ax.set_xlabel("Sample")
    ax.set_ylabel("SNV frequency (%)")
    ax.tick_params(axis="x", labelrotation=45)
    ax.legend_.remove()
    fig = ax.get_figure()
    fig.set_size_inches(width, len(df["Sample"].unique()) * barHeight)
    fig.savefig(f"{resultsDir}/{ref}/graphs/{source}.frequencyPerRun.png", dpi=config.figDPI, bbox_inches="tight") 
    fig.clf()

def _graphHistograms(df, ref, source=""):
    """It plots a distribution histogram with the 6 possible mutations 
    and their reverse complement.

    First, it splits the locations of the SNVs into bins, and counts
    the number of times an SNV falls into each bin. Then, it draws the 
    bar plot in a way that the mutation goes upwards, and the reverse
    complement goes downwards.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (str, optional): _description_. Defaults to "".
    """
    for chrom in df.CHROM.unique():
        dfChrom = df[df.CHROM == chrom]
        maxDf = dfChrom["Position"].max()
        binSize = config.figDistributionBinSize
        maxBin = int(maxDf / binSize) 
        if maxDf % binSize != 0: 
            maxBin += 1
        keys = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
        dfChrom["Bin"] = dfChrom["Position"] // binSize
        countsDf = dfChrom[["Mutation", "Bin"]].groupby(["Mutation", "Bin"]).size().reset_index(name='Count')
        counts = {key: {bin: 0 for bin in np.arange(0, maxBin, 1)} for key in keys}
        for i in range(len(countsDf)):
            row = countsDf.iloc[i]
            counts[row.Mutation][row.Bin] = row.Count
        
        for key in ["TC", "GA", "TG", "TA", "GT", "GC"]:
            for bin in counts[key]:
                counts[key][bin] *= -1

        plt.figure()
        fig, axs = plt.subplots(3, 2)
        pairs = [
            ("AG", "TC", "#6f3777", "#bf5fcc", axs[0,1]), 
            ("CT", "GA", "#774f1b", "#cc882f", axs[2,1]),
            ("AT", "TA", "#121b37", "#2f468c", axs[1,0]), 
            ("CA", "GT", "#7f89aa", "#becdff", axs[1,1]), 
            ("AC", "TG", "#34543c", "#69a978", axs[0,0]), 
            ("CG", "GC", "#984533", "#ed6b50", axs[2,0])
        ] 
        for pair in pairs:
            fig.tight_layout(pad=1)
            height = config.figDistributionHeight
            fig.set_figheight(height)
            width = config.figDistributionWidth
            fig.set_figwidth(width)
            numberOfYTicks = config.figDistributionTicksY
            mut1 = pair[0]
            mut2 = pair[1]
            color1 = pair[2]
            color2 = pair[3]
            ax = pair[4]
            barDf = pd.concat([pd.Series(counts[mut1], name=mut1), pd.Series(counts[mut2], name=mut2)], axis=1)
            sns.barplot(y=barDf[mut1].values, x=barDf.index, color=color1, alpha=0.4, label=mut1, ax=ax)
            sns.barplot(y=barDf[mut2].values, x=barDf.index, color=color2, alpha=0.4, label=mut2, ax=ax)
            leg1 = mlines.Line2D([], [], color=color1, marker='s', ls='', label=mut1)
            leg2 = mlines.Line2D([], [], color=color2, marker='s', ls='', label=mut2)
            ax.legend(bbox_to_anchor=(0, 1.02, 3.6/width, 0.2), loc="lower left", mode="expand", ncol=2, 
                handles=[leg1, leg2])
            ax.xaxis.set_major_locator(mticker.FixedLocator(ax.get_xticks()))
            ax.tick_params(axis="x", labelrotation=90)
            ax.set_xticklabels([str(int(x * binSize)) for x in ax.get_xticks()])
            ax.set_xticklabels([str(int(x * binSize) + 1) + "-" + str(int((x + 1) * binSize)) for x in ax.get_xticks()])
            ax.yaxis.set_major_locator(mticker.MaxNLocator(numberOfYTicks))
            ax.set_yticklabels([str(abs(y)) for y in ax.get_yticks()])
            ax.set_ylabel('SNV count')

        plt.savefig(f"{resultsDir}/{ref}/graphs/{chrom}.histogram.png", dpi=config.figDPI)

def _graphRegression(df, ref, source=""):
    """It plots a regression graph that compares the total reads
    and the number of SNVs.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        source (str, optional): _description_. Defaults to "".
    """
    height = config.figRegressionHeight
    width = config.figRegressionWidth
    for chrom in df.CHROM.unique():
        fig, ax = plt.subplots(figsize=(width, height))
        sns.regplot(ax=ax, x="TotalReads", y="AltReads", data=df[df.CHROM == chrom][["TotalReads", "AltReads"]], scatter_kws={"s": 100})
        ax.set_xlabel('Total reads', fontsize=25)
        ax.set_ylabel('Alternative reads', fontsize=25)

    plt.savefig(f"{resultsDir}/{ref}/graphs/{chrom}.regression.png", dpi=config.figDPI)

def _aggregateByPositionPerMutation(group):
    genes = ', '.join(group['GeneName'].unique())
    gene_starts = group['GenePos'].min(), group['GenePos'].max()
    result = pd.Series({'GeneName': genes, 'GenePosMin': gene_starts[0], 'GenePosMax': gene_starts[1]})
    result = pd.concat([result, group.iloc[0, 1:1]], axis=0)
    result = pd.concat([result, group.iloc[0, 3:-1]], axis=0)
    return result

def _aggregateByPositionPerRun(group):
    genes = ', '.join(group['GeneName'].unique())
    gene_starts = group['GenePos'].min(), group['GenePos'].max()
    result = pd.Series({'GeneName': genes, 'GenePosMin': gene_starts[0], 'GenePosMax': gene_starts[1]})
    return result

def _graphCircosFrequencyPerMutation_plotGroup(df, chrom, group, mutations, maxFreq, mutationFontSize, colorBarFontSize, colorBarTickSize, pltSize, color, geneFontSize, titleFontSize, blankDegrees, centerSize, positionFontSize, cm, ref, run):
    fig, ax = plt.subplots()
    ax.axis("equal")
    for i in range(len(mutations)):
        dfTemp = []    
        mutation = mutations[i]
        dfTemp.append(df.loc[(df["Mutation"] == mutation) & (df["Group"] == group) & (df["CHROM"] == chrom), :])
        dfTemp = pd.concat(dfTemp)
        dfTemp = dfTemp.drop_duplicates().reset_index()
        dfTemp = dfTemp.groupby(["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax", "Sample"]).sum().reset_index()
        dfTemp["Degree"] = (360 - blankDegrees) / len(dfTemp)
        dfTemp = pd.concat([dfTemp, 
            pd.DataFrame([{"CHROM": "", "Position": "", "GeneName": "", 
                            "GenePosMin": (df["Position"].max() + 1), 
                            "GenePosMax": (df["Position"].max() + 1), 
                            "Sample": "", "Mutation": "",
                            "Frequency": 0.0, "Degree": blankDegrees}])], 
            ignore_index=True)
        dfTemp["Intensity"] = (((dfTemp["Frequency"]/maxFreq)**0.5)*255).astype(int).apply(hex)
        dfTemp["Intensity"] = dfTemp["Intensity"].str.partition("x")[2]
        dfTemp["Color"] = color
        dfTemp["Color"] = dfTemp["Color"].mask(dfTemp["Intensity"].str.len() == 2, color + dfTemp["Intensity"])
        dfTemp["Color"] = dfTemp["Color"].mask(dfTemp["Intensity"].str.len() != 2, color + "0" + dfTemp["Intensity"])
        shift = i // 3 * 1/3
        radius = pltSize - (pltSize - centerSize)/len(mutations) * i - shift
        ax.text(-0.2, radius - (pltSize - centerSize)/len(mutations)/2, mutations[i], style='italic', ha='right', va='center', fontsize=mutationFontSize)
        dfTemp = dfTemp.drop(columns=["Mutation", "Intensity"])
        wedges, texts = ax.pie(list(dfTemp["Degree"]),               
                            radius=radius, 
                            counterclock=False, startangle=-270,
                colors=list(dfTemp["Color"]),
                wedgeprops=dict(width=(pltSize - centerSize)/len(mutations), edgecolor="#000000FF", linewidth=5)
            )
        wedges[-1].set_visible(False)

    wedges, texts = ax.pie(list(dfTemp["Degree"]),
        radius=pltSize + 2, counterclock=False, startangle=-270,
        labels=list(dfTemp[["Position", "GeneName", "Sample"]].drop_duplicates()["Position"]),
        labeldistance=0.98, rotatelabels=True,
        colors=["#FFFFFF00"],
        wedgeprops=dict(width=(len(mutations)-1), alpha=0.0, linewidth=2),
        textprops=dict(rotation_mode="anchor", ha="center", va="center", fontsize=positionFontSize)
    )

    genes = (dfTemp[["GeneName", "GenePosMin", "GenePosMax", "Degree"]]
                .groupby(["GeneName", "GenePosMin", "GenePosMax"]).sum()
                .sort_values(["GenePosMin", "GenePosMax"]).reset_index())
    wedges, texts = ax.pie(list(genes["Degree"]),
            radius= centerSize - 1.2, startangle=-270, counterclock=False,
            labels=list(genes["GeneName"]),
            labeldistance=0.95, rotatelabels=True,
            colors=["#AAAAAA", "#333333"],
            wedgeprops=dict(width=0.15, edgecolor="white", linewidth=5),
            textprops=dict(rotation_mode="anchor", va="center", fontsize=geneFontSize)
        )
    wedges[-1].set_visible(False)
    for t in texts:
        if t.get_ha() == "right":
            t.set_ha("left")
        elif t.get_ha() == "left":
            t.set_ha("right")

    ax.text(0, -0.3, f"{chrom}", style='italic', ha='center', va='center', fontsize=titleFontSize)
    ax.text(0, 0.3, f"GROUP {int(group)}", style='italic', ha='center', va='center', fontsize=titleFontSize)
    
    img = ax.imshow(np.array([np.arange(0, 100, 0.1)]), cmap=cm, extent=[0, 0, 0, 0])
    ax.set_axis_off()
    cax = ax.inset_axes([-0.25, 0, 1.5, 0.2])
    cbar = fig.colorbar(img, orientation="horizontal", cax=cax)
    cbar.set_label('Frequency (0 is white)', fontsize=colorBarFontSize)
    cbar.set_ticks(range(0, 120, 20))
    cbar.ax.tick_params(labelsize=colorBarTickSize)
    
    fileName = f"{chrom}_{int(group)}"
    fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/circos/{fileName}.png", bbox_inches='tight')
    plt.clf()

    with open(f"{resultsDir}/{ref}/{run}/graphs/circos/_groups.txt", "a") as f:
        for gene in genes["GeneName"].unique():
            if gene.strip() != "":
                f.write("\n" + chrom + "\t" + gene + "\t" + str(int(group)))

def _graphCircosFrequencyPerMutation_plotGeneGroup(df, chrom, gene, group, mutations, maxFreq, mutationFontSize, colorBarFontSize, colorBarTickSize, pltSize, color, titleFontSize, blankDegrees, centerSize, positionFontSize, cm, ref, run):
    fig, ax = plt.subplots()
    ax.axis("equal")
    for i in range(len(mutations)):
        dfTemp = []    
        mutation = mutations[i]
        dfTemp.append(df.loc[(df["Mutation"] == mutation) & (df["GeneName"] == gene) & (df["Group"] == group) & (df["CHROM"] == chrom), :])
        dfTemp = pd.concat(dfTemp)
        dfTemp = dfTemp.drop_duplicates().reset_index().drop(columns=["index"])
        dfTemp = dfTemp.groupby(["CHROM", "Position", "GeneName", "Sample"]).sum().reset_index()
        dfTemp["Degree"] = (360 - blankDegrees) / len(dfTemp)
        dfTemp = pd.concat([dfTemp, 
            pd.DataFrame([{"CHROM": "", "Position": "", "GeneName": "", "Sample": "", "Mutation": "",
                            "Frequency": 0.0, "Degree": blankDegrees}])], 
            ignore_index=True)
        dfTemp["Intensity"] = (((dfTemp["Frequency"]/maxFreq)**0.5)*255).astype(int).apply(hex)
        dfTemp["Intensity"] = dfTemp["Intensity"].str.partition("x")[2]
        dfTemp["Color"] = color
        dfTemp["Color"] = dfTemp["Color"].mask(dfTemp["Intensity"].str.len() == 2, color + dfTemp["Intensity"])
        dfTemp["Color"] = dfTemp["Color"].mask(dfTemp["Intensity"].str.len() != 2, color + "0" + dfTemp["Intensity"])
        shift = i // 3 * 1/3
        radius = pltSize - (pltSize - centerSize)/len(mutations) * i - shift
        ax.text(-0.2, radius - (pltSize - centerSize)/len(mutations)/2, mutations[i], style='italic', ha='right', va='center', fontsize=mutationFontSize)
        dfTemp = dfTemp.drop(columns=["Mutation", "Intensity"])
        wedges, texts = ax.pie(list(dfTemp["Degree"]),               
                            radius=radius, 
                            counterclock=False, startangle=-270,
                colors=list(dfTemp["Color"]),
                wedgeprops=dict(width=(pltSize - centerSize)/len(mutations), edgecolor="#000000FF", linewidth=5)
            )
        wedges[-1].set_visible(False)

    wedges, texts = ax.pie(list(dfTemp["Degree"]),
        radius=pltSize + 2, counterclock=False, startangle=-270,
        labels=list(dfTemp[["Position", "GeneName", "Sample"]].drop_duplicates()["Position"]),
        labeldistance=0.98, rotatelabels=True,
        colors=["#FFFFFF00"],
        wedgeprops=dict(width=(len(mutations)-1), alpha=0.0, linewidth=2),
        textprops=dict(rotation_mode="anchor", ha="center", va="center", fontsize=positionFontSize)
    )

    ax.text(0, -0.3, f"{chrom}", style='italic', ha='center', va='center', fontsize=titleFontSize)
    ax.text(0, 0.3, f"{gene} - GROUP {int(group)}", style='italic', ha='center', va='center', fontsize=titleFontSize)
    
    img = ax.imshow(np.array([np.arange(0, 100, 0.1)]), cmap=cm, extent=[0, 0, 0, 0])
    ax.set_axis_off()
    cax = ax.inset_axes([-0.25, 0, 1.5, 0.2])
    cbar = fig.colorbar(img, orientation="horizontal", cax=cax)
    cbar.set_label('Frequency (0 is white)', fontsize=colorBarFontSize)
    cbar.set_ticks(range(0, 120, 20))
    cbar.ax.tick_params(labelsize=colorBarTickSize)
    
    fileName = f"{gene}_{int(group)}"
    fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/circos/{fileName}.png", bbox_inches='tight')
    plt.clf()

def _graphCircosFrequencyPerMutation(df, ref, run):
    """It plots a circos with each position on the x axis and the
    mutations on the Y axis, where the pie is a heatmap, each cell
    is colored with a different intensity depending on the frequency
    of that mutation in that position.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        run (_type_): _description_.
        source (_type_): _description_
    """
    pltSize = config.figCircosSize
    centerSize = config.figCircosCenterSize
    titleFontSize = config.figCircosTitleSize
    mutationFontSize = config.figCircosMutationLabelSize
    positionFontSize = config.figCircosPositionLabelSize
    colorBarFontSize = config.figCircosColorBarLabelSize
    blankDegrees = config.figCircosMutationsBlankDegrees
    minIntensity = config.figPositionGraphMinIntensity
    colorBarTickSize = config.figCircosColorBarTickSize
    geneFontSize = config.figCircosGeneLabelSize
    color = config.figCircosMutationsColor
    groupSize = config.figPositionsPerCircos

    minFreq = 0
    maxFreq = 100
    minColor = color + hex(int(((2/100)**0.5)*255))[-2:]
    maxColor = color + hex(255)[-2:]
    cm = mcolors.LinearSegmentedColormap.from_list(np.arange(0, 100, 0.1), [minColor, maxColor], 
        N=len(np.arange(0, 100, 0.1)))
    nonZero = df.Frequency > 0
    df.loc[nonZero, 'Frequency'] = minIntensity + (df.loc[nonZero, 'Frequency'] - minFreq) / (maxFreq - minFreq) * (1 - minIntensity)
    mutations = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    for mutation in mutations:
        df[mutation] = 0.0
        df[mutation].mask(df["Mutation"] == mutation, df["Frequency"], inplace=True)
    columns = ["CHROM", "Position", "GeneName", "Sample"] + mutations
    df = df[columns].drop_duplicates()
    genePosDf = df[["GeneName", "Position"]]
    genePosDf = genePosDf.merge(genePosDf.groupby("GeneName").min()
                                .rename(columns={"Position": "GenePos"}), 
                                    left_on="GeneName", right_index=True)
    genePosDf = genePosDf[["GeneName", "GenePos"]].drop_duplicates()
    genePosDf = genePosDf.set_index("GeneName")
    df = df.groupby(["CHROM", "Position", "Sample", "GeneName"]).sum().reset_index()
    df = df.set_index("GeneName")
    df = df.merge(genePosDf, left_index=True, right_index=True)
    df = df.reset_index()
    df = (df.groupby(["CHROM", "Position"]).apply(_aggregateByPositionPerMutation)
          .reset_index())
    counts = (df[["CHROM", "GeneName", "Position", "GenePosMin", "GenePosMax"]]
              .groupby(["CHROM", "GenePosMin", "GenePosMax", "GeneName"]).count()
              .rename(columns={"Position": "Count"}))
    counts = counts.sort_index()
    countsLessThanX = counts[counts["Count"] < groupSize]
    countsLessThanX["Sum"] = countsLessThanX.groupby("CHROM").cumsum()
    countsLessThanX["Group"] = countsLessThanX.Sum.floordiv(groupSize) + 1
    countsLessThanX["Group"] = countsLessThanX["Group"].astype(int)
    df = df.drop(columns=["GenePosMin", "GenePosMax"])
    df = pd.melt(df, id_vars=["CHROM", "Position", "GeneName", "Sample"], 
                 var_name="Mutation", value_name="Frequency")
    df = df.set_index(["CHROM", "GeneName"])
    df = df.merge(countsLessThanX, left_index=True, right_index=True, how="outer")
    df = df.reset_index()
    dfPerGene = df[df.Count.isna()]
    df = df[~df.Count.isna()]
    df = df.set_index("GeneName")
    df = df.sort_values(["Group", "GenePosMin", "GenePosMax"])
    maxFreq = 1

    if len(df) > 0:
        filePath = f"{resultsDir}/{ref}/{run}/graphs/circos/_groups.txt"
        with open(filePath, "w") as f:
            f.write("Chromosome\tGeneName\tGroup")

    with multiprocessing.Pool(processes=max(1, numProcesses // 4)) as pool:
        pool.starmap(_graphCircosFrequencyPerMutation_plotGroup, [(df, chrom, group, mutations, maxFreq, mutationFontSize, colorBarFontSize, colorBarTickSize, pltSize, color, geneFontSize, titleFontSize, blankDegrees, centerSize, positionFontSize, cm, ref, run)
            for chrom in df.CHROM.unique()
            for group in df[df.CHROM == chrom].Group.unique()])

        if len(dfPerGene) > 0:
            df = dfPerGene.drop(columns=["Count", "Sum", "Group"])
            dfAux = df[["CHROM", "GeneName", "Position"]].drop_duplicates()
            dfAux["Count"] = 1
            dfAux["Sum"] = dfAux.groupby(["CHROM", "GeneName"])["Count"].cumsum()
            dfAux["Group"] = dfAux["Sum"].floordiv(groupSize) + 1
            dfAux["Group"] = dfAux["Group"].astype(int)
            dfAux = dfAux.set_index(["CHROM", "GeneName", "Position"])
            df = df.set_index(["CHROM", "GeneName", "Position"])
            df = df.merge(dfAux, left_index=True, right_index=True)
            df = df.reset_index()

        pool.starmap(_graphCircosFrequencyPerMutation_plotGeneGroup, [(df, chrom, gene, group, mutations, maxFreq, mutationFontSize, colorBarFontSize, colorBarTickSize, pltSize, color, titleFontSize, blankDegrees, centerSize, positionFontSize, cm, ref, run)
            for chrom in df.CHROM.unique()
            for gene in dfPerGene.GeneName.unique()
            for group in df[(df.GeneName == gene) & (df.CHROM == chrom)].Group.unique()])
                    
def _graphCircosPresencePerRun_plotGroup(df, chrom, group, samples, pltSize, color, sampleFontSize, geneFontSize, titleFontSize, blankDegrees, centerSize, positionFontSize, ref):
        fig, ax = plt.subplots()
        ax.axis("equal")
        for i in range(len(samples)):
            dfTemp = []    
            sample = samples[i]
            dfTemp.append(df.loc[(df["Sample"] == sample) & (df["Group"] == group) & (df["CHROM"] == chrom), :])
            dfTemp = pd.concat(dfTemp)
            dfTemp = dfTemp.drop_duplicates().reset_index().drop(columns=["index"])
            dfTemp = dfTemp.groupby(["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax", "Sample"]).sum().reset_index()
            dfTemp["Degree"] = (360 - blankDegrees) / len(dfTemp) # Generate the degrees in the pie chart for each position
            dfTemp = pd.concat([dfTemp, 
                pd.DataFrame([{"CHROM": "", "Position": "", "GeneName": "", 
                                "GenePosMin": (df["Position"].max() + 1), 
                                "GenePosMax": (df["Position"].max() + 1), "Sample": "", 
                                "Presence": 0, "Degree": blankDegrees}])], 
                ignore_index=True)
            dfTemp["Intensity"] = (((dfTemp["Presence"])**0.5)*255).astype(int).apply(hex)
            dfTemp["Intensity"] = dfTemp["Intensity"].str.partition("x")[2]
            dfTemp["Color"] = ""
            dfTemp["Color"] = dfTemp["Color"].mask(dfTemp["Intensity"].str.len() == 2, color + dfTemp["Intensity"])
            dfTemp["Color"] = dfTemp["Color"].mask(dfTemp["Intensity"].str.len() != 2, color + "0" + dfTemp["Intensity"])
            shift = i // 3 * 1/3
            radius = pltSize - (pltSize - centerSize)/len(samples) * i - shift
            ax.text(-0.2, radius - (pltSize - centerSize)/len(samples)/2, samples[i], style='italic', ha='right', va='center', fontsize=sampleFontSize)
            wedges, texts = ax.pie(list(dfTemp["Degree"]),               
                                radius=radius, 
                                counterclock=False, startangle=-270,
                    colors=list(dfTemp["Color"]),
                    wedgeprops=dict(width=(pltSize - centerSize)/len(samples), edgecolor="#000000FF", linewidth=5)
                )
            wedges[-1].set_visible(False)
            
        wedges, texts = ax.pie(list(dfTemp["Degree"]),
            radius=pltSize + 2, counterclock=False, startangle=-270,
            labels=list(dfTemp[["Position", "GeneName", "Sample"]].drop_duplicates()["Position"]),
            labeldistance=0.98, rotatelabels=True,
            colors=["#FFFFFF00"],
            wedgeprops=dict(width=(len(samples)-1), alpha=0.0, linewidth=2),
            textprops=dict(rotation_mode="anchor", ha="center", va="center", fontsize=positionFontSize)
        )

        genes = (dfTemp[["GeneName", "GenePosMin", "GenePosMax", "Degree"]]
                    .groupby(["GeneName", "GenePosMin", "GenePosMax"]).sum()
                    .sort_values(["GenePosMin", "GenePosMax"]).reset_index())
        wedges, texts = ax.pie(list(genes["Degree"]),
                radius= centerSize - 1.2, startangle=-270, counterclock=False,
                labels=list(genes["GeneName"]),
                labeldistance=0.95, rotatelabels=True,
                colors=["#AAAAAA", "#333333"],
                wedgeprops=dict(width=0.15, edgecolor="white", linewidth=5),
                textprops=dict(rotation_mode="anchor", va="center", fontsize=geneFontSize)
            )
        wedges[-1].set_visible(False)
        for t in texts:
            if t.get_ha() == "right":
                t.set_ha("left")
            elif t.get_ha() == "left":
                t.set_ha("right")

        ax.text(0, -0.3, f"{chrom}", style='italic', ha='center', va='center', fontsize=titleFontSize)
        ax.text(0, 0.3, f"GROUP {int(group)}", style='italic', ha='center', va='center', fontsize=titleFontSize)
    
        fileName = f"{chrom}_{int(group)}"
        fig.savefig(f"{resultsDir}/{ref}/graphs/circos/{fileName}.png", bbox_inches='tight')
        plt.clf()

        with open(f"{resultsDir}/{ref}/graphs/circos/_groups.txt", "a") as f:
            for gene in genes["GeneName"].unique():
                if gene.strip() != "":
                    f.write("\n" + chrom + "\t" + gene + "\t" + str(int(group)))

def _graphCircosPresencePerRun_plotGeneGroup(df, chrom, gene, group, samples, pltSize, blankDegrees, color, centerSize, sampleFontSize, positionFontSize, titleFontSize, ref):
        fig, ax = plt.subplots()
        ax.axis("equal")
        for i in range(len(samples)):
            dfTemp = []    
            sample = samples[i]
            dfTemp.append(df.loc[(df["Sample"] == sample) & (df["GeneName"] == gene) & (df["Group"] == group) & (df["CHROM"] == chrom), :])
            dfTemp = pd.concat(dfTemp)
            dfTemp = dfTemp.drop_duplicates().reset_index().drop(columns=["index"])
            dfTemp = dfTemp.groupby(["CHROM", "Position", "GeneName", "Sample"]).sum().reset_index()
            dfTemp["Degree"] = (360 - blankDegrees) / len(dfTemp) # Generate the degrees in the pie chart for each position
            dfTemp = pd.concat([dfTemp, 
                pd.DataFrame([{"CHROM": "", "Position": "", "GeneName": "", "Sample": "", "Presence": 0, "Degree": blankDegrees}])], 
                ignore_index=True)
            dfTemp["Intensity"] = (((dfTemp["Presence"])**0.5)*255).astype(int).apply(hex)
            dfTemp["Intensity"] = dfTemp["Intensity"].str.partition("x")[2]
            dfTemp["Color"] = ""
            dfTemp["Color"] = dfTemp["Color"].mask(dfTemp["Intensity"].str.len() == 2, color + dfTemp["Intensity"])
            dfTemp["Color"] = dfTemp["Color"].mask(dfTemp["Intensity"].str.len() != 2, color + "0" + dfTemp["Intensity"])
            shift = i // 3 * 1/3
            radius = pltSize - (pltSize - centerSize)/len(samples) * i - shift
            ax.text(-0.2, radius - (pltSize - centerSize)/len(samples)/2, samples[i], style='italic', ha='right', va='center', fontsize=sampleFontSize)
            wedges, texts = ax.pie(list(dfTemp["Degree"]),               
                                radius=radius, 
                                counterclock=False, startangle=-270,
                    colors=list(dfTemp["Color"]),
                    wedgeprops=dict(width=(pltSize - centerSize)/len(samples), edgecolor="#000000FF", linewidth=5)
                )
            wedges[-1].set_visible(False)
            
        wedges, texts = ax.pie(list(dfTemp["Degree"]),
            radius=pltSize + 2, counterclock=False, startangle=-270,
            labels=list(dfTemp[["Position", "GeneName", "Sample"]].drop_duplicates()["Position"]),
            labeldistance=0.98, rotatelabels=True,
            colors=["#FFFFFF00"],
            wedgeprops=dict(width=(len(samples)-1), alpha=0.0, linewidth=2),
            textprops=dict(rotation_mode="anchor", ha="center", va="center", fontsize=positionFontSize)
        )

        ax.text(0, -0.3, f"{chrom}", style='italic', ha='center', va='center', fontsize=titleFontSize)
        ax.text(0, 0.3, f"{gene} - GROUP {int(group)}", style='italic', ha='center', va='center', fontsize=titleFontSize)
        
        fileName = f"{gene}_{int(group)}"
        fig.savefig(f"{resultsDir}/{ref}/graphs/circos/{fileName}.png", bbox_inches='tight')
        plt.clf()

def _graphCircosPresencePerRun(df, ref):
    """It plots a circos with each position on the x axis and the
    samples on the Y axis, where the pie is a heatmap, each cell
    is colored if the position shows mutations for that sample.

    Args:
        df (_type_): _description_
        ref (_type_): _description_
        run (_type_): _description_.
        source (_type_): _description_
    """
    pltSize = config.figCircosSize
    centerSize = config.figCircosCenterSize
    titleFontSize = config.figCircosTitleSize
    sampleFontSize = config.figCircosSampleLabelSize
    positionFontSize = config.figCircosPositionLabelSize
    geneFontSize = config.figCircosGeneLabelSize
    blankDegrees = config.figCircosSamplesBlankDegrees
    color = config.figCircosSamplesColor
    
    minColor = color + hex(int(((2/100)**0.5)*255))[-2:]
    maxColor = color + hex(255)[-2:]
    cm = mcolors.LinearSegmentedColormap.from_list(np.arange(0, 100, 0.1), [minColor, maxColor], 
        N=len(np.arange(0, 100, 0.1)))
    groupSize = config.figPositionsPerCircos
    columns = ["CHROM", "Position", "GeneName", "Sample"] 
    df = df[columns].drop_duplicates()
    genePosDf = df[["CHROM", "GeneName", "Position"]]
    genePosDf = (genePosDf.merge(genePosDf.groupby(["CHROM", "GeneName"]).min()
                                .rename(columns={"Position": "GenePos"}), 
                                left_on=["CHROM", "GeneName"], right_index=True))
    genePosDf = genePosDf[["GeneName", "GenePos"]].drop_duplicates()
    genePosDf = genePosDf.set_index("GeneName")
    df = df.set_index(["CHROM", "GeneName"])
    df = df.merge(genePosDf, left_index=True, right_index=True)
    df = df.reset_index()
    dfAggregated = (df.groupby(["CHROM", "Position"])
        .apply(_aggregateByPositionPerRun)
        .reset_index())
    df = pd.merge(df.drop(columns=["GeneName", "GenePos"]), dfAggregated, on=["CHROM", 'Position'], how='left')
    df = df.drop_duplicates()
    auxDf = (df[["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax"]]
            .value_counts().loc[lambda x: x > 1]
            .reset_index(name="Count")
            .drop(columns=["Count"]))
    df = df.set_index(["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax"])
    auxDf = auxDf.set_index(["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax"])
    df = df.merge(auxDf, left_index=True, right_index=True)
    df = df.reset_index()
    samples = list(df["Sample"].drop_duplicates())
    for sample in samples:
        df[sample] = 0
        df.loc[df["Sample"] == sample, sample] = 1
    df = (df.drop(columns="Sample").drop_duplicates()
        .groupby(["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax"])
        .sum().reset_index())
    counts = (df[["CHROM", "GeneName", "GenePosMin", "GenePosMax", "Position"]]
            .groupby(["CHROM", "GenePosMin", "GenePosMax", "GeneName"]).count()
            .rename(columns={"Position": "Count"}))
    counts = counts.sort_index()
    countsLessThanX = counts[counts["Count"] < groupSize]
    countsLessThanX["Sum"] = countsLessThanX.groupby("CHROM").cumsum()
    countsLessThanX["Group"] = countsLessThanX.Sum.floordiv(groupSize) + 1
    countsLessThanX["Group"] = countsLessThanX["Group"].astype(int)

    df = df.drop(columns=["GenePosMin", "GenePosMax"])
    df = pd.melt(df, id_vars=["CHROM", "Position", "GeneName"], var_name="Sample", 
                            value_name="Presence")
    df = df.set_index(["CHROM", "GeneName"])
    df = df.merge(countsLessThanX, left_index=True, right_index=True, how="outer")
    df = df.reset_index()
    dfPerGene = df[df.Count.isna()]
    df = df[~df.Count.isna()]
    df["Presence"] = np.minimum(df["Presence"], 1)
    df = df.sort_values(["Group", "GenePosMin", "GenePosMax"])

    samples = list(df["Sample"].unique()[:20]) #Limit to 20 rows

    if len(df) > 0:
        filePath = f"{resultsDir}/{ref}/graphs/circos/_groups.txt"
        with open(filePath, "w") as f:
            f.write("Chromosome\tGeneName\tGroup")

    with multiprocessing.Pool(processes=max(1, numProcesses // 4)) as pool:
        pool.starmap(_graphCircosPresencePerRun_plotGroup, [(df, chrom, group, samples, pltSize, color, sampleFontSize, geneFontSize, titleFontSize, blankDegrees, centerSize, positionFontSize, ref)
            for chrom in df.CHROM.unique()
            for group in df[df.CHROM == chrom].Group.unique()])

        if len(dfPerGene) > 0:
            df = dfPerGene.drop(columns=["Count", "Sum", "Group"])
            dfAux = df[["CHROM", "GeneName", "Position"]].drop_duplicates()
            dfAux["Count"] = 1
            dfAux["Sum"] = dfAux.groupby(["CHROM", "GeneName"])["Count"].cumsum()
            dfAux["Group"] = dfAux["Sum"].floordiv(groupSize) + 1
            dfAux["Group"] = dfAux["Group"].astype(int)
            dfAux = dfAux.set_index(["CHROM", "GeneName", "Position"])
            df = df.set_index(["CHROM", "GeneName", "Position"])
            df = df.merge(dfAux, left_index=True, right_index=True)
            df = df.reset_index()

        pool.starmap(_graphCircosPresencePerRun_plotGeneGroup, [(df, chrom, gene, group, samples, pltSize, blankDegrees, color, centerSize, sampleFontSize, positionFontSize, titleFontSize, ref)
            for chrom in df.CHROM.unique()
            for gene in dfPerGene.GeneName.unique()
            for group in df[(df.GeneName == gene) & (df.CHROM == chrom)].Group.unique()])

def _graphHeatmapFrequencyPerMutation_plotGroup(df, chrom, group, mutations, groupedMutations, pltSize, cmap, tickSize, geneSize, titleSize, titlePadding, minColor, maxColor, colorBarLabelSize, ref, run):
    dfGroup = df.loc[(df.CHROM == chrom) & (df.Group == group)].reset_index()
    dfPivot = dfGroup[["Position", "Mutation", "Frequency", "GeneName"]].pivot(index="Position", columns='Mutation', values='Frequency')
    dfPivot = dfPivot.reindex(columns=mutations, fill_value=0)
    dfPivot = dfPivot.filter(mutations)
    dfPivot = dfPivot.fillna(0)
    size = len(dfGroup.Position.unique())
    fig, axes = plt.subplots(1, len(groupedMutations), figsize=(pltSize, size * pltSize / 20))
    for i in range(len(groupedMutations)):
        groupDf = dfPivot.iloc[:, i * 3:(i + 1) * 3]
        ax = axes[i]
        data = groupDf.values
        cax = ax.pcolormesh(data, cmap=cmap, vmin=0, vmax=1, edgecolor="gainsboro", linewidths=1)
        for _, spine in ax.spines.items():
            spine.set_visible(True)
        ax.set_xticklabels(groupDf.columns)
        ax.set_xticks(np.arange(data.shape[1]) + 0.5)
        ax.set_yticklabels(groupDf.index)
        ax.set_yticks(np.arange(data.shape[0]) + 0.5)
        ax.tick_params(axis="x", bottom=False, top=False, labelbottom=True, labeltop=True, labelsize=tickSize)
        ax.tick_params(axis="y", left=False, labelleft=True, labelsize=tickSize)
        if i != 0:
            ax.set_ylabel("")
            ax.tick_params(axis="y", labelleft=False)
    geneGroups = dfPivot.merge(dfGroup[["Position", "GeneName", "GenePosMin", "GenePosMax"]].drop_duplicates(), right_on="Position", left_index=True)
    geneGroups = geneGroups.groupby(["GeneName", "GenePosMin", "GenePosMax"]).count()
    geneGroups = geneGroups.sort_values(["GenePosMin", "GenePosMax"])["Position"]
    rowCount = 0
    geneCount = 0
    colors = ["lightgray", "darkgray"]
    maxTextWidth = 0
    for i, row in geneGroups.items():
        color = colors[geneCount % 2]
        rect_x = 3.2 
        rect_y = rowCount
        rect_width = 0.25 
        rect_height = row
        rectangle = Rectangle((rect_x, rect_y), rect_width, rect_height, linewidth=1, edgecolor='none', facecolor=color, clip_on=False)
        ax.add_patch(rectangle)
        text_y = rect_y + (rect_height / 2)
        text_x = rect_x + rect_width + 0.25
        text = ax.text(text_x, text_y, i[0], ha='left', va='center', fontsize=geneSize)
        text_extent = text.get_window_extent(renderer=fig.canvas.get_renderer())
        maxTextWidth = max(maxTextWidth, text_extent.width)
        rowCount += row
        geneCount += 1
    common_ax = fig.add_subplot(111, frameon=False)
    top_ax = common_ax.twiny()
    top_ax.set_xticks([])
    top_ax.set_yticks([])
    top_ax.set_frame_on(False)
    bottom_ax = common_ax
    bottom_ax.set_xticks([])
    bottom_ax.set_yticks([])
    bottom_ax.set_frame_on(False)
    top_ax.set_xlabel("Mutation", labelpad=tickSize + 10, fontsize=tickSize)
    bottom_ax.set_xlabel("Mutation", labelpad=tickSize + 10, fontsize=tickSize)
    plt.title(f"{chrom} GROUP {int(group)}", fontsize=titleSize, pad=titlePadding)
    
    top_label_x_pixels = top_ax.xaxis.label.get_transform().transform((0.5, 0))[0]
    top_label_x_pixels = top_label_x_pixels
    colorBarXPosition = top_label_x_pixels / fig.dpi / fig.get_figwidth() - 0.3
    
    row_height = 1 / size
    colorbar_height = row_height / 2

    title_bbox = plt.gca().title.get_tightbbox(fig.canvas.get_renderer())
    title_bbox_fig = title_bbox.transformed(fig.transFigure.inverted())
    colorBarYPosition = title_bbox_fig.y0 - colorbar_height * 1.1

    cm = mcolors.LinearSegmentedColormap.from_list(np.arange(0, 100, 0.1), [minColor, maxColor], N=len(np.arange(0, 100, 0.1)))
    img = common_ax.imshow(np.array([np.arange(0, 100, 0.1)]), cmap=cm, extent=[0, 0, 0, 0])
    cax = fig.add_axes([colorBarXPosition, colorBarYPosition, 0.6, colorbar_height]) 
    cbar = plt.colorbar(img, cax=cax, orientation="horizontal")
    cbar.set_label('Frequency (0 is white)', fontsize=colorBarLabelSize)
    cbar.set_ticks(range(0, 120, 20))
    fig.subplots_adjust(wspace=0.1)
    fileName = f"{chrom}_{int(group)}"
    fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/heatmap/{fileName}.png", bbox_inches='tight')
    plt.clf()
    with open(f"{resultsDir}/{ref}/{run}/graphs/heatmap/_groups.txt", "a") as f:
        for i, row in geneGroups.items():
            f.write("\n" + chrom + "\t" + i[0] + "\t" + str(int(group)))

def _graphHeatmapFrequencyPerMutation_plotGeneGroup(df, chrom, gene, group, groupedMutations, pltSize, cmap, tickSize, titleSize, titlePadding, minColor, maxColor, colorBarLabelSize, ref, run):
    dfGroup = df.loc[(df.CHROM == chrom) & (df.Group == group) & (df.GeneName == gene)]
    dfPivot = dfGroup.pivot(index='Position', columns='Mutation', values='Frequency')
    dfPivot = dfPivot.filter(['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG'])
    dfPivot = dfPivot.fillna(0)
    size = len(dfGroup.Position.unique())
    fig, axes = plt.subplots(1, len(groupedMutations), figsize=(pltSize, size * pltSize / 20))
    for i in range(len(groupedMutations)):
        groupDf = dfPivot.iloc[:, i * 3:(i + 1) * 3]
        ax = axes[i]
        data = groupDf.values
        cax = ax.pcolormesh(data, cmap=cmap, vmin=0, vmax=1, edgecolor="gainsboro", linewidths=1)
        for _, spine in ax.spines.items():
            spine.set_visible(True)
        ax.set_xticklabels(groupDf.columns)
        ax.set_xticks(np.arange(data.shape[1]) + 0.5)
        ax.set_yticklabels(groupDf.index)
        ax.set_yticks(np.arange(data.shape[0]) + 0.5)
        ax.tick_params(axis="x", bottom=False, top=False, labelbottom=True, labeltop=True)
        ax.tick_params(axis="y", left=False, labelleft=True)
        if i != 0:
            ax.set_ylabel("")
            ax.tick_params(axis="y", labelleft=False)
    
    common_ax = fig.add_subplot(111, frameon=False)
    top_ax = common_ax.twiny()
    top_ax.set_xticks([])
    top_ax.set_yticks([])
    top_ax.set_frame_on(False)
    bottom_ax = common_ax
    bottom_ax.set_xticks([])
    bottom_ax.set_yticks([])
    bottom_ax.set_frame_on(False)
    top_ax.set_xlabel("Mutation", labelpad=20, fontsize=tickSize)
    bottom_ax.set_xlabel("Mutation", labelpad=20, fontsize=tickSize)
    
    plt.title(f"{gene} GROUP {int(group)}", fontsize=titleSize, pad=titlePadding)
        
    top_label_x_pixels = top_ax.xaxis.label.get_transform().transform((0.5, 0))[0]
    top_label_x_pixels = top_label_x_pixels
    colorBarXPosition = top_label_x_pixels / fig.dpi / fig.get_figwidth() - 0.3
    
    row_height = 1 / size
    colorbar_height = row_height / 2

    title_bbox = plt.gca().title.get_tightbbox(fig.canvas.get_renderer())
    title_bbox_fig = title_bbox.transformed(fig.transFigure.inverted())
    colorBarYPosition = title_bbox_fig.y0 - colorbar_height * 1.3

    cm = mcolors.LinearSegmentedColormap.from_list(np.arange(0, 100, 0.1), [minColor, maxColor], N=len(np.arange(0, 100, 0.1)))
    img = common_ax.imshow(np.array([np.arange(0, 100, 0.1)]), cmap=cm, extent=[0, 0, 0, 0])
    cax = fig.add_axes([colorBarXPosition, colorBarYPosition, 0.6, colorbar_height]) 
    cbar = plt.colorbar(img, cax=cax, orientation="horizontal")
    cbar.set_label('Frequency (0 is white)', fontsize=colorBarLabelSize)
    cbar.set_ticks(range(0, 120, 20))

    fig.subplots_adjust(wspace=0.1)
    fileName = f"{gene}_{int(group)}"
    fig.savefig(f"{resultsDir}/{ref}/{run}/graphs/heatmap/{fileName}.png", bbox_inches='tight')
    plt.clf()

def _graphHeatmapFrequencyPerMutation(df, ref, run):
    minIntensity = config.figPositionGraphMinIntensity
    groupSize = config.figPositionsPerHeatmap
    pltSize = config.figHeatmapSize
    titlePadding = config.figHeatmapTitlePadding
    titleSize = config.figHeatmapTitleSize
    colorBarLabelSize = config.figHeatmapColorBarLabelSize
    tickSize = config.figHeatmapTickSize
    geneSize = config.figHeatmapGeneLabelSize
    color = config.figHeatmapMutationsColor

    cmap = sns.blend_palette(["white", "#D6806D"], as_cmap=True)
    minColor = color + hex(int(((2/100)**0.5)*255))[-2:]
    maxColor = color + hex(255)[-2:]

    nonZero = df.Frequency > 0
    minFreq = 0
    maxFreq = 1
    columns = ["CHROM", "Position", "GeneName", "Sample", "Mutation", "Frequency"]
    df = df[columns].drop_duplicates()
    df.loc[nonZero, 'Frequency'] = minIntensity + (df.loc[nonZero, 'Frequency'] - minFreq) / (maxFreq - minFreq) * (1 - minIntensity)
    mutations = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
    for mutation in mutations:
        df[mutation] = 0.0
        df[mutation].mask(df["Mutation"] == mutation, df["Frequency"], inplace=True)
    groupedMutations = [['AC', 'AG', 'AT'], ['CA', 'CG', 'CT'], ['GA', 'GC', 'GT'], ['TA', 'TC', 'TG']]
    columns = ["CHROM", "Position", "GeneName", "Sample"] + mutations
    df = df[columns].drop_duplicates()
    genePosDf = df[["GeneName", "Position"]]
    genePosDf = genePosDf.merge(genePosDf.groupby("GeneName").min()
                                .rename(columns={"Position": "GenePos"}), 
                                    left_on="GeneName", right_index=True)
    genePosDf = genePosDf[["GeneName", "GenePos"]].drop_duplicates()
    genePosDf = genePosDf.set_index("GeneName")
    df = df.groupby(["CHROM", "Position", "Sample", "GeneName"]).sum().reset_index()
    df = df.set_index("GeneName")
    df = df.merge(genePosDf, left_index=True, right_index=True)
    df = df.reset_index()
    df = (df.groupby(["CHROM", "Position"]).apply(_aggregateByPositionPerMutation)
          .reset_index())
    counts = (df[["CHROM", "GeneName", "Position", "GenePosMin", "GenePosMax"]]
              .groupby(["CHROM", "GenePosMin", "GenePosMax", "GeneName"]).count()
              .rename(columns={"Position": "Count"}))
    counts = counts.sort_index()
    countsLessThanX = counts[counts["Count"] < groupSize]
    countsLessThanX["Sum"] = countsLessThanX.groupby("CHROM").cumsum()
    countsLessThanX["Group"] = countsLessThanX.Sum.floordiv(groupSize) + 1
    countsLessThanX["Group"] = countsLessThanX["Group"].astype(int)
    df = df.drop(columns=["GenePosMin", "GenePosMax"])
    df = pd.melt(df, id_vars=["CHROM", "Position", "GeneName", "Sample"], 
                 var_name="Mutation", value_name="Frequency")
    df = df.set_index(["CHROM", "GeneName"])
    df = df.merge(countsLessThanX, left_index=True, right_index=True, how="outer")
    df = df.reset_index()
    dfPerGene = df[df.Count.isna()]
    df = df[~df.Count.isna()]
    df = df.set_index("GeneName")
    df = df.sort_values(["Group", "GenePosMin", "GenePosMax"])
    
    if len(df) > 0:
        filePath = f"{resultsDir}/{ref}/{run}/graphs/heatmap/_groups.txt"
        with open(filePath, "w") as f:
            f.write("Chromosome\tGeneName\tGroup")

    with multiprocessing.Pool(processes=max(1, numProcesses // 4)) as pool:
        pool.starmap(_graphHeatmapFrequencyPerMutation_plotGroup, [(df, chrom, group, mutations, groupedMutations, pltSize, cmap, tickSize, geneSize, titleSize, titlePadding, minColor, maxColor, colorBarLabelSize, ref, run)
            for chrom in df.CHROM.unique()
            for group in df[df.CHROM == chrom].Group.unique()])

        if len(dfPerGene) > 0:
            df = dfPerGene.drop(columns=["Count", "Sum", "Group"])
            dfAux = df[["CHROM", "GeneName", "Position"]].drop_duplicates()
            dfAux["Count"] = 1
            dfAux["Sum"] = dfAux.groupby(["CHROM", "GeneName"])["Count"].cumsum()
            dfAux["Group"] = dfAux["Sum"].floordiv(groupSize) + 1
            dfAux["Group"] = dfAux["Group"].astype(int)
            dfAux = dfAux.set_index(["CHROM", "GeneName", "Position"])
            df = df.set_index(["CHROM", "GeneName", "Position"])
            df = df.merge(dfAux, left_index=True, right_index=True)
            df = df.reset_index()
        
        pool.starmap(_graphHeatmapFrequencyPerMutation_plotGeneGroup, [(df, chrom, gene, group, groupedMutations, pltSize, cmap, tickSize, titleSize, titlePadding, minColor, maxColor, colorBarLabelSize, ref, run)
            for chrom in df.CHROM.unique()
            for gene in dfPerGene.GeneName.unique()
            for group in df[(df.GeneName == gene) & (df.CHROM == chrom)].Group.unique()])

def _graphHeatmapPresencePerRun_plotGroup(df, chrom, group, samples, pltSize, cmap, tickSize, geneSize, titleSize, titlePadding, minColor, maxColor, colorBarLabelSize, ref):
    dfGroup = []
    dfGroup.append(df.loc[(df["Group"] == group) & (df["CHROM"] == chrom), :])
    dfGroup = pd.concat(dfGroup)
    dfPivot = dfGroup[["Position", "Sample", "GeneName", "Presence"]].pivot(index="Position", columns='Sample', values='Presence')
    dfPivot = dfPivot.reindex(columns=samples, fill_value=0)
    dfPivot = dfPivot.filter(samples)
    dfPivot = dfPivot.fillna(0)
    size = len(dfGroup.Position.unique())
    blockSize = 5
    numBlocks = math.ceil(len(samples) / blockSize)
    fig, axes = plt.subplots(figsize=(pltSize, size * pltSize / (len(samples) * 2)))
    ratios = [blockSize] * numBlocks
    ratios[-1] = len(dfPivot.columns) % blockSize
    gs = gridspec.GridSpec(1, numBlocks, width_ratios=ratios) 
    for i in range(numBlocks):
        groupDf = dfPivot.iloc[:, i * blockSize:(i + 1) * blockSize].sort_values("Position", ascending=False)
        ax = plt.subplot(gs[i])
        data = groupDf.values
        cax = ax.pcolormesh(data, cmap=cmap, vmin=0, vmax=1, edgecolor="gainsboro", linewidths=1)
        for _, spine in ax.spines.items():
            spine.set_visible(True)
        ax.set_xticks(np.arange(len(groupDf.columns)) + 0.5)
        ax.set_xticklabels(groupDf.columns)
        ax.set_yticks(np.arange(len(groupDf)) + 0.5)
        ax.set_yticklabels(groupDf.index)
        ax.tick_params(axis="x", bottom=False, top=False, labelbottom=True, labeltop=True, labelsize=tickSize, rotation=90)
        ax.tick_params(axis="y", left=False, labelleft=True, labelsize=tickSize)
        if i != 0:
            ax.set_ylabel("")
            ax.tick_params(axis="y", labelleft=False)

    geneGroups = dfPivot.merge(dfGroup[["Position", "GeneName", "GenePosMin", "GenePosMax"]].drop_duplicates(), right_on="Position", left_index=True)
    geneGroups = geneGroups.groupby(["GeneName", "GenePosMin", "GenePosMax"]).count()
    geneGroups = geneGroups.sort_values(["GenePosMin", "GenePosMax"])["Position"]
    rowCount = 0
    geneCount = 0
    colors = ["lightgray", "darkgray"]
    maxTextWidth = 0
    for i, row in geneGroups.items():
        color = colors[geneCount % 2]
        rect_x = len(groupDf.columns) + 0.2 
        rect_y = rowCount
        rect_width = 0.25 
        rect_height = row
        rectangle = Rectangle((rect_x, rect_y), rect_width, rect_height, linewidth=1, edgecolor='none', facecolor=color, clip_on=False)
        ax.add_patch(rectangle)
        text_y = rect_y + (rect_height / 2)
        text_x = rect_x + rect_width + 0.25
        text = ax.text(text_x, text_y, i[0], ha='left', va='center', fontsize=geneSize)
        text_extent = text.get_window_extent(renderer=fig.canvas.get_renderer())
        maxTextWidth = max(maxTextWidth, text_extent.width)
        rowCount += row
        geneCount += 1

    common_ax = fig.add_subplot(111, frameon=False)
    top_ax = common_ax.twiny()
    top_ax.set_xticks([])
    top_ax.set_yticks([])
    top_ax.set_frame_on(False)
    bottom_ax = common_ax
    bottom_ax.set_xticks([])
    bottom_ax.set_yticks([])
    bottom_ax.set_frame_on(False)
    
    plt.title(f"{chrom} GROUP {int(group)}", fontsize=titleSize, pad=titlePadding)
        
    top_label_x_pixels = top_ax.xaxis.label.get_transform().transform((0.5, 0))[0]
    top_label_x_pixels = top_label_x_pixels
    colorBarXPosition = top_label_x_pixels / fig.dpi / fig.get_figwidth() - 0.3
    
    row_height = 1 / size
    colorbar_height = row_height / 2

    title_bbox = plt.gca().title.get_tightbbox(fig.canvas.get_renderer())
    title_bbox_fig = title_bbox.transformed(fig.transFigure.inverted())
    colorBarYPosition = title_bbox_fig.y0 - colorbar_height * 1.3

    cm = mcolors.LinearSegmentedColormap.from_list(np.arange(0, 100, 0.1), [minColor, maxColor], N=len(np.arange(0, 100, 0.1)))
    img = common_ax.imshow(np.array([np.arange(0, 100, 0.1)]), cmap=cm, extent=[0, 0, 0, 0])
    cax = fig.add_axes([colorBarXPosition, colorBarYPosition, 0.6, colorbar_height]) 
    cbar = plt.colorbar(img, cax=cax, orientation="horizontal")
    cbar.set_label('Frequency (0 is white)', fontsize=colorBarLabelSize)
    cbar.set_ticks(range(0, 120, 20))

    fig.subplots_adjust(wspace=0.1)
    fileName = f"{chrom}_{int(group)}"
    fig.savefig(f"{resultsDir}/{ref}/graphs/heatmap/{fileName}.png", bbox_inches='tight')
    plt.clf()

    with open(f"{resultsDir}/{ref}/graphs/heatmap/_groups.txt", "a") as f:
        for i, row in geneGroups.items():
            f.write("\n" + chrom + "\t" + i[0] + "\t" + str(int(group)))

def _graphHeatmapPresencePerRun_plotGeneGroup(df, chrom, gene, group, samples, pltSize, cmap, tickSize, titleSize, titlePadding, minColor, maxColor, colorBarLabelSize, ref):
    dfGroup = []
    dfGroup.append(df.loc[(df.CHROM == chrom) & (df.Group == group) & (df.GeneName == gene)])
    dfGroup = pd.concat(dfGroup)
    dfPivot = dfGroup[["Position", "Sample", "GeneName", "Presence"]].pivot(index="Position", columns='Sample', values='Presence')
    dfPivot = dfPivot.reindex(columns=samples, fill_value=0)
    dfPivot = dfPivot.filter(samples)
    dfPivot = dfPivot.fillna(0)
    size = len(dfGroup.Position.unique())
    blockSize = 5
    numBlocks = math.ceil(len(samples) / blockSize)
    fig, axes = plt.subplots(figsize=(pltSize, size * pltSize / (len(samples) * 2)))
    ratios = [blockSize] * numBlocks
    ratios[-1] = len(dfPivot.columns) % blockSize
    gs = gridspec.GridSpec(1, numBlocks, width_ratios=ratios) 
    for i in range(numBlocks):
        groupDf = dfPivot.iloc[:, i * blockSize:(i + 1) * blockSize].sort_values("Position", ascending=False)
        ax = plt.subplot(gs[i])
        data = groupDf.values
        cax = ax.pcolormesh(data, cmap=cmap, vmin=0, vmax=1, edgecolor="gainsboro", linewidths=1)
        for _, spine in ax.spines.items():
            spine.set_visible(True)
        ax.set_xticks(np.arange(len(groupDf.columns)) + 0.5)
        ax.set_xticklabels(groupDf.columns)
        ax.set_yticks(np.arange(len(groupDf)) + 0.5)
        ax.set_yticklabels(groupDf.index)
        ax.tick_params(axis="x", bottom=False, top=False, labelbottom=True, labeltop=True, labelsize=tickSize, rotation=90)
        ax.tick_params(axis="y", left=False, labelleft=True, labelsize=tickSize)
        if i != 0:
            ax.set_ylabel("")
            ax.tick_params(axis="y", labelleft=False)

    common_ax = fig.add_subplot(111, frameon=False)
    top_ax = common_ax.twiny()
    top_ax.set_xticks([])
    top_ax.set_yticks([])
    top_ax.set_frame_on(False)
    bottom_ax = common_ax
    bottom_ax.set_xticks([])
    bottom_ax.set_yticks([])
    bottom_ax.set_frame_on(False)
    
    plt.title(f"{gene} GROUP {int(group)}", fontsize=titleSize, pad=titlePadding)
        
    top_label_x_pixels = top_ax.xaxis.label.get_transform().transform((0.5, 0))[0]
    top_label_x_pixels = top_label_x_pixels
    colorBarXPosition = top_label_x_pixels / fig.dpi / fig.get_figwidth() - 0.3
    
    row_height = 1 / size
    colorbar_height = row_height / 2

    title_bbox = plt.gca().title.get_tightbbox(fig.canvas.get_renderer())
    title_bbox_fig = title_bbox.transformed(fig.transFigure.inverted())
    colorBarYPosition = title_bbox_fig.y0 - colorbar_height * 1.3

    cm = mcolors.LinearSegmentedColormap.from_list(np.arange(0, 100, 0.1), [minColor, maxColor], N=len(np.arange(0, 100, 0.1)))
    img = common_ax.imshow(np.array([np.arange(0, 100, 0.1)]), cmap=cm, extent=[0, 0, 0, 0])
    cax = fig.add_axes([colorBarXPosition, colorBarYPosition, 0.6, colorbar_height]) 
    cbar = plt.colorbar(img, cax=cax, orientation="horizontal")
    cbar.set_label('Frequency (0 is white)', fontsize=colorBarLabelSize)
    cbar.set_ticks(range(0, 120, 20))

    fig.subplots_adjust(wspace=0.1)
    fileName = f"{gene}_{int(group)}"
    fig.savefig(f"{resultsDir}/{ref}/graphs/heatmap/{fileName}.png", bbox_inches='tight')
    plt.clf()

def _graphHeatmapPresencePerRun(df, ref):
    pltSize = config.figHeatmapSize
    titlePadding = config.figHeatmapTitlePadding
    titleSize = config.figHeatmapTitleSize
    colorBarLabelSize = config.figHeatmapColorBarLabelSize
    tickSize = config.figHeatmapTickSize
    geneSize = config.figHeatmapGeneLabelSize
    color = config.figHeatmapSamplesColor

    minColor = color + hex(int(((2/100)**0.5)*255))[-2:]
    maxColor = color + hex(255)[-2:]
    cm = mcolors.LinearSegmentedColormap.from_list(np.arange(0, 100, 0.1), [minColor, maxColor], 
        N=len(np.arange(0, 100, 0.1)))
    groupSize = config.figPositionsPerHeatmap
    cmap = sns.blend_palette(["white", "#D6806D"], as_cmap=True)
    columns = ["CHROM", "Position", "GeneName", "Sample"] 
    df = df[columns].drop_duplicates()
    genePosDf = df[["CHROM", "GeneName", "Position"]]
    genePosDf = (genePosDf.merge(genePosDf.groupby(["CHROM", "GeneName"]).min()
                                .rename(columns={"Position": "GenePos"}), 
                                left_on=["CHROM", "GeneName"], right_index=True))
    genePosDf = genePosDf[["GeneName", "GenePos"]].drop_duplicates()
    genePosDf = genePosDf.set_index("GeneName")
    df = df.set_index(["CHROM", "GeneName"])
    df = df.merge(genePosDf, left_index=True, right_index=True)
    df = df.reset_index()
    dfAggregated = (df.groupby(["CHROM", "Position"])
        .apply(_aggregateByPositionPerRun)
        .reset_index())
    df = pd.merge(df.drop(columns=["GeneName", "GenePos"]), dfAggregated, on=["CHROM", 'Position'], how='left')
    df = df.drop_duplicates()
    auxDf = (df[["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax"]]
            .value_counts().loc[lambda x: x > 1]
            .reset_index(name="Count")
            .drop(columns=["Count"]))
    df = df.set_index(["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax"])
    auxDf = auxDf.set_index(["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax"])
    df = df.merge(auxDf, left_index=True, right_index=True)
    df = df.reset_index()
    samples = list(df["Sample"].drop_duplicates())
    for sample in samples:
        df[sample] = 0
        df.loc[df["Sample"] == sample, sample] = 1
    df = (df.drop(columns="Sample").drop_duplicates()
        .groupby(["CHROM", "Position", "GeneName", "GenePosMin", "GenePosMax"])
        .sum().reset_index())
    counts = (df[["CHROM", "GeneName", "GenePosMin", "GenePosMax", "Position"]]
            .groupby(["CHROM", "GenePosMin", "GenePosMax", "GeneName"]).count()
            .rename(columns={"Position": "Count"}))
    counts = counts.sort_index()
    countsLessThanX = counts[counts["Count"] < groupSize]
    countsLessThanX["Sum"] = countsLessThanX.groupby("CHROM").cumsum()
    countsLessThanX["Group"] = countsLessThanX.Sum.floordiv(groupSize) + 1
    countsLessThanX["Group"] = countsLessThanX["Group"].astype(int)

    df = df.drop(columns=["GenePosMin", "GenePosMax"])
    df = pd.melt(df, id_vars=["CHROM", "Position", "GeneName"], var_name="Sample", 
                            value_name="Presence")
    df = df.set_index(["CHROM", "GeneName"])
    df = df.merge(countsLessThanX, left_index=True, right_index=True, how="outer")
    df = df.reset_index()
    dfPerGene = df[df.Count.isna()]
    df = df[~df.Count.isna()]
    df["Presence"] = np.minimum(df["Presence"], 1)
    df = df.sort_values(["Group", "GenePosMin", "GenePosMax"])

    samples = list(df["Sample"].unique()) 

    filePath = f"{resultsDir}/{ref}/graphs/heatmap/_groups.txt"
    with open(filePath, "w") as f:
        f.write("Chromosome\tGeneName\tGroup")

    with multiprocessing.Pool(processes=max(1, numProcesses // 4)) as pool:
        pool.starmap(_graphHeatmapPresencePerRun_plotGroup, [(df, chrom, group, samples, pltSize, cmap, tickSize, geneSize, titleSize, titlePadding, minColor, maxColor, colorBarLabelSize, ref)
            for chrom in df.CHROM.unique()
            for group in df[df.CHROM == chrom].Group.unique()])

        if len(dfPerGene) > 0:
            df = dfPerGene.drop(columns=["Count", "Sum", "Group"])
            dfAux = df[["CHROM", "GeneName", "Position"]].drop_duplicates()
            dfAux["Count"] = 1
            dfAux["Sum"] = dfAux.groupby(["CHROM", "GeneName"])["Count"].cumsum()
            dfAux["Group"] = dfAux["Sum"].floordiv(groupSize) + 1
            dfAux["Group"] = dfAux["Group"].astype(int)
            dfAux = dfAux.set_index(["CHROM", "GeneName", "Position"])
            df = df.set_index(["CHROM", "GeneName", "Position"])
            df = df.merge(dfAux, left_index=True, right_index=True)
            df = df.reset_index()

        pool.starmap(_graphHeatmapPresencePerRun_plotGeneGroup, [(df, chrom, gene, group, samples, pltSize, cmap, tickSize, titleSize, titlePadding, minColor, maxColor, colorBarLabelSize, ref)
            for chrom in df.CHROM.unique()
            for gene in dfPerGene.GeneName.unique()
            for group in df[(df.GeneName == gene) & (df.CHROM == chrom)].Group.unique()])