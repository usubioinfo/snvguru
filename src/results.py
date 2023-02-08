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
import glob

log = logger.logger

resultsDir = config.workPath + "/results"

def mergeCalling(sras):
    callingDir = config.workPath + "/calling"
    depthsDir = config.workPath + "/depths"
    util.makeDirectory(resultsDir)
    for fullRef in config.virusReferencePath:
        globalMain = None
        recReditools = None
        recJacusa = None
        globalRec = None
        ref = fullRef[fullRef.rfind("/") + 1:]
        ref = ref[:ref.rfind(".")]        
        util.makeDirectory(f"{resultsDir}/{ref}")
        for sra in sras:
            run = sra[0][0]
            run = run[:run.rfind(".")]
            type = sra[1]
            if type == "paired":
                run = run.replace("_1", "")
            util.makeDirectory(f"{resultsDir}/{ref}/{run}")
            files = glob.glob(f"{callingDir}/{ref}/{run}.calling.txt")
            if len(files) == 0:
                log.error(f"REDitools2 calling for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            log.info(f"Analyzing calling results from {run} vs {ref}...")
            reditools = pd.read_csv(f"{callingDir}/{ref}/{run}.calling.txt", sep='\t', header=0)
            reditools = reditools[reditools["Frequency"] >= config.minFrequency]
            reditools["BaseCount[A,C,G,T]"] = reditools["BaseCount[A,C,G,T]"].str[1:-1].str.replace(" ", "")
            reditools[["A", "C", "G", "T"]] = reditools["BaseCount[A,C,G,T]"].str.split(",", expand=True).apply(pd.to_numeric)
            reditools["Depth"] = reditools["A"] + reditools["C"] + reditools["G"] + reditools["T"]
            recReditools = reditools.copy()
            recReditools["FilteredSubs"] = recReditools.apply(lambda x: _filterMutationCt(x, 1), axis=1) # Just for formatting
            reditools = reditools[reditools["Depth"] >= config.minSNVCoverage]
            reditools["FilteredSubs"] = reditools.apply(lambda x: _filterMutationCt(x, config.minMainReadSupport), axis=1)
            reditoolsMain = reditools[reditools["AllSubs"] != ""]
            reditoolsMain = reditools.drop(columns=["gBaseCount[A,C,G,T]", "gAllSubs", "gFrequency", "gCoverage-q30", "gMeanQ"])

            
            files = glob.glob(f"{callingDir}/{ref}/{run}.vcf")
            if len(files) == 0:
                log.error(f"JACUSA calling for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            jacusa = _readJacusaVcf(f"{callingDir}/{ref}/{run}.vcf")
            jacusa[["GT", "BC", "DP"]] = jacusa["INFO2"].str.split(":", expand=True)
            jacusa[["A", "C", "G", "T"]] = jacusa["BC"].str.split(",", expand=True).apply(pd.to_numeric)
            jacusa["DP"] = jacusa["DP"].apply(pd.to_numeric)
            jacusa["Frequency"] = jacusa.apply(_calculateFrequencyJacusa, axis=1)
            jacusa = jacusa[jacusa["Frequency"] >= config.minFrequency]
            recJacusa = jacusa.copy()
            recJacusa["ALT"] = recJacusa.apply(lambda x: _filterMutationCt(x, 1), axis=1) # Just for formatting
            recJacusa = recJacusa.rename(columns={"A": "jacA", "C": "jacC", "G": "jacG", "T": "jacT"})
            jacusa = jacusa[jacusa["DP"] >= config.minSNVCoverage]
            jacusa["ALT"] = jacusa.apply(lambda x: _filterMutationCt(x, config.minMainReadSupport), axis=1)
            jacusa = jacusa.rename(columns={"A": "jacA", "C": "jacC", "G": "jacG", "T": "jacT"})
            jacusaMain = jacusa[jacusa["ALT"] != ""]

            files = glob.glob(f"{depthsDir}/{ref}/{run}_filtered.csv")
            if len(files) == 0:
                log.error(f"AS_StrandOddsRatio filter for run {run} against virus reference {ref} not found.")
                util.stopProgram()
            sor = pd.read_csv(f"{depthsDir}/{ref}/{run}_filtered.csv")
            reditoolsMain = pd.merge(reditoolsMain, sor["Position"], left_on="Position", right_on="Position", how="inner")
            jacusaMain = pd.merge(jacusaMain, sor["Position"], left_on="Position", right_on="Position", how="inner")
            
            main = pd.merge(reditoolsMain, jacusaMain[["Position", "ALT", "jacA", "jacC", "jacG", "jacT"]], left_on="Position", right_on="Position", how="inner")
            rec = pd.merge(recReditools, recJacusa[["Position", "ALT", "jacA", "jacC", "jacG", "jacT"]], left_on="Position", right_on="Position", how="inner")
            main["MergedSubs"] = main.apply(_checkSameMutations, axis=1)
            rec["MergedSubs"] = rec.apply(_checkSameMutations, axis=1)
            main = main.rename(columns={"ALT": "jacSubs"})
            rec = rec.rename(columns={"ALT": "jacSubs"})
            main = main[main["MergedSubs"] != ""]
            rec = rec[rec["MergedSubs"] != ""]
            main["RUN"] = run
            rec["RUN"] = run

            if globalMain is None:
                globalMain = main.copy()
            else:
                globalMain = pd.concat([globalMain, main])

            if globalRec is None:
                globalRec = rec.copy()
            else:
                globalRec = pd.concat([globalRec, rec])

            jacusaMain.to_excel(f"{resultsDir}/{ref}/{run}/jacusa.xlsx")
            reditoolsMain.to_excel(f"{resultsDir}/{ref}/{run}/reditools.xlsx")
            main.to_excel(f"{resultsDir}/{ref}/{run}/runCommon.xlsx")
            rec.to_excel(f"{resultsDir}/{ref}/{run}/runRecurrent.xlsx")

            main = main.assign(Sub=main['MergedSubs'].str.split(',')).explode('Sub')
            main["Mutation"] = main["Reference"] + main["Sub"]

            rec = rec.assign(Sub=rec['MergedSubs'].str.split(',')).explode('Sub')
            rec["Mutation"] = rec["Reference"] + rec["Sub"]

            _graphCounts(main, ref, run)
            _graphCounts(rec, ref, run, main=False)

        globalMain.to_excel(f"{resultsDir}/{ref}/globalCommon.xlsx", index=False)

        globalRec = globalRec.groupby("Position").filter(lambda x: len(x) > 1)
        globalRec = globalRec.groupby(["Position"]).apply(_checkSameMutationsGroup)
        globalRec = globalRec[globalRec["RecSubs"] != ""]
        globalRec.to_excel(f"{resultsDir}/{ref}/globalRecurrent.xlsx", index=False)

        globalRec = globalRec[globalRec["Frequency"] < 0.2]
        globalRec = globalRec.assign(Sub=globalRec['MergedSubs'].str.split(',')).explode('Sub')
        globalRec["Mutation"] = globalRec["Reference"] + globalRec["Sub"]
        globalRec = globalRec[["Depth", "Mutation"]].drop_duplicates()

        globalMain = globalMain[globalMain["Frequency"] < 0.2]
        globalMain = globalMain.assign(Sub=globalMain['MergedSubs'].str.split(',')).explode('Sub')
        globalMain["Mutation"] = globalMain["Reference"] + globalMain["Sub"]

        _graphCounts(globalMain, ref)
        _graphCounts(globalRec, ref, main=False)
        _graphFrequency(globalMain, ref)
        _graphHistograms(globalMain, ref)

def _checkSameMutations(row):
    reditools = row["FilteredSubs"].split(",")
    jacusa = row["ALT"].split(",")
    same = []
    for b in jacusa:
        if b in reditools:
            same.append(b)
    return ",".join(same)

def _checkSameMutationsGroup(rows):
    lists = rows["MergedSubs"].str.split(",")
    common = set(lists.iloc[0])
    for ls in lists:
        common = common.intersection(ls)
    common = ",".join(common)
    rows["RecSubs"] = common
    return rows

def _joinAllMutations(row):
    reditools = row["AllSubs"].split(",")
    jacusa = row["ALT"].split(",")
    bases = []
    if "A" in reditools or "A" in jacusa:
        bases.append("A")
    if "C" in reditools or "C" in jacusa:
        bases.append("C")
    if "G" in reditools or "G" in jacusa:
        bases.append("G")
    if "T" in reditools or "T" in jacusa:
        bases.append("T")
    return ",".join(bases)

def _filterMutationCt(row, min):
    bases = []
    if row["Reference"] != "A" and row["A"] >= min:
        bases.append("A")
    if row["Reference"] != "C" and row["C"] >= min:
        bases.append("C")
    if row["Reference"] != "G" and row["G"] >= min:
        bases.append("G")
    if row["Reference"] != "T" and row["T"] >= min:
        bases.append("T")
    return ",".join(bases)

def _calculateFrequencyJacusa(row):
    alts = row["ALT"].split(",")
    maxAlts = 0
    for alt in alts:
        ct = row[alt]
        maxAlts = max(ct, maxAlts)
    return maxAlts / (row[row["Reference"]] + maxAlts)

def _readJacusaVcf(path):
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
    ).rename(columns={'#CHROM': 'CHROM', "11": "INFO2", "REF": "Reference", "POS": "Position"})
    
    return csv

def _graphCounts(df, ref, run=None, main=True):
    counts = {"AC": 0, "AG": 0, "AT": 0, "CA": 0, "CG": 0, "CT": 0, "GA": 0, "GC": 0, "GT": 0, "TA": 0, "TC": 0, "TG": 0}
    for _, row in df.iterrows():
        counts[row["Mutation"]] += 1
    mean = df["Depth"].mean()
    mean = "{:.2f}".format(mean)
    median = df["Depth"].median()
    median = "{:.2f}".format(median)
    ax = sns.barplot(x=list(counts.keys()), y=list(counts.values()))
    if run is None:
        ax.set_title(f"Global SNV count per mutation for {ref}")
    else:
        ax.set_title(f"SNV count per mutation for {run} vs {ref}")
    ax.set_xlabel("Mutation")
    ax.set_ylabel("SNV count")
    ax.text(0.02, -0.2, f"Mean coverage: {mean}\nMedian coverage: {median}", ha="left", va="top", transform=ax.transAxes)
    fig = ax.get_figure()
    fig.tight_layout(pad=2)
    if not main:
        fileName = "recurrent.count"
    else:
        fileName = "common.count"
    if run is None:
        fig.savefig(f"{resultsDir}/{ref}.{fileName}.png") 
    else:
        fig.savefig(f"{resultsDir}/{ref}/{run}.{fileName}.png") 
    fig.clf()

def _graphFrequency(df, ref):
    df["Mutation"] = df["Reference"] + df["Sub"]
    df["Frequency"] *= 100
    ax = sns.stripplot(x="Mutation", y="Frequency", data=df, hue="Mutation", size=2,
        order = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"])
    ax.set_title(f"SNV frequency (%) for {ref}")
    ax.set_xlabel("Mutation")
    ax.set_ylabel("SNV frequency (%)")
    ax.legend_.remove()
    fig = ax.get_figure()
    fig.tight_layout(pad=2)
    fig.savefig(f"{resultsDir}/{ref}.frequency.png") 
    fig.clf()

def _graphHistograms(df, ref):
    maxDf = df["Position"].max()
    binSize = config.figDistributionBinSize
    maxBin = int(maxDf / binSize) 
    if maxDf % binSize != 0: 
        maxBin += 1
    keys = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
    counts = {key: {bin: 0 for bin in np.arange(0, maxBin, 1)} for key in keys}
    for index, row in df.iterrows():
        mutation = row["Mutation"]
        bin = int(row["Position"] / binSize)
        counts[mutation][bin] += 1
    
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
        fig.set_figheight(12)
        width = config.figDistributionWidth
        fig.set_figwidth(width)
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
        ax.xaxis.set_major_locator(plt.MaxNLocator(config.figDistributionTicksX))
        ax.xaxis.set_major_locator(mticker.FixedLocator(ax.get_xticks()))
        ax.set_xticklabels([str(int(x * binSize)) for x in ax.get_xticks()])
        yticks = np.arange(barDf[mut2].min(), barDf[mut1].max() + 1, 1)
        ax.yaxis.set_major_locator(mticker.FixedLocator(yticks))
        ax.set_yticklabels([str(int(abs(x))) for x in yticks])
        ax.set_ylabel('SNV count')

    plt.savefig(f"{resultsDir}/{ref}.distribution.png", dpi=config.figDPI)