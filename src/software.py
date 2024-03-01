"""This module downloads the required tools.
"""

import util
import logger
import config
from fileinput import FileInput

log = logger.logger

def checkTools():
    """Checks if all the needed tools are available.
    """
    tool = ""
    try:
        tool = "SRAtoolkit"
        util.execCmd(f"{config.sratoolkitPath}fastq-dump")
        tool = "Java"
        util.execCmd(f"{config.javaPath} -version")
        tool = "SAMtools"
        util.execCmd(f"{config.samtoolsPath}")
        tool = "Qualimap"
        util.execCmd(f"{config.qualimapPath} --version")
        tool = "FastQC"
        util.execCmd(f"{config.fastqcPath}")
        tool = "SnpEff"
        util.execCmd(f"{config.javaPath} -jar {config.snpEffPath} -version")
        tool = "BCFtools"
        util.execCmd(f"{config.bcftoolsPath}")

        if config.cropSoftware == "trimmomatic":
            tool = "Trimmomatic"
            util.execCmd(f"{config.javaPath} -jar {config.trimmomaticPath}")
        elif config.cropSoftware == "trimgalore":
            tool = "Trim Galore"
            util.execCmd(f"{config.trimGalorePath}")

        if config.callingSoftware in ["reditools", "both"]:
            tool = "REDitools2"
            util.execCmd(f"{config.reditoolsCommand}")
        elif config.callingSoftware in ["jacusa", "both"]:
            tool = "JACUSA 2"
            util.execCmd(f"{config.javaPath} -jar {config.jacusaPath}")

        if "bwa" in [config.alignmentSoftwareHost, config.alignmentSoftwarePathogen]:
            tool = "BWA"
            util.execCmd(f"{config.bwaPath}")
        elif "gmap" in [config.alignmentSoftwareHost, config.alignmentSoftwarePathogen]:
            tool = "GMAP"
            util.execCmd(f"{config.gmapPath}gmap --version")
        elif "hisat2" in [config.alignmentSoftwareHost, config.alignmentSoftwarePathogen]:
            tool = "Hisat2"
            util.execCmd(f"{config.hisat2Path}hisat2 --version")
        elif "magicblast" in [config.alignmentSoftwareHost, config.alignmentSoftwarePathogen]:
            tool = "Magic-BLAST"
            util.execCmd(f"{config.magicblastPath}magicblast -version")
        elif "minimap" in [config.alignmentSoftwareHost, config.alignmentSoftwarePathogen]:
            tool = "Minimap2"
            util.execCmd(f"{config.minimapPath} --version")
        elif "star" in [config.alignmentSoftwareHost, config.alignmentSoftwarePathogen]:
            tool = "STAR"
            util.execCmd(f"{config.starPath}")

    except Exception as e:
        print(f"{tool} not found. You can run the pipeline with the --download parameter to download the required software.")
        util.stopProgram()

def downloadTools():
    """Downloads and configures all the tools needed by the pipeline.

    It downloads the following tools:
        * Java 13
        * FastQC
        * Trimmomatic
        * Hisat2
        * Magic-BLAST
        * Minimap2
        * Qualimap
        * REDItools2
        * JACUSA 2
        * SnpEff
    """
    softwareDir = config.toolsPath
    util.makeDirectory(softwareDir)

    print("~~~> Testing SRA Toolkit...")
    output = util.execCmd(f"fastq-dump --version")
    if "fastq-dump : " in output[0].split("\n")[1].strip():
        print("~~~ Test successful.")
    else:
        print("~~~ Test failed.")
        print(output[0])
        util.stopProgram()

    print("~~~> Downloading Java 1.8...")
    util.execCmd("wget https://download.java.net/openjdk/jdk8u42/ri/openjdk-8u42-b03-linux-x64-14_jul_2022.tar.gz")
    util.execCmd(f"tar -xf openjdk-8u42-b03-linux-x64-14_jul_2022.tar.gz -C {softwareDir}")
    util.execCmd("rm openjdk-8u42-b03-linux-x64-14_jul_2022.tar.gz")
    print("~~~> Testing Java...")
    output = util.execCmd(f"{softwareDir}/java-se-8u42-ri/bin/java -version")
    if output[1].split("\n")[0].strip() == "openjdk version \"1.8.0_42\"":
        print("~~~ Test successful.")
    else:
        print("~~~ Test failed.")
        print(output[1])
        util.stopProgram()

    print("~~~> Downloading Java 13...")
    util.execCmd("wget https://download.java.net/java/GA/jdk13.0.2/d4173c853231432d94f001e99d882ca7/8/GPL/openjdk-13.0.2_linux-x64_bin.tar.gz")
    util.execCmd(f"tar -xf openjdk-13.0.2_linux-x64_bin.tar.gz -C {softwareDir}")
    util.execCmd("rm openjdk-13.0.2_linux-x64_bin.tar.gz")
    print("~~~> Testing Java...")
    output = util.execCmd(f"{softwareDir}/jdk-13.0.2/bin/java -version")
    if output[1].split("\n")[0].strip().startswith("openjdk version"):
        print("~~~ Test successful.")
    else:
        print("~~~ Test failed.")
        print(output[1])
        util.stopProgram()


    print("~~~> Downloading FastQC...")
    util.execCmd("wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip")
    util.execCmd(f"unzip -o -q fastqc_v0.11.9.zip -d {softwareDir}")
    with FileInput(files=[f"{softwareDir}/FastQC/fastqc"], inplace=True, backup='.bak') as f:
        for line in f:
            if line == "my $java_bin = \"java\";\n":
                print("my $java_bin = \"$FindBin::Bin/../jdk-13.0.2/bin/java\";")
            else:
                print(line, end="")
    util.execCmd(f"chmod +x {softwareDir}/FastQC/fastqc")
    util.execCmd("rm fastqc_v0.11.9.zip")
    print("~~~> Testing FastQC...")
    output = util.execCmd(f"{softwareDir}/FastQC/fastqc --version")
    if output[0].strip() != "FastQC v0.11.9":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading Trimmomatic...")
    util.execCmd("wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip")
    util.execCmd(f"unzip -o -q Trimmomatic-0.39.zip -d {softwareDir}")
    util.execCmd("rm Trimmomatic-0.39.zip")
    print("~~~> Testing Trimmomatic...")
    output = util.execCmd(f"{softwareDir}/jdk-13.0.2/bin/java -jar {softwareDir}/Trimmomatic-0.39/trimmomatic-0.39.jar -version")
    if output[0].strip() != "0.39":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading Hisat2...")
    util.execCmd("wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download")
    util.execCmd(f"unzip -o -q download -d {softwareDir}")
    util.execCmd("rm download")
    print("~~~> Testing Hisat2...")
    output = util.execCmd(f"{softwareDir}/hisat2-2.1.0/hisat2 --version")
    if "hisat2-2.1.0/hisat2-align-s version 2.1.0" not in output[0].split("\n")[0].strip():
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Testing BWA...")
    output = util.execCmd(f"bwa")
    if not output[1].split("\n")[2].strip().startswith("Version: "):
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Testing STAR...")
    output = util.execCmd(f"STAR")
    if not output[0].split("\n")[3].strip().startswith("STAR version="):
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading Magic-BLAST...")
    util.execCmd("wget https://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/1.6.0/ncbi-magicblast-1.6.0-x64-linux.tar.gz")
    util.execCmd(f"tar -xf ncbi-magicblast-1.6.0-x64-linux.tar.gz -C {softwareDir}")
    util.execCmd("rm ncbi-magicblast-1.6.0-x64-linux.tar.gz")
    print("~~~> Testing Magic-BLAST...")
    output = util.execCmd(f"{softwareDir}/ncbi-magicblast-1.6.0/bin/magicblast -version")
    if output[0].split("\n")[0].strip() != "magicblast: 1.6.0":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading Minimap2...")
    util.execCmd("wget https://github.com/lh3/minimap2/releases/download/v2.20/minimap2-2.20_x64-linux.tar.bz2")
    util.execCmd(f"tar -xf minimap2-2.20_x64-linux.tar.bz2 -C {softwareDir}")
    util.execCmd("rm minimap2-2.20_x64-linux.tar.bz2")
    print("~~~> Testing Minimap2...")
    output = util.execCmd(f"{softwareDir}/minimap2-2.20_x64-linux/minimap2 --version")
    if output[0].split("\n")[0].strip() != "2.20-r1061":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Testing GMAP...")
    output = util.execCmd(f"gmap --version")
    if "GMAP version 2020" not in output[1]:
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Testing samtools...")
    output = util.execCmd(f"samtools --version")
    if output[0].split("\n")[0].strip() != "samtools 1.9":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading Qualimap 2...")
    util.execCmd("wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip")
    util.execCmd(f"unzip -o -q qualimap_v2.2.1.zip -d {softwareDir}")
    util.execCmd("rm qualimap_v2.2.1.zip")
    with FileInput(files=[f"{softwareDir}/qualimap_v2.2.1/qualimap"], inplace=True, backup='.bak') as f:
        for line in f:
            if line == "java $java_options -classpath \"$QUALIMAP_HOME\"/qualimap.jar:\"$QUALIMAP_HOME\"/lib/* org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain \"${ARGS[@]}\"\n":
                print(f"{softwareDir}/jdk-13.0.2/bin/java $java_options -classpath \"$QUALIMAP_HOME\"/qualimap.jar:\"$QUALIMAP_HOME\"/lib/* org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain \"${{ARGS[@]}}\"")
            else:
                print(line, end="")
    print("~~~> Testing Qualimap 2...")
    output = util.execCmd(f"{softwareDir}/qualimap_v2.2.1/qualimap --version")
    if output[0].split("\n")[3].strip() != "QualiMap v.2.2.1":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading REDItools2...")
    util.execCmd(f"rm -r -f {softwareDir}/REDItools2")
    util.execCmd("git clone https://github.com/BioinfoUNIBA/REDItools2.git")
    util.execCmd(f"mv REDItools2 {softwareDir}/")
    util.makeDirectory(f"{softwareDir}/REDItools2/env")
    util.execCmd(f"virtualenv -p python2 {softwareDir}/REDItools2/env")
    util.execCmd(f"{softwareDir}/REDItools2/env/bin/pip install --upgrade pip")
    with FileInput(files=[f"{softwareDir}/REDItools2/requirements.txt"], inplace=True, backup='.bak') as f:
        for line in f:
            if line == "mpi4py\n":
                print("3to2")
            else:
                print(line, end="")
    util.execCmd(f"{softwareDir}/REDItools2/env/bin/pip install -r {softwareDir}/REDItools2/requirements.txt")
    util.execCmd(f"{softwareDir}/REDItools2/env/bin/pip install --no-use-pep517 mpi4py")
    print("~~~> Testing REDItools2...")
    output = util.execCmd(f"{softwareDir}/REDItools2/env/bin/python {softwareDir}/REDItools2/src/cineca/reditools.py")
    print(output)
    if output[0].strip() != "[ERROR] An input bam file is mandatory. Please, provide one (-f|--file)":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading JACUSA 2...")
    util.execCmd("wget https://github.com/dieterich-lab/JACUSA2/releases/download/v2.0.2-RC/JACUSA_v2.0.2-RC.jar")
    util.makeDirectory(f"{softwareDir}/jacusa")
    util.execCmd(f"mv JACUSA_v2.0.2-RC.jar {softwareDir}/jacusa/")
    print("~~~> Testing JACUSA 2...")
    output = util.execCmd(f"{softwareDir}/jdk-13.0.2/bin/java -jar {softwareDir}/jacusa/JACUSA_v2.0.2-RC.jar")
    if output[1].split("\n")[0].strip() != "":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading SnpEff...")
    util.execCmd("wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip")
    util.execCmd(f"unzip -o -q snpEff_latest_core.zip -d {softwareDir}")
    util.execCmd("rm snpEff_latest_core.zip")
    print("~~~> Testing SnpEff...")
    output = util.execCmd(f"{softwareDir}/jdk-13.0.2/bin/java -jar {softwareDir}/snpEff/snpEff.jar -version")
    if not output[0].strip().startswith("SnpEff"):
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")