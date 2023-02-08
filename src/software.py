import util
import distro
import logger
from fileinput import FileInput
import os

log = logger.logger

def downloadTools():
    softwareDir = "../tools"
    util.makeDirectory(softwareDir)

    print("~~~> Finding Linux distro name...")
    output = util.execCmd("hostnamectl")[0]
    if "CentOS" in output or "Rocky" in output:
        fileDistro = "centos_linux64"
        distro = "CentOS"
    elif "Debian" in output or "Ubuntu" in output or "Linux Mint" in output:
        fileDistro = "ubuntu64"
        distro = "Ubuntu"
    else:
        # Try to make
        print(f"~~~> There was an error trying to build SRA Toolkit.")
        util.stopProgram()
    print("~~~> Downloading SRA Toolkit for ...")
    util.execCmd(f"wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-{fileDistro}.tar.gz")
    util.execCmd(f"tar -xf sratoolkit.3.0.2-centos_linux64.tar.gz -C {softwareDir}")
    util.execCmd("rm sratoolkit.3.0.2-centos_linux64.tar.gz")
    print("~~~> Testing SRA Toolkit...")
    output = util.execCmd(f"{softwareDir}/sratoolkit.3.0.2-{fileDistro}/bin/fastq-dump --version")
    if output[0].split("\n")[1].strip() == f"{softwareDir}/sratoolkit.3.0.2-{fileDistro}/bin/fastq-dump : 3.0.2":
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


    print("~~~> Downloading FastQC...")
    util.execCmd("wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip")
    util.execCmd(f"unzip -o -q fastqc_v0.11.9.zip -d {softwareDir}")
    with FileInput(files=[f"{softwareDir}/FastQC/fastqc"], inplace=True, backup='.bak') as f:
        for line in f:
            if line == "my $java_bin = \"java\";\n":
                print("my $java_bin = \"$FindBin::Bin/../java-se-8u42-ri/bin/java\";")
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
    output = util.execCmd(f"{softwareDir}/java-se-8u42-ri/bin/java -jar {softwareDir}/Trimmomatic-0.39/trimmomatic-0.39.jar -version")
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

    print("~~~> Downloading BWA...")
    util.execCmd("wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2")
    util.execCmd(f"tar -xf bwa-0.7.17.tar.bz2 -C {softwareDir}")
    util.execCmd("rm bwa-0.7.17.tar.bz2")
    os.system(f"cd {softwareDir}/bwa-0.7.17 ; make" )
    print("~~~> Testing BWA...")
    output = util.execCmd(f"{softwareDir}/bwa-0.7.17/bwa")
    if output[1].split("\n")[2].strip() != "Version: 0.7.17-r1188":
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading STAR...")
    util.execCmd("wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.zip")
    util.execCmd(f"unzip -o -q 2.7.9a.zip -d {softwareDir}")
    util.execCmd("rm 2.7.9a.tar.gz")
    output = util.execCmd("grep avx /proc/cpuinfo")
    if output[0].strip() == "":
        cmd = f"cd {softwareDir}/STAR-2.7.9a/source ; make STAR CXXFLAGS_SIMD=sse --quiet"
    else: 
        cmd = f"cd {softwareDir}/STAR-2.7.9a/source ; make STAR --quiet"
    os.system(cmd)
    print("~~~> Testing STAR...")
    output = util.execCmd(f"{softwareDir}/STAR-2.7.9a/bin/Linux_x86_64_static/STAR")
    if output[0].split("\n")[3].strip() != "STAR version=2.7.9a":
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

    print("~~~> Downloading GMAP...")
    util.execCmd("wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2019-12-01.tar.gz")
    util.execCmd(f"tar -xf gmap-gsnap-2019-12-01.tar.gz -C {softwareDir}")
    util.execCmd("rm gmap-gsnap-2019-12-01.tar.gz")
    currPath = os.path.abspath(os.getcwd())
    os.system(f"cd {softwareDir}/gmap-2019-12-01 ; ./configure --prefix={currPath}/{softwareDir}/gmap ; make install -s ; make clean ")
    util.execCmd("rm -r gmap-2019-12-01")
    print("~~~> Testing GMAP...")
    output = util.execCmd(f"{softwareDir}/gmap/bin/gmap --version")
    if not output[1].startswith("GMAP version 2019-12-01"):
        print("~~~> Test failed.")
        print(output[1])
        util.stopProgram()
    else:
        print("~~~> Test successful.")

    print("~~~> Downloading samtools...")
    util.execCmd("wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2")
    util.execCmd(f"tar -xf samtools-1.9.tar.bz2 -C {softwareDir}")
    util.execCmd("rm samtools-1.9.tar.bz2")
    currPath = os.path.abspath(os.getcwd())
    os.system(f"cd {softwareDir}/samtools-1.9 ; ./configure --prefix={currPath}/{softwareDir}/samtools ; make ; make install -s ; make clean ")
    util.execCmd("rm -r samtools-1.9")
    print("~~~> Testing samtools...")
    output = util.execCmd(f"{softwareDir}/samtools/bin/samtools --version")
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
                print(f"{softwareDir}/java-se-8u42-ri/bin/java $java_options -classpath \"$QUALIMAP_HOME\"/qualimap.jar:\"$QUALIMAP_HOME\"/lib/* org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain \"${{ARGS[@]}}\"")
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



downloadTools()