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
        if config.callingSoftware in ["reditools", "both"]:
            tool = "REDitools2"
            util.execCmd(f"{config.reditoolsCommand}")

    except Exception as e:
        print(f"{tool} not found. You can run the pipeline with the --download parameter to download the required software.")
        util.stopProgram()

def downloadTools():
    """Downloads and configures all the tools needed by the pipeline.

    It downloads the following tools:
        * REDItools2
    """
    softwareDir = config.toolsPath
    util.makeDirectory(softwareDir)

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