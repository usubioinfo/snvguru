"""This module contains utility functions.
"""

import os
import subprocess
import logger
import sys
import config
import re
from waiting import wait

log = logger.logger

def stopProgram():
    """Terminates the program.
    """
    handlers = log.handlers[:]
    for handler in handlers:
        log.removeHandler(handler)
        handler.close()
    sys.exit(0)

def makeDirectory(dir):
    """Creates a directory given a path.

    Args:
        dir (str): Path of the directory.
    """
    if not os.path.exists(dir):
        os.makedirs(dir)

def runCommand(command, jobName="", jobs=None, dep="", outFile=None):
    """Determines whether a command has to be run on SLURM or locally
    and executes the command using the appropriate function.

    It runs locally or on SLURM depending on the configuration in the
    main.config file, found in the config directory, or if the -hs 
    flag is active.

    Args:
        command (str): Command to be executed.
        jobName (str, optional): Name for the job in case it is executed on SLURM. Defaults to "".
        jobs (list, optional): List of jobs where the job ID will be appended in case it is executed on SLURM. Defaults to None.
        dep (str, optional): List of jobs that must be finished before this job is executed. The list must be a comma-separated string. Only applies for SLURM. Defaults to "".
        outFile (str, optional): Path to the output file. Defaults to None.

    Returns:
        str/tuple(str, str): Job ID if it is executed on SLURM or it is executed locally and the output is sent to a file. Standard output and error output if it is executed locally and no file is given.
    """
    if config.slurm:
        if outFile is not None:
            command += f" > {outFile}"
        job = runSlurm(jobName, command, dep=dep)
        if jobs != None:
            jobs.append(job)
        return job
    else:
        return execCmd(command, file=outFile)

def runSlurm(jobName, command, dep=""):
    """Runs the given command on SLURM.
    
    Args:
        jobName (str): Name for the job.
        command (str): Command to be executed.
        dep (str, optional): List of jobs that must be finished before this job is executed. The list must be a comma-separated string. Defaults to "".

    Returns:
        str: Job ID.
    """
    try:
        if dep != "":
            dep = f'--dependency=afterok:{dep} --kill-on-invalid-dep=yes'
        sbatchCmd = f"sbatch -J {jobName} -o {config.workPath}/logs/slurm/{jobName}-%j.out -e {config.workPath}/logs/slurm/{jobName}-%j.err -t {config.slurmTime}:00:00  --mem={config.slurmMem} --cpus-per-task={config.slurmCpus} --wrap='{command}' {dep}"
        output = subprocess.getoutput(sbatchCmd)
        jobId = output.split(' ')[-1].strip()
        log.info(f"===> Job {jobId}: {command}")
    except Exception as e:
        log.error(f"Job submission failed: {e}")
    return jobId

def _checkStatus(jobId):
    """Checks whether a job has been completed.

    Args:
        jobId (str): Job ID.

    Returns:
        bool: True if the job has been finished, False if not.
    """
    output = subprocess.check_output(f'squeue', shell=True, universal_newlines=True)
    if jobId in output:
        return False
    return True

def waitForJobs(jobs):
    """Waits until all the jobs in a given list are finished running.

    Args:
        jobs (list): List of jobs.
    """
    for job in jobs:
        wait(lambda: _checkStatus(job))

def execCmd(cmd, file=None, mode="w"):
    """Executes a command locally.

    Args:
        cmd (str): Command to be executed.
        file (str, optional): Path to the output file. Defaults to None.
        mode (str, optional): "w" if the file should be rewritten with the output, "a" if the output should be appended to an existing file. Defaults to "w".

    Returns:
        str/tuple(str, str): Empty string if the output is written to a file. Standard output and error output if no file is given.
    """
    if file != None:
        strOut = ""
        if mode == "w": 
            strOut = " > " + file
        elif mode == "a": 
            strOut = " >> " + file
        log.info(f"===> {cmd}{strOut}")
        with open(file, mode) as f:
            subprocess.run(cmd.split(), stdout=f)
            return ""
    else:
        log.info(f"===> {cmd}")
        output = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outString = output.stdout.decode("utf-8")
        errString = output.stderr.decode("utf-8")
        # print(outString)
        # print(errString)
        if "err" in outString or "Error:" in outString or "command not found" in outString or "No such file or directory" in outString:
            log.error(outString)
            stopProgram()
        elif "err" in errString or "Error:" in errString or "command not found" in errString or "No such file or directory" in errString:
            log.error(errString)
            stopProgram()
        return outString, errString