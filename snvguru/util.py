"""This module contains utility functions.
"""

import os
import subprocess
from snvguru import logger
import sys
from snvguru import config
from Bio import SeqIO
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

def _get_fasta_ids(fasta_file):
    """Extracts chromosome IDs from the FASTA headers."""
    return {rec.id for rec in SeqIO.parse(fasta_file, "fasta")}

def _fix_genbank(fasta_ids, gbk_in, gbk_out):
    """
    Sync GenBank LOCUS names to match FASTA headers.

    Rules:
      * If FASTA headers include version numbers (e.g. 'NC_000001.11') and LOCUS does not,
        append version to LOCUS.
      * If FASTA headers omit version numbers (e.g. 'NC_000001') and LOCUS includes them,
        strip version from LOCUS.
    """
    # Detect whether FASTA IDs use version numbers
    has_version_in_fasta = any(re.search(r"\.\d+$", fid) for fid in fasta_ids)

    out_records = []
    for rec in SeqIO.parse(gbk_in, "genbank"):
        locus = rec.name or ""                      # LOCUS comes in as record.name
        locus_has_version = bool(re.search(r"\.\d+$", locus))

        if has_version_in_fasta and not locus_has_version:
            # prefer explicit accession + sequence_version from annotations
            acc = None
            if rec.annotations.get("accessions"):
                acc = rec.annotations["accessions"][0]
            seqver = rec.annotations.get("sequence_version")

            if acc and seqver is not None:
                new_locus = f"{acc}.{seqver}"
            elif re.search(r"\.\d+$", rec.id):
                # fallback: use record.id if it already contains a version
                new_locus = rec.id
            else:
                # nothing reliable to append — keep existing locus
                new_locus = locus

            if new_locus and new_locus != locus:
                rec.name = new_locus

        elif (not has_version_in_fasta) and locus_has_version:
            rec.name = re.sub(r"\.\d+$", "", locus)

        out_records.append(rec)

    SeqIO.write(out_records, gbk_out, "genbank")

def _fix_gff_or_gtf_fast(fasta_ids, ann_in, ann_out):
    """
    Rewrite GFF/GTF seqid column to match FASTA headers.
    fasta_ids: set of FASTA headers
    """
    # Build mapping: base ID (no version) -> canonical FASTA ID
    fasta_map = {}
    for fid in fasta_ids:
        base = fid.split(".")[0]
        fasta_map[base] = fid

    with open(ann_in) as f_in, open(ann_out, "w") as f_out:
        for line in f_in:
            if line.startswith("#"):
                f_out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                f_out.write(line)
                continue

            seqid = parts[0]
            base_seqid = seqid.split(".")[0]

            if base_seqid in fasta_map:
                parts[0] = fasta_map[base_seqid]
            else:
                # Optional: warn user about unmatched seqid
                print(f"Warning: seqid '{seqid}' not found in FASTA headers. Keeping original.")

            f_out.write("\t".join(parts) + "\n")

def sync_annotation_to_fasta(fasta_file, annotation_file, output_file):
    """
    Synchronize an annotation file with a reference FASTA so that sequence IDs match.

    Args:
        fasta_file (str): Path to the reference genome FASTA file. The chromosome IDs in this
                        file are treated as canonical.
        annotation_file (str): Path to the annotation file. Supported formats are GenBank (.gb, .gbk),
                            GFF/GFF3 (.gff, .gff3), or GTF (.gtf).
        output_file (str): Path where the corrected annotation file will be written.
    """
    fasta_ids = _get_fasta_ids(fasta_file)
    ext = os.path.splitext(annotation_file)[1].lower()

    if ext in [".gb", ".gbk"]:
        _fix_genbank(fasta_ids, annotation_file, output_file)
    elif ext in [".gff", ".gff3", ".gtf"]:
        _fix_gff_or_gtf(fasta_ids, annotation_file, output_file)
    elif ext == ".refSeq":
        execCmd(f"cp {annotation_file} {output_file}")
    else:
        raise ValueError(f"Unsupported annotation format: {annotation_file}")