import os
import subprocess
import logger
import sys

log = logger.logger

def stopProgram():
    handlers = log.handlers[:]
    for handler in handlers:
        log.removeHandler(handler)
        handler.close()
    sys.exit(0)

def makeDirectory(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def execCmd(cmd, file=None, mode="w"):
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
        if "err" in outString or "Error" in outString or "command not found" in outString:
            log.error(outString)
            stopProgram()
        elif "err" in errString or "Error" in errString or "command not found" in errString:
            log.error(errString)
            stopProgram()
        return outString, errString