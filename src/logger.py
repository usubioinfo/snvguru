import logging
import config
import os
from datetime import datetime

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s", datefmt='%H:%M:%S')

logger_blocklist = [
    "numexpr",
    "matplotlib",
]

for module in logger_blocklist:
    logging.getLogger(module).setLevel(logging.CRITICAL)

if config.logFile:
    dir = f"{config.workPath}/logs"
    if not os.path.exists(dir):
        os.makedirs(dir)
    file = f"{dir}/pipeline.log"
    if os.path.exists(file):
        prevCreatedDate = os.path.getctime(file)
        createdDateStr = datetime.fromtimestamp(prevCreatedDate).strftime("%Y-%m-%d_%H:%M")
        os.rename(file, f"{dir}/{createdDateStr}.log")
    fileHandler = logging.FileHandler(file)
    fileHandler.setFormatter(logFormatter)
    fileHandler.setLevel(logging.DEBUG)
    logger.addHandler(fileHandler)

if config.logConsole:
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    consoleHandler.setLevel(logging.DEBUG)
    logger.addHandler(consoleHandler)