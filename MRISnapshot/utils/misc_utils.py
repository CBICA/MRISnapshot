#!/usr/bin/env python

#### Check if the input file exists
def checkFile(fileName):
    if not os.path.exists(fileName):
        sys.exit("Input file (" + fileName + ") does not exist")

### Write  text to log file
def writeLog(fileName, msg):
    lfp = open(fileName, 'a')
    lfp.write(msg + '\n')
    lfp.close()
