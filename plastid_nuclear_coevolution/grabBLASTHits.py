'''Make sure that the BLAST result file, the database fasta file, and the 
BLASTHits.py file are all in the same directory. 

To run this in a unix environment move into the working directory containing 
your sequence list, the assembly file in fasta format, and the 
grabBLASTHits.py program. Then type the following:
    
    python grabBLASTHits.py <fastaFile> <hitFile> numHits

Where the <fastaFile> is the name of your database fasta file, the <hitFile> is 
the name of your BLAST result file, and numHits is the number of hits that you 
wish to extract from the database (if no number is included for this option, 
all hits will be extracted).

To run this in a python environment, start python, and type:
    
    import grabBLASTHits
    grabBLASTHits.makeFastaFile(<fastaFile>,<textFile>,numHits)
    
***Note: You will receive an error when importing this module into python. 
Disregard that error message, the module loaded properly.

This program is provided as is and will not be updated to reflect different 
input formats. If you would like some help getting it to work, or if you would 
like a program that accepts more flexible inputs, please contact Joel 
Sharbrough at jsharbro@rams.colostate.edu

Copyright 2016 All Rights Reserved
This program is not intended for commercial use and no profit may be made as a 
result of its use. Sharing, editing, and free use are heartily encouraged.'''
    
import sys

def makeSeqList(Infile,numHits=False):
    #Takes a text or xml BLAST output file as input and returns a list of contigs
    infile = open(Infile,'r')
    if Infile[-3:] == 'xml' or Infile[-3:] == 'XML':
        fileType = 'XML'
    else:
        fileType = 'TXT'
    seqList = []
    if fileType == 'TXT':
        seqName = ''
        for line in infile:
            if line[0] == '>':
                seqName = line
                while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                    seqName = seqName[0:-1]
                seqList.append(seqName)
            else:
                seqName = ''
    else:
        seqName = ''
        for line in infile:
            if 'Hit_def' in line:
                lineSplit = line.split('>')
                seqName = lineSplit[1]
                seqName = '>' + seqName[0:-9]
                #new starting from here
                #seqName = seqName.split('<')[0]
                #seqName = seqName.split(' ')[0]
            else:
                if seqName != '':
                    seqList.append(seqName)
                seqName = ''    
    infile.close()
    if type(numHits) == int:
        return seqList[0:numHits]
    else:
        return seqList

def makeFastaFile(fastaFile,hitFile,numHits=False):
    seqList = makeSeqList(hitFile)
    infile = open(fastaFile,'r') # Put the name of your input fasta file in this line
    fileName = hitFile.split('.')
    newFileName = ''
    for item in fileName[0:-1]:
        newFileName += item + '.'
    outfile = open(newFileName + 'fasta','w')  # Put the name of your output fasta file here
    seqDict = {}
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if currSeq != '':
                seqDict[seqName] = currSeq
            lineSplit = line.split(' ')
            seqName = lineSplit[0]
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\r' or currSeq[-1] == '\n' or currSeq[-1] == '\t':
                currSeq = currSeq[0:-1]
    seqDict[seqName] = currSeq
    if numHits == False:
        numHits = len(seqList)
    for item in seqList[0:numHits]:
        outfile.write(item + '\n')
        outfile.write(seqDict[item] + '\n')
    infile.close()
    outfile.close()

if len(sys.argv) == 4:
    makeFastaFile(sys.argv[1],sys.argv[2],int(sys.argv[3]))

else:
    makeFastaFile(sys.argv[1],sys.argv[2],)
        
