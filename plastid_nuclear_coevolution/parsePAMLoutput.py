#parsePAMLoutput.py
#Function to pull out the nwk formatted trees from PAML txt output
#Alissa Williams
#March 14, 2017


import sys
import re

def main():

	pamlFilename = sys.argv[1] #name of paml output file
	outputTreeFilename = sys.argv[2] #what you want output file of just tree to be called
	
	pamlOutputFile = open(pamlFilename, "r")
	
	for line in pamlOutputFile:
		if line.startswith("("):
			tree = line #this should update twice, the second time being the tree with species names
			
	pamlOutputFile.close()
	
	#open outputTreeFilename and write tree to it
	outputTreeFile = open(outputTreeFilename, "w")
	outputTreeFile.write(tree)
	outputTreeFile.close()
	


if __name__ == "__main__":
	main()