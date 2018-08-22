# add_codon_info_final.py
# final version to color only angiosperm branches based on RNA editing state

from Bio import SeqIO

#fasta file of codons
fastaFilename = "/home/williams/scratch_dir/RNAediting/angiosperm_seqs_only.trimmed.codon187ONLY.fas"
fastaFilename = "/home/williams/scratch_dir/RNAediting/angiosperm_seqs_only.trimmed.codon28ONLY.fas"
#make dictionary of keys (names) and values (codon state)
states = {}
#loop through and add key and value to dictionary
for record in SeqIO.parse(fastaFilename, "fasta"):
	states[record.id] = str(record.seq)

#read in tree as a big string
treeFile = open("/home/williams/scratch_dir/RNAediting/clpP1_output_tree.nex", "r")
treeString = treeFile.read()
treeFile.close()

#now find each dictionary key in the treeString and reformat its info 
# (so that it will include both a branch length and an RNA editing state)
for key in states:
	startLoc = treeString.find(key) #starting location of that name in treeString
	partOne = treeString[:startLoc+len(key)+1] #first part = genus_species + colon
	partTwo = treeString[startLoc+len(key)+1:] #second part = everything else
	#find first non-number after branch length
	commaLoc = partTwo.find(",",1)
	paraLoc = partTwo.find(")",1)
	colonLoc = partTwo.find(":",1)
	firstNonNum = min(commaLoc, paraLoc, colonLoc)
	partTwo_1 = partTwo[:firstNonNum]
	partTwo_2 = partTwo[firstNonNum:]
	#figure out which integer to encode after "&codon" (depends on codon at that site for that species)
	if states[key] in ("CAT", "CAC"): #if editing site has been retained
		inte = 1
	elif states[key] in ("TAT", "TAC"): #hard-coded
		inte = 2
	elif states[key].startswith("C"): #still a C, no RNA editing
		inte = 3
	else: #no RNA editing
		inte = 4
	#okay, now reformat.
	partTwo_1_new = partTwo_1 + "[&codon=" + str(inte) + "]"
	treeString = partOne + partTwo_1_new + partTwo_2 #reassign treeString
	
outputFile = open("/home/williams/scratch_dir/RNAediting/clpP1_output_tree_with&codon187.nex", "w")
outputFile = open("/home/williams/scratch_dir/RNAediting/clpP1_output_tree_with&codon28.nex", "w")
outputFile.write(treeString)
outputFile.close()




## for catalytic sites  ##
characterFile = open("catalytic_site_loss.txt", "r")
states = {}
#put catalytic site states in dictionary
for line in characterFile:
	line = line.strip()
	line = line.split("\t")
	states[line[0]] = line[1]

characterFile.close()

#read in tree as a big string
treeFile = open("/home/williams/scratch_dir/RNAediting/clpP1_output_tree.nex", "r")
treeString = treeFile.read()
treeFile.close()

#now find each dictionary key in the treeString and reformat its info 
# (so that it will include both a branch length and an RNA editing state)
for key in states:
	startLoc = treeString.find(key) #starting location of that name in treeString
	partOne = treeString[:startLoc+len(key)+1] #first part = genus_species + colon
	partTwo = treeString[startLoc+len(key)+1:] #second part = everything else
	#find first non-number after branch length
	commaLoc = partTwo.find(",",1)
	paraLoc = partTwo.find(")",1)
	colonLoc = partTwo.find(":",1)
	firstNonNum = min(commaLoc, paraLoc, colonLoc)
	partTwo_1 = partTwo[:firstNonNum]
	partTwo_2 = partTwo[firstNonNum:]
	#figure out which integer to encode after "&codon" (depends on codon at that site for that species)
	if states[key] == "his":
		inte = 1
	elif states[key] == "asp":
		inte = 2
	elif states[key] == "serandasp":
		inte = 3
	elif states[key] == "hisandasp":
		inte = 4
	elif states[key] == "all":
		inte = 5
	#okay, now reformat.
	partTwo_1_new = partTwo_1 + "[&state=" + str(inte) + "]"
	treeString = partOne + partTwo_1_new + partTwo_2 #reassign treeString
	
outputFile = open("/home/williams/scratch_dir/character_states/catalytic_site_clpp1_tree.nex", "w")
outputFile.write(treeString)
outputFile.close()



## for duplications ##

characterFile = open("duplications.txt","r")
states = {}
#put catalytic site states in dictionary
for line in characterFile:
	line = line.strip()
	line = line.split("\t")
	states[line[0]] = line[1]

characterFile.close()

#read in tree as a big string
treeFile = open("/home/williams/scratch_dir/RNAediting/clpP1_output_tree.nex", "r")
treeString = treeFile.read()
treeFile.close()

for key in states:
	startLoc = treeString.find(key) #starting location of that name in treeString
	partOne = treeString[:startLoc+len(key)+1] #first part = genus_species + colon
	partTwo = treeString[startLoc+len(key)+1:] #second part = everything else
	#find first non-number after branch length
	commaLoc = partTwo.find(",",1)
	paraLoc = partTwo.find(")",1)
	colonLoc = partTwo.find(":",1)
	firstNonNum = min(commaLoc, paraLoc, colonLoc)
	partTwo_1 = partTwo[:firstNonNum]
	partTwo_2 = partTwo[firstNonNum:]
	#okay, now reformat.
	partTwo_1_new = partTwo_1 + "[&state=" + states[key] + "]"
	treeString = partOne + partTwo_1_new + partTwo_2 #reassign treeString

outputFile = open("/home/williams/scratch_dir/character_states/duplications_clpp1_tree.nex", "w")
outputFile.write(treeString)
outputFile.close()



## for intron losses ##
characterFile = open("intron_losses.txt","r")
states = {}
#put catalytic site states in dictionary
for line in characterFile:
	print(line)
	line = line.strip()
	line = line.split("\t")
	states[line[0]] = line[1]

characterFile.close()

#read in tree as a big string
treeFile = open("/home/williams/scratch_dir/RNAediting/clpP1_output_tree.nex", "r")
treeString = treeFile.read()
treeFile.close()

for key in states:
	if key in treeString:
		startLoc = treeString.find(key) #starting location of that name in treeString
		partOne = treeString[:startLoc+len(key)+1] #first part = genus_species + colon
		partTwo = treeString[startLoc+len(key)+1:] #second part = everything else
		#find first non-number after branch length
		commaLoc = partTwo.find(",",1)
		paraLoc = partTwo.find(")",1)
		colonLoc = partTwo.find(":",1)
		firstNonNum = min(commaLoc, paraLoc, colonLoc)
		partTwo_1 = partTwo[:firstNonNum]
		partTwo_2 = partTwo[firstNonNum:]
		#figure out which integer to encode after "&codon" (depends on codon at that site for that species)
		if states[key] == "first":
			inte = 1
		elif states[key] == "second":
			inte = 2
		elif states[key] == "both":
			inte = 3
		#okay, now reformat.
		partTwo_1_new = partTwo_1 + "[&state=" + str(inte) + "]"
		treeString = partOne + partTwo_1_new + partTwo_2 #reassign treeString
	else: #I want to see which ones weren't in the tree string
		print(key)

outputFile = open("/home/williams/scratch_dir/character_states/intron_losses_land_plants_clpp1_tree2.nex", "w")
outputFile.write(treeString)
outputFile.close()























