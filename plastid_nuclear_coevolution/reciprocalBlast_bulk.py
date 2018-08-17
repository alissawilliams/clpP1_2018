#reciprocalBlast_bulk.py

import os
import sys
from Bio import SeqIO

def singleReciprocalBlast(query_seq, query_db, blast_output_file, reference_db, reciprocal_blast_file, compiled_results_file, summary_file):
				
	#we need to remember the sequence ID of the original fasta file
	querySeq = next(SeqIO.parse(query_seq, 'fasta'))

	#run blastp and save output to an xml file
	command = "blastp -query " + query_seq + " -db " + query_db + " -outfmt 5 -out " + blast_output_file
	os.system(command)

	#now pull out the first sequence based on the sequence IDs found by the blastp command (this line creates a fasta file with the first hit in it)
	command = "python grabBLASTHits.py " + query_db + " " + blast_output_file + " 1"
	os.system(command)

	#we need to use the output of grabBLASTHits.py, so we need a handle for that file 
	firstFastaFilename = blast_output_file
	firstFastaFilename = firstFastaFilename.replace(".xml", ".fasta")

	#now we can Blast that top result against Arabidopsis
	command = "blastp -query " + firstFastaFilename + " -db " + reference_db + " -outfmt 5 -out " + reciprocal_blast_file 
	os.system(command)

	#now pull out the first hit of the reciprocal Blast
	command = "python grabBLASTHits.py " + reference_db + " " + reciprocal_blast_file + " 1"
	os.system(command)

	#again, need a file handle
	secondFastaFilename = reciprocal_blast_file
	secondFastaFilename = secondFastaFilename.replace(".xml", ".fasta")

	reciprocalHit = next(SeqIO.parse(secondFastaFilename, 'fasta'))
	
	#check to see whether IDs are the same; if they are, add the sequence to a compilation of sequences
	#regardless, write results in summary file
	
	#we need this name regardless:
	genome = query_db.split('/')[-1].split('.')[0] 
	genome = genome.split('_')[1] + '_' + genome.split('_')[2] #genus_species
	
	#take off transcript variant portion of id
	justQueryGene = querySeq.id.split(".")[0]
	justFoundGene = reciprocalHit.id.split(".")[0]
	
	#if querySeq.id == reciprocalHit.id:
	if justQueryGene == justFoundGene:
		#add descriptive header to compiled results file so we know what genome each sequence came from 
		#gene = query_seq.split('/')[-1]
		topHit = SeqIO.parse(open(firstFastaFilename, "rU"), "fasta").next()
		#newHeader = ">" + topHit.id + " genome:" + genome + " gene:" + gene
		newHeader = ">" + genome 
		command = "echo " + "'" + newHeader + "'" + " >> " + compiled_results_file
		os.system(command)
		#also add sequence
		hitSequence = str(topHit.seq)
		command = "echo " + "'" + hitSequence + "'" + " >> " + compiled_results_file
		os.system(command)
		#finally, add line to summary file
		summaryLine = "PASS: For " + genome + ", the reciprocal blast was correct."
		command = "echo " + "\"" + summaryLine + "\"" + " >> " + summary_file
		os.system(command)
	else:
		summaryLine = "FAIL: For " + genome + ", the reciprocal blast was incorrect."
		command = "echo " + "\"" + summaryLine + "\"" + " >> " + summary_file
		os.system(command)



def main():
	#call script with the following arguments:
	#genesDirectory: directory where all of your query sequences are (each in a separate file)
	#genomesDirectory: directory where all of your genomes of interest are
	#genomeEndsWith: what each query genome ends with 
	#referenceGenome: original genome from which query sequences came; needs to be in genomesDirectory
	#outputPath: where to put the summary and compiled sequence files

	genesDirectory = sys.argv[1]
	genomesDirectory = sys.argv[2]
	genomeEndsWith = sys.argv[3]
	referenceGenome = sys.argv[4]
	outputPath = sys.argv[5] #this name needs to end with a /
	
	#pull out names of files
	genes = os.listdir(genesDirectory)
	genomes = os.listdir(genomesDirectory)

	#we only want some of the files in the genomes directory
	genomeFiles = []
	for name in genomes:
		if name.endswith(genomeEndsWith): genomeFiles.append(name)
	if referenceGenome.endswith(genomeEndsWith):
		genomeFiles.remove(referenceGenome)
	#if referenceGenome.endswith(genomeEndsWith):
		#genomeFiles.remove(referenceGenome.split('/')[-1]) #so we don't query the reference

	#make files for results using gene names
	for gene in genes:
		geneName = gene.split('.')[0] #get rid of extension
		compiledFileName = outputPath + geneName + '_compiledSeqs.txt'
		summaryFileName = outputPath + geneName + '_runSummary.txt'
		
		#now use my reciprocal Blast function
		for genome in genomeFiles:

			#parsing will need to be changed for different file naming conventions
			genomeName = genome.split(".")[0]
			genomeNameFinal = genomeName.split("_")[1] + "_" +  genomeName.split("_")[2]
			firstBlastXMLname = geneName + '_against_' + genomeNameFinal + '.xml'
			reciprocalBlastXMLname = geneName + '_' + genomeNameFinal + '_reciprocal.xml'
			singleReciprocalBlast(genesDirectory+gene, genomesDirectory+genome, genomesDirectory+firstBlastXMLname, referenceGenome, genomesDirectory+reciprocalBlastXMLname, compiledFileName, summaryFileName)
				

if __name__ == "__main__":
	main()
