#makeTreeThings.R
#Script to 1.make a pretty tree from a .nwk file and 2.produce a tree and tables for use in
#decoding which nodes are which
 
 
library(ape)
args=(commandArgs(TRUE))  

#call script with the following arguments: 
treeFilename = args[1] #.nwk file name
outputPrefix = args[2] #prefix to use on all output files

#read in tree
tree = read.tree(treeFilename)

#make the good looking tree
pdfName = paste(outputPrefix, '.pdf', sep="")
pdf(pdfName)
plot(tree)
dev.off()

#make a tree with all the nodes labeled
pdfName = paste(outputPrefix, '_nodeNumbers.pdf', sep="")
pdf(pdfName)
plot(tree)
nodelabels()
tiplabels()
dev.off()

#make a table of names with corresponding numbers
names = c(tree$tip.label)
namesDf = data.frame(names)
write.table(namesDf, file = paste(outputPrefix,"_names_numbers.txt", sep=""), quote=FALSE, sep="\t", col.names = FALSE)

#make a table of node numbers with corresponding branch lengths
nodeCoorsTable = data.frame(tree$edge)
nodeCoorsTable = cbind(nodeCoorsTable, tree$edge.length)
colnames(nodeCoorsTable) = c("node1","node2", "edge_len")
write.table(nodeCoorsTable, file = paste(outputPrefix,"_edge_lengths.txt", sep=""), quote=FALSE, sep="\t", col.names=NA)

            
            
            