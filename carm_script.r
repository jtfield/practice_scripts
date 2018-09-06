# Welcome to the phylogenetics lab!
# This lab will walk you through a few steps of running phylogenetic analysis in the R language
# This lab will first take multiple morphological measurements that you will provide
# and construct a phylogenetic tree!

# first, set the working directory of this program to wherever you put all the files for this lab
setwd('/home/vortacs/Documents/evolution_fall_2018/phylo-lab/')

# next, load all the packages you'll need to have available to run various analyses
install.packages("ape")
install.packages("phangorn")
install.packages("seqinr")

library(ape)
library(phangorn)
library(seqinr)

# Next we need to load in the file of morphological traits you constructed
# this file should be a .csv file with the names of animals as the column headers
# and the trait names as the row names
data = read.csv(file="carm_test.csv", header=TRUE, sep=",")

# Now we need to transform the .csv file into a data matrix that the program can understand
morphData = as.phyDat(data, type="USER", levels = c(0, 1))
morphData

# This command will give you a rough parsimony tree
# what is Parsimony and how can it provide phylogenetic information?
trees = bab(morphData, tree = NULL, trace = 1)
plot(trees)

# This command will give you the Parsimony score of your tree
pscore = parsimony(trees, morphData, method="fitch")
pscore

# Ok, this is the New Age of Genomics, lets use some of that data to construct a tree
# Load your sequence data into R
# sequence data is normally in the form of a fasta file
seqdata = read.dna("caminalcules_seqs.fasta", format="fasta")

# Take a look at the data
# this is a good way to understand how data is formatted
seqdata
# what do you think the numbers under the nucleotide labels mean? Do those numbers add up to 1?

# Ok now put this data in another data matrix
# Notice that some of the descriptive "flags" have changed
camin_phyDat <- phyDat(seqdata, type = "DNA", levels = NULL)

# First things first, we need a nucleotide substitution model that fits our data
mt <- modelTest(camin_phyDat)

# Lets look at what we get
print(mt)

# Ok, we're picking the simplest model for our analyses to keep things simple
# the JC69 model
# now create a distance matrix. What do you think a distance matrix is?
dna_dist <- dist.ml(camin_phyDat, model="JC69")
# Do a little internet digging and find out what the JC69 model entails

# ok now lets use some different types of analyses to look at our data
# UPGMA and Neighbor Joining (NJ) are both algorithm methods of phylogenetic analyses
# this is just a fancy way of saying that they optimize finding the best tree and do it once
# if you think you found the best tree on the first try, why could this be a problem?
camin_UPGMA <- upgma(dna_dist)

camin_NJ  <- NJ(dna_dist)

plot(camin_UPGMA, main="UPGMA")

plot(camin_NJ, main = "Neighbor Joining")

parsimony(camin_UPGMA, camin_phyDat)

parsimony(camin_NJ, camin_phyDat)

camin_optim <- optim.parsimony(camin_NJ, camin_phyDat)

camin_pratchet <- pratchet(camin_phyDat)

plot(camin_optim)

plot(camin_pratchet)

fit <- pml(camin_NJ, camin_phyDat)

print(fit)

fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")

logLik(fitJC)

bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))

plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
