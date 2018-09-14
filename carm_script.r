# Welcome to the phylogenetics lab!
# This lab will walk you through a few steps of running phylogenetic analysis in the R language
# This lab will first take multiple morphological measurements that you will provide
# and construct a phylogenetic tree!

# first, set the working directory of this program to wherever you put all the files for this lab
setwd('/home/vortacs/git-repos/practice_scripts/')

# next, load all the packages you'll need to have available to run various analyses
install.packages("ape")
install.packages("phangorn")
install.packages("seqinr")

library(ape)
library(phangorn)
library(seqinr)

# Next we need to load in the file of morphological character information you constructed
# this file should be a .csv file with the names of animals as the row names
# and the character names as the column names
morph_data = rbind(A=c(1,2,5,0,7,5,2,2,6,0),
                   B=c(2,4,0,3,1,1,1,3,0,0),
                   C=c(2,1,5,1,0,1,5,2,6,7),
                   D=c(0,3,0,4,1,0,2,2,1,0), 
                   E=c(2,1,1,1,0,5,6,0,1,0), 
                   I=c(0,4,3,3,7,1,2,3,1,0), 
                   G=c(1,1,1,1,1,1,1,1,1,0), 
                   H=c(0,0,0,0,0,0,0,0,0,0))
head(morph_data)

# Before you change this data and make it your own, make a tree using the commands below
# then try changing a few numbers in various columns and observe what happens to the tree

upg = upgma(dist.gene(morph_data))
plot(upg, main="UPGMA")

# Ok back on track! zero (0) out the values in the dataset and record your own measurements
# If you don't, you'll get the wrong tree!
# then import your data into a phyDat format. This is necessary for some of the commands later on

phymorph_data <- as.phyDat(morph_data, type = "USER", levels = c(0:7), ambiguity = c('-'), header = "TRUE")


# ok Lets just make sure that worked and run a quick neighbor joining tree
# first we have to construct a distance matrix to show the phylogenetic distance
# between all the organisms

morph_dist = dist.gene(morph_data)
morph_dist

# what do you think this distance matrix is showing you? What is the max number of character information that we input?
# with that in mind, what do the numbers look like they mean in this distance matrix?

# now we can use a few different methods
# lets start with UPGMA and Neighbor Joining trees
# these are algorithmic methods that focus on getting the best tree ###########expand################

upg = upgma(dist.gene(morph_data))
plot(upg, main="UPGMA")

stree = nj(dist.gene(morph_data))
plot(stree, main="NJ")

# Do these trees have the same topology?

# a parsimony check will tell you if your data is in the phyDat format
pars = parsimony(upg, phymorph_data)
pars

pscore = parsimony(stree, phymorph_data, method="fitch")
pscore

# Awesome, you're on your way to being an evolutionary biologist!
# Ok, lets try some other stuff

# This command will give you a rough parsimony tree
# what is Parsimony and how can it provide phylogenetic information?
trees = bab(phymorph_data, tree = NULL)
plot(trees)

# did you hit return and see all the trees this made?
# so many rearrangements!

# This command will give you the Parsimony score of your tree
pscore = parsimony(trees, phymorph_data, method="fitch")
pscore


mbs = boot.phylo(stree, morph_dist, FUN=function(nj), B=100, trees=TRUE)


# Ok, this is the New Age of Genomics, lets use some of that data to construct a tree
# Load your sequence data into R
# sequence data is normally in the form of a fasta file
#seqdata = read.dna("caminalcules_seqs.fasta", format="fasta")
seqdata = read.dna("cytb_align.fasta", format="fasta")

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
#dna_dist <- dist.ml(camin_phyDat)

# Do a little internet digging and find out what the JC69 model entails

# ok now lets use some different types of analyses to look at our data
# UPGMA and Neighbor Joining (NJ) are both algorithm methods of phylogenetic analyses
# this is just a fancy way of saying that they optimize finding the best tree and do it once
# if you think you found the best tree on the first try, why could this be a problem?


camin_UPGMA <- upgma(dna_dist)

camin_NJ  <- NJ(dna_dist)

plot(camin_UPGMA, main="UPGMA")

plot(camin_NJ, main = "Neighbor Joining")

# Great! now lets check the parsimony scores of this new, inredible dataset we have

seqpars = parsimony(camin_UPGMA, camin_phyDat, method="fitch")
seqpars

seqpscore = parsimony(camin_NJ, camin_phyDat, method="fitch")
seqpscore

# Lets take a look at the various parsimony tree arrangements that are considered equally good
seq_trees = bab(camin_phyDat, tree = camin_NJ)
plot(seq_trees)

# Now we can begin using the nucleotide substitution model we identified to calculate the maximum likelihood of a topology
# We'll need to do this a few times to give us a measure of our certainty of this arrangement

fit <- pml(camin_NJ, camin_phyDat)

print(fit)

fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")

bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))

# Ok, take a look at the tree. See those numbers?
# If we recalculated the tree 100 times, what do those numbers mean?
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")


# This command will write your tree to a file
# if you're in Rstudio, you can use the EXPORT button as well
write.tree(bs, file="bootstrap_example.tre")







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
