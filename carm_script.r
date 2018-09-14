# Welcome to the phylogenetics lab!
# This lab will walk you through a few steps of running phylogenetic analysis in the R language
# This lab will first take multiple morphological measurements that you will provide
# and construct a phylogenetic tree!

###############
## Section 1 ##
###############


# first, set the working directory of this program to wherever you put all the files for this lab
setwd('/home/vortacs/git-repos/practice_scripts/')

# next, load all the packages you'll need to have available to run various analyses
install.packages("ape")
install.packages("phangorn")
install.packages("seqinr")
install.packages("ips")

library(ape)
library(phangorn)
library(seqinr)
library(ips)

# Next we need to load in the file of morphological character information you constructed
# this file should be a .csv file with the names of animals as the row names
# and the character names as the column names
morph_data = rbind(A=c(1,2,5,0,7,5,2,2,6,1,6,3,2,0,1),
                   B=c(2,4,0,3,1,1,1,3,0,2,6,0,2,0,1),
                   C=c(2,1,5,1,0,1,5,2,6,7,0,6,4,7,1),
                   D=c(0,3,0,4,1,0,2,2,1,0,2,1,0,2,1), 
                   E=c(2,1,1,1,0,5,6,0,1,0,0,0,0,1,1), 
                   I=c(0,4,3,3,7,1,2,3,1,0,0,0,0,0,1), 
                   G=c(1,1,1,1,1,1,1,1,1,0,0,0,1,1,1), 
                   H=c(0,2,4,5,6,0,7,3,5,2,1,0,1,0,7))
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

phydat_morph_dist = dist.hamming(phymorph_data)
phydat_morph_dist

# what do you think this distance matrix is showing you? What is the max number of character information that we input?
# with that in mind, what do the numbers look like they mean in this distance matrix?

# now we can use a few different methods
# lets start with UPGMA and Neighbor Joining trees
# UPGMA and Neighbor Joining (NJ) are both algorithm methods of phylogenetic analyses
# this is just a fancy way of saying that they optimize finding the best tree and do it once
# if you think you found the best tree on the first try, why could this be a problem?

upg = upgma(dist.gene(morph_data))
plot(upg, main="UPGMA")

stree = nj(dist.gene(morph_data))
plot(stree, main="NJ")

# We are using two different phylogenetic packages here and they specify the same functions
# in different ways
# plug this in to call NJ methods on your data with the Phangorn package, which we will use more

nj_morph = NJ(phydat_morph_dist)

# Ok, lets calculate parsimony changes and plot a tree with that data.

treepars_morph = optim.parsimony(nj_morph, phymorph_data, method = "fitch" )

# Take a look at the different parsimony trees that are considered equally parsimonious
# what is Parsimony and how can it provide phylogenetic information?

trees = bab(phymorph_data, tree = NULL)
plot(trees)

# This command roots a specific group as the OUTGROUP. Change this value to change which species
# is the outgroup

treepars_morph.rooted = root(treepars_morph, outgroup = "I")

# Plot your tree, noting that you can change the kind of tree you display
# check the Phangorn documentation for more options

plot(treepars_morph.rooted , "phylogram")

# Ok, the big moment. Bootstrap your analysis to perform is multiple times
# and get information about how often we find a particular topology

bs_morph.set = bootstrap.phyDat(phymorph_data, pratchet, bs = 100)

# Plot the resuts using a special bootstrap plotting function
# what do those numbers mean?

plotBS(treepars_morph, bs_morph.set, p = 10, type = "phylogram", bs.col="black")

# set a root for these bootstrapped trees
# Here we construct the consensus tree from all of our bootstraps and display
# information supporting particular topologies

bs_morph.set.rooted = root(bs_morph.set, outgroup = "I")

plotBS(treepars_morph.rooted, bs_morph.set.rooted, p = 10, type = "phylogram")

# Awesome, you're on your way to being an evolutionary biologist!
# Ok, lets try some other stuff

###############
## Section 2 ##
###############

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

#ape_bs = boot.phylo(phy=camin_NJ, x=seqdata, NJ, trees=TRUE)


# Ok, take a look at the consensus tree. See those numbers?
# If we recalculated the tree 100 times, what do those numbers mean?
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")


# This command will write your tree to a file
# if you're in Rstudio, you can use the EXPORT button as well
write.tree(bs, file="bootstrap_example.tre")









seqdata = read.dna("cytb_align.fasta", format="fasta")

# Take a look at the data
# this is a good way to understand how data is formatted

seqdata

# Count parsimony informative sites

pis(seqdata, what ="absolute", use.ambiguities = FALSE)

camin_phyDat <- phyDat(seqdata, type = "DNA", levels = NULL)

dm = dist.ml(camin_phyDat)

treeNJ = NJ(dm)

plot(treeNJ, "phylogram", main = "NJ")

treepars = optim.parsimony(treeNJ, camin_phyDat, method = "fitch" )

treepars.rooted = root(treepars, outgroup = "Balaenoptera")
plot(treepars.rooted , "phylogram")

bs.set = bootstrap.phyDat(camin_phyDat, pratchet, bs = 100)

bs.set.rooted = root(bs.set, outgroup = "Carcharodon")
treepars.rooted = root(treepars, outgroup = "Carcharodon")

summary(bs.set)


plotBS(treepars.rooted, bs.set.rooted, type = "phylogram")


bs.set.root = superTree(bs.set, rooted = TRUE)

plotBS(treepars, bs.set.root, type = "phylogram", bs = 30)

plotBS(treepars, bs.set.root, bs = 50, type = "phylogram")









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
