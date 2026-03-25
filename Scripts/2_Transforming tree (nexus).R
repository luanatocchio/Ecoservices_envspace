library(ape)

# Loading the phylogenetic tree

load("Data\\Primary Data\\final_tree.RData")

# Writing the nexus file

write.nexus(final_tree, file = "final_tree.nex")
