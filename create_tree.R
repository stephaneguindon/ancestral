# The script below was written by A. Oliva

library(TreeSim)
library(ape)

#Recover the argument from the Python Script to determine the number of tips
nb_of_tips<-commandArgs(trailingOnly = TRUE)
nb_of_tips=as.integer(nb_of_tips)

#This was the two ways we decided the depth of our trees
#Born_sup=(nb_of_tips/2)*10^-1
#Born_inf=(nb_of_tips/2)*10^-2
Born_inf=0.1
Born_sup=1.0
Tree_Depth=runif(1,Born_inf,Born_sup)

#Using "TreeSim" to create a rooted tree
Arbre1=sim.bd.taxa.age(n = nb_of_tips,numbsim = 1, lambda = 10.0, mu = 0.5, frac = 1, age = Tree_Depth, mrca = TRUE)
#unroot it thanks to "ape"
unrooted_tr=unroot(Arbre1[[1]])

#Output them
write.tree(unrooted_tr, file = "Unrooted_Tree.tree")
write.tree(Arbre1[[1]], file = "Rooted_Tree.tree")
