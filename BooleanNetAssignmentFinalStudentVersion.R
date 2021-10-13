#Boolean Networks Practical

#########################################################################################
#                                                                                       #
#   Loading packages and defining required functions. Just run this and continue (:     #
#                                                                                       #
#########################################################################################



if(!require('pacman'))install.packages('pacman')
library(pacman)
p_load(BoolNet, XML, igraph, graphsim)
set.seed(52525)

#article on BoolNet
#https://academic.oup.com/bioinformatics/article/26/10/1378/193238

#article on model for gene regulatory network for plant root stem cell
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2972269/


#functions. Do not worry about these, they are simply there for your convenience and used at appropriate places in the script.
#No need to look through them.

#function to define the indices of genes to fix in the reordered network.
funcSetNames <- function(net, geneNames = c("AUXIN")) {
  which(names(net$genes) %in% geneNames)
}
#function to change the network order so that it is the same as in the paper
funcReorderModelToPaperOrder <- function(net){
  
  paperOrder = c("PLT", "AUXIN", "ARF", "AUXIAA", "SHR", "SCR", "JKD", "MGP", "WOX5", "CLEX")
  netReorder = net; netReorder$genes = paperOrder; names(netReorder$genes) = paperOrder
  geneNumberAndName = seq(1,10,1); names(geneNumberAndName) = net$genes
  geneNumberAndNameNewOrder = geneNumberAndName[paperOrder]
  
  #need to change numeric indexes of genes per interaction, and list ordering of interaction list.
  for (gene in names(netReorder$interactions)) {
    
    newNumbers = c()
    for (number in netReorder$interactions[[gene]]$input) {
      numberToAdd = which(geneNumberAndNameNewOrder == number)
      newNumbers = c(newNumbers,numberToAdd)
      
    }
    netReorder$interactions[[gene]]$input = newNumbers
    
  }
  
  
  newList = netReorder$interactions[names(netReorder$genes)]
  netReorder$interactions = newList
  
  return(netReorder)
}



#########################################################################################
#                                                                                       #
#   Practical start. First walk through these functions and see what they do.           #
#                                                                                       #
#########################################################################################



#Load in the network. Corresponds to model B from the plant root SCN publication.
#Model was downloaded from the biomodels database here: https://www.ebi.ac.uk/biomodels/MODEL1504170002
net <- loadSBML("MODEL_Azpeitaetal.xml",symbolic=FALSE)

#reorder genes to paper order
netReorder <- funcReorderModelToPaperOrder(net)

#look at the network. Shows involved genes and representations of their Boolean rules.
print(netReorder)

#plot how the network is wired
plotNetwork <- plotNetworkWiring(netReorder, simplify = FALSE)
#You can also upload the model to BooleSim (https://rumo.biologie.hu-berlin.de/boolesim/) in Google Chrome and toy with it there.
saveNetwork(netReorder, "plantGRNBoolean.txt") #you can now upload this .txt file to BooleSim in Google Chrome
#We can plot the positive and negative interactions in R, but then we lose the self-loops.
#Ignore the warnings, they are for the self-loops that can't be drawn.
vectorInhibitionOrActivation <- c(1, 1, -1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1)
plot_directed(plotNetwork, cex.label = 0.5, state = vectorInhibitionOrActivation,
              fill.node = "lightblue", cex.node = 2, arrow_clip = 0.15)
#be sure to zoom in on the plot to see it better.

#List attractors and their basin of attraction. Plot as a table or print in the console.

att <-getAttractors(netReorder)
plotAttractors(att)
print(att)

#plots attractors and their basins of attraction
#each state is a node, each transition to next state an edge
#different attractors and their basins get different colors
#attractor states are indicated by having clearer outline
#of their nodes, cyclic attractors can be seen to have clearer
#transition arrows between their states
#from nodes you cannot see which genes are on or off
plotStateGraph(att)

#To fix certain genes on or off in the network, use this syntax:
someNetWithFixations <- fixGenes(netReorder, fixIndices = funcSetNames(netReorder, c("WOX5", "CLEX")),
                                 values = c(0, 1))
print(someNetWithFixations) #as you can see, WOX5 is now fixed at 0, and CLEX at 1.
attSomeNetWithFixations <- getAttractors(someNetWithFixations)
plotAttractors(attSomeNetWithFixations)
#Note that to answer question 4 you need to change this code!

#########################################################################################
#                                                                                       #
#                          Write your own altered code below                            #
#                                                                                       #
#########################################################################################


#Question 4



#Question 8



#Question 9
data("cellcycle")


#Question 11



#Question 12



#Question 13



#Question 14  perturb the network


#########################################################################################
#                                                                                       #
#                             Optional final questions                                  #
#                                                                                       #
#########################################################################################

#Only do these if you want a small challenge!

#Question 16: Perturb the network randomly 70 times and save the network and its attractors.
set.seed(666999)
#make empty lists
perturbedNetList = list()
perturbedNetAttList = list()
  for (perturbationNumber in seq(1, 70, 1)) {
    #Your Code Here!
  }

#Question 17