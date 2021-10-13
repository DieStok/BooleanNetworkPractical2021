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
net <- loadSBML("../MODEL_Azpeitaetal.xml",symbolic=FALSE)

#reorder genes to paper order
netReorder <- funcReorderModelToPaperOrder(net)

#look at the network. Shows involved genes and representations of their Boolean rules.
print(netReorder)

#plot how the network is wired
plotNetwork <- plotNetworkWiring(netReorder, simplify = FALSE)

#You can also upload the model to BooleSim (https://rumo.biologie.hu-berlin.de/boolesim/) in Google Chrome and toy with it there.
saveNetwork(netReorder, "plantGRNBoolean.txt") #you can now upload this .txt file to BooleSim

#We can plot the positive and negative interactions in R, but then we lose the self-loops.
#Ignore the warnings, they are for the self-loops that can't be drawn.
vectorInhibitionOrActivation <- c(1, 1, -1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1)
plot_directed(plotNetwork, cex.label = 0.5, state = vectorInhibitionOrActivation,
              fill.node = "lightblue", cex.node = 2, arrow_clip = 0.15)

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
#from nodes you can not see which genes are on or off
plotStateGraph(att)

#To fix certain genes on or off in the network, use this syntax:
someNetWithFixations <- fixGenes(netReorder, fixIndices = funcSetNames(netReorder, c("WOX5", "CLEX")),
                                 values = c(0, 1))
print(someNetWithFixations) #as you can see, WOX5 is now fixed at 0, and CLEX at 1.
attSomeNetWithFixations <- getAttractors(someNetWithFixations)
plotAttractors(attSomeNetWithFixations)



#########################################################################################
#                                                                                       #
#                          Write your own altered code below                            #
#                                                                                       #
#########################################################################################



#Question 1
#plot network
plotNetwork <- plotNetworkWiring(netReorder, simplify = FALSE)

#Question 2
#simply use print(att) to answer this question
print(att)

#Question 4
fixedNetAuxinSHRZero <- fixGenes(netReorder, fixIndices = funcSetNames(netReorder, c("AUXIN", "SHR")),
                     values = c(0, 0))
newAttractorsKnockoutAuxinAndSCR <- getAttractors(fixedNetAuxinSHRZero)
print(newAttractorsKnockoutAuxinAndSCR)

#Question 6
#plot state space
plotStateGraph(att)

#Question 7
print(att) #look in here for the attractors corresponding to those of the paper.

#Question 8
print(netReorder) #look at the rule for Wox5
#can check with fixGenes to be sure
fixSHROnSCROnCLEXOff <- fixGenes(netReorder,
                                 funcSetNames(netReorder,
                                              geneNames = c("SHR", "SCR", "CLEX")),
                                 values = c(1, 1, 0))
attFixedSHROnSCRONCLEXOff <- getAttractors(fixSHROnSCROnCLEXOff)
plotAttractors(attFixedSHROnSCRONCLEXOff)

#Question 9
data("cellcycle")
print(cellcycle)
plotNetworkWiring(cellcycle)
cellCycleAtt <- getAttractors(cellcycle)
plotAttractors(cellCycleAtt)
print(cellCycleAtt)
plotStateGraph(cellCycleAtt)

#Question 11
cellCycleAsynchronousAtt <- getAttractors(cellcycle, type = "asynchronous", method = "random", startStates = 800)
plotAttractors(cellCycleAsynchronousAtt)


#Question 12 Asynchronous updating plant Boolean network
attAsynchronous <- getAttractors(network = netReorder,type="asynchronous", startStates=1024, method = "random")
plotAttractors(attAsynchronous)
plotAttractors(att)
print(attAsynchronous)
print(att)


#Question 13 
fixedSHROff <- fixGenes(netReorder, fixIndices = funcSetNames(netReorder, "SHR"), values = c(0))
attSHROff   <- getAttractors(fixedSHROff, type = "synchronous")
plotAttractors(attSHROff, drawLegend = FALSE, title = "Network SHR knocked out")
plotAttractors(att, drawLegend = FALSE, title = "Normal network B")
print(attSHROff)


#Question 14: perturb the network
#picks a random gene and a random function of that gene (for example 0101 --> 1, now becomes 0101 --> 0, only for that gene)
#If you do it 3 times, you will get a different result each time.
perturbedNet <- perturbNetwork(netReorder, perturb = "functions", method = "bitflip", readableFunctions = TRUE)
plotNetworkWiring(perturbedNet)
perturbedNet2 <- perturbNetwork(netReorder, perturb = "functions", method = "bitflip", readableFunctions = TRUE)
plotNetworkWiring(perturbedNet2)
perturbedNet3 <- perturbNetwork(netReorder, perturb = "functions", method = "bitflip", readableFunctions = TRUE)
plotNetworkWiring(perturbedNet3)


#########################################################################################
#                                                                                       #
#                             Optional final questions                                  #
#                                                                                       #
#########################################################################################

#Question 16: Perturb the network randomly 70 times and save the network and its attractors.
set.seed(666999)
#make empty lists
perturbedNetList = list()
perturbedNetAttList = list()
  for (perturbation in seq(1, 70, 1)) {
    #generate perturbed network and its attractors
    perturbedNet    <- perturbNetwork(netReorder, perturb = "functions", method = "bitflip", readableFunctions = TRUE)
    attPerturbedNet <- getAttractors(perturbedNet)
    
    #store in respective lists
    perturbedNetList[[as.character(perturbation)]]    = perturbedNet
    perturbedNetAttList[[as.character(perturbation)]] = attPerturbedNet
  }

#sapply applies a certain function to each element of a list, and returns the most simple representation of the results (here: vector)
#the function here simply gets the amount of attractors by seeing how many attractors are listed in the stateInfo$attractorAssignment
attractorsPerPerturbation = sapply(perturbedNetAttList,
                                   FUN = function(currentAttractor) {max(currentAttractor$stateInfo$attractorAssignment) }
                                   ) 

#Doing the sapply with a loop in 3 different ways
attractorsPerPerturbationLoop1 = c()
attractorsPerPerturbationLoop2 = c()
attractorsPerPerturbationLoop3 = c()

for (attractor in perturbedNetAttList) {
  
  amountOfAttractors = max(attractor$stateInfo$attractorAssignment)
  #other options for amountOfAttractors, there are multiple ways to get it from the object:
  amountOfAttractors2 = length(unique(attractor$stateInfo$attractorAssignment))
  amountOfAttractors3 = length(attractor$attractors)
  
  attractorsPerPerturbationLoop1 = c(attractorsPerPerturbationLoop1, amountOfAttractors)
  attractorsPerPerturbationLoop2 = c(attractorsPerPerturbationLoop2, amountOfAttractors2)
  attractorsPerPerturbationLoop3 = c(attractorsPerPerturbationLoop3, amountOfAttractors3)
}

print(attractorsPerPerturbationLoop1)
print(attractorsPerPerturbationLoop2)
print(attractorsPerPerturbationLoop3)

#tabulate results
print(table(attractorsPerPerturbation))

#Question 18: if 2 attractors only, what has happened?
#in these cases, SHR function has been changed
indicesTwoAttractors = which(attractorsPerPerturbation == 2)
for(networkTwoAttractors in indicesTwoAttractors) {
  print(perturbedNetList[[networkTwoAttractors]]$interactions$SHR)
  
}
#original function for comparison
print(netReorder$interactions$SHR)

#look at the attractors
getAttractors(perturbedNetList[[indicesTwoAttractors[1]]])



#Question 17. What has happened if you get 12 attractors?
indexTwelveAttractors = which(attractorsPerPerturbation == 12)
#look only at one example
if(length(indexTwelveAttractors) > 1) {
  indexTwelveAttractors = indexTwelveAttractors[1]
}
print(perturbedNetList[[indexTwelveAttractors]]) #see CLEX is just dependent on CLEX now.
print(perturbedNetList[[indexTwelveAttractors]]$interactions$CLEX)
plotAttractors(perturbedNetAttList[[indexTwelveAttractors]])

#compare with normal Net
print(netReorder)
print(netReorder$interactions$CLEX)

