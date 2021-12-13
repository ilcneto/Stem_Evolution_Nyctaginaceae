# Character dependent diversification for Nyctaginaceae
# using Hidden State Speciation and Extinction (HiSSE; Beaulieu & O'Meara 2016)

install.packages("hisse")
install.packages("deSolve")
install.packages("GenSA")
install.packages(c("subplex", "nloptr"))

lapply(X=c("hisse","ape","geiger", "deSolve", "GenSA","subplex","nloptr"),FUN = library, character.only=TRUE)
suppressWarnings(library(diversitree))

# Trait-dependent diversification in Nyctaginaceae

setwd("~/Documents/Proyectos/Nyctaginaceae/")

# reading tree
tree<-read.tree("nycta_mcc_newick.tree")
tree

# reading data file
data<-read.csv("Nycta_character_matrix_ILCN.csv")
data

##### habit #####
habit<-data.frame(data$taxon, data$habit, row.names = data$taxon)
habit

# check if tree tip names and data names match
name.check(tree, data = habit) #OK

# Model 1. Dull-null: turnover and extinction fraction equal for the two states
# no hidden character.
turnover <- c(1,1) #equal turnover rate
extinction.fraction<-c(1,1) #equal extinction fraction

#sampling fraction: proportion of species sampled that have 0 and 1 (in this order),
#considering the TOTAL existing species that have each state.
sum(habit == 0) # 82. So proportion is 82/552=0.15
sum(habit == 1) # 22; 22/45=0.49
f<-c(0.15,0.49)

# generating a transition matrix (q)
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)
#    (0) (1)
#(0)  NA   2
#(1)   1  NA

# running Model 1
dull.null_habit <- hisse(phy=tree, data=habit, f=f, turnover=turnover, 
                   eps=extinction.fraction, hidden.states=FALSE,
                   trans.rate=trans.rates.bisse)


# Model 2. BiSSE: turnover different; extinction fraction equal for the two states.
turnover <- c(1,2) #different turnover rate
extinction.fraction<-c(1,1) #equal extinction fraction
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0) # different transition rates

# running Model 2
BiSSE_habit <- hisse(phy=tree, data=habit, f=f, turnover=turnover,
               eps=extinction.fraction, hidden.states=FALSE,
               trans.rate=trans.rates.bisse)
# warning messages, OK
BiSSE_habit


# Model 3. HiSSE: turnover different for all the states (observed and hidden).
turnover <-c(1,2,3,4)
extinction.fraction<- c(1,1,1,1)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits = 1)
print(trans.rate.hisse)

# running Model 3
HiSSE_habit <- hisse (phy=tree, data=habit, f=f, turnover=turnover,
                eps=extinction.fraction, hidden.states=TRUE,
                trans.rate=trans.rate.hisse)

HiSSE_habit

# Model 4. CID2: turnover different for the hidden states and equal for the observed states,
# that is, character-independent diversification.
turnover<- c(1,1,2,2)
extinction.fraction<-c(1,1,1,1)
trans.rate<-TransMatMakerHiSSE(hidden.traits = 1, make.null = TRUE) #Three rates
trans.rate
#     (0A) (1A) (0B) (1B)
#(0A)   NA    2    3   NA
#(1A)    1   NA   NA    3
#(0B)    3   NA   NA    2
#(1B)   NA    3    1   NA


CID2_habit<-hisse(phy=tree, data=habit, f=f, turnover=turnover, eps=extinction.fraction, 
           hidden.states=TRUE,trans.rate=trans.rate)
CID2_habit



# Model 5. CID4: turnover different for the hidden states but there are two additional
# states (C and D), equal turnover for the observed states, it has the same number of 
# parameters than the Model 3 (HiSSE) so this evaluates the character-independent 
# model as a null model to compare the HiSSE model.

turnover<-c(1,1,2,2,3,3,4,4)
extinction.fraction<-rep(1,8)
trans.rate<-TransMatMakerHiSSE(hidden.traits = 3, make.null = TRUE)
trans.rate
#      (0A) (1A) (0B) (1B) (0C) (1C) (0D) (1D)
#(0A)   NA    2    3   NA    3   NA    3   NA
#(1A)    1   NA   NA    3   NA    3   NA    3
#(0B)    3   NA   NA    2    3   NA    3   NA
#(1B)   NA    3    1   NA   NA    3   NA    3
#(0C)    3   NA    3   NA   NA    2    3   NA
#(1C)   NA    3   NA    3    1   NA   NA    3
#(0D)    3   NA    3   NA    3   NA   NA    2
#(1D)   NA    3   NA    3   NA    3    1   NA

CID4_habit<-hisse(phy=tree, data=habit, f=f, turnover=turnover, eps=extinction.fraction, 
            hidden.states=TRUE,trans.rate=trans.rate)
CID4_habit


# Summarize the results of the models to compare the likelihoods and select models
library(tidyr)

models_habit<- as.data.frame(rbind(Null=c(dull.null_habit$loglik, dull.null_habit$AICc),
                             BiSSE=c(BiSSE_habit$loglik, BiSSE_habit$AICc),
                             HiSSE=c(HiSSE_habit$loglik, HiSSE_habit$AICc),
                             CID2=c(CID2_habit$loglik, CID2_habit$AICc),
                             CID4=c(CID4_habit$loglik, CID4_habit$AICc)))
models_habit                             
colnames(models_habit) <- c("loglik", "AICc")                             
                    
# order models in the dataframe
models_habit<-models_habit[order(models_habit$AICc), ]
models_habit
#loglik     AICc
#CID2  -355.4073 723.6806
#HiSSE -350.7401 723.8457
#CID4  -354.1289 725.7736
#Null  -366.5628 741.5296
#BiSSE -366.5258 743.6639

# obtaining Akaike values
aic <- models_habit$AICc

# selecting the best model (i.e., minimum value of Akaike)
best <- min(aic)

# obtaining the difference between Akaike values
delta <- aic - best

#obtaining Akaike weights by first calculating the relative likelihood of the models using delta
sumDelta <- sum(exp(-0.5 * delta))

# The Akaike weight for a model is this value divided by the sum of these values across all models.
# http://brianomeara.info/aic.html
w_models_habit <- round((exp(-0.5 * delta)/sumDelta), 3) %>%
  as.matrix()
colnames(w_models_habit) <- c("AICw")
w_models_habit

# saving into a matrix
models_habit_2 <- cbind(models_habit, w_models_habit) %>%
  as.data.frame() %>%
  write.csv("hisse_models_habit.csv")



##### eustele type #####

#eustele data
eustele<-data.frame(data$taxon, data$eustele_type, row.names = data$taxon)
head(eustele)

name.check(tree, eustele)
#OK

# Model 1. Null: turnover and extinction fraction equal for the two states
# no hidden character.
turnover <- c(1,1) #equal turnover rate
extinction.fraction<-c(1,1) #equal extinction fraction

#sampling fraction: proportion of species sampled from the total existing species having
# 0 and 1 states.
sum(eustele == 0) # 13, so 13/167=0.07
sum(eustele == 1) # 91, 91/430=0.21
f<-c(0.07,0.21)

trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)
#    (0) (1)
#(0)  NA   2
#(1)   1  NA

#running Model 1. Null
Null_eustele <- hisse(phy=tree, data=eustele, f=f, turnover=turnover, 
                   eps=extinction.fraction, hidden.states=FALSE,
                   trans.rate=trans.rates.bisse)
Null_eustele

# Model 2. BiSSE: turnover different; extinction fraction equal for the two states.
turnover <- c(1,2) #different turnover rate
extinction.fraction<-c(1,1) #equal extinction fraction
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0) # different transition rates

# running Model 2
BiSSE_eustele <- hisse(phy=tree, data=eustele, f=f, turnover=turnover,
               eps=extinction.fraction, hidden.states=FALSE,
               trans.rate=trans.rates.bisse)
BiSSE_eustele

# Model 3. HiSSE: turnover different for all the states (observed and hidden).
turnover <-c(1,2,3,4)
extinction.fraction<- c(1,1,1,1)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits = 1)
print(trans.rate.hisse)

# running Model 3
HiSSE_eustele <- hisse (phy=tree, data=eustele, f=f, turnover=turnover,
                eps=extinction.fraction, hidden.states=TRUE,
                trans.rate=trans.rate.hisse)


# Model 4. CID2: turnover different for the hidden states and equal for the observed states,
# that is, character-independent diversification.
turnover<- c(1,1,2,2)
extinction.fraction<-c(1,1,1,1)
trans.rate<-TransMatMakerHiSSE(hidden.traits = 1, make.null = TRUE) #Three rates
trans.rate

CID2_eustele<-hisse(phy=tree, data=eustele, f=f, turnover=turnover, eps=extinction.fraction, 
            hidden.states=TRUE,trans.rate=trans.rate)

CID2_eustele

# Model 5. CID4: turnover different for the hidden states but there are two additional
# states (C and D), equal turnover for the observed states, it has the same number of
# parameters than the Model 3 (HiSSE) so this evaluates the character-independent 
# model as a null model to compare the HiSSE model.

turnover<-c(1,1,2,2,3,3,4,4)
extinction.fraction<-rep(1,8)
trans.rate<-TransMatMakerHiSSE(hidden.traits = 3, make.null = TRUE)
trans.rate

CID4_eustele<-hisse(phy=tree, data=eustele, f=f, turnover=turnover, eps=extinction.fraction, 
            hidden.states=TRUE,trans.rate=trans.rate)



# Summarize the results of the models to compare the likelihoods and select models
models_eustele<- as.data.frame(rbind(Null=c(Null_eustele$loglik, Null_eustele$AICc),
                             BiSSE=c(BiSSE_eustele$loglik, BiSSE_eustele$AICc),
                             HiSSE=c(HiSSE_eustele$loglik, HiSSE_eustele$AICc),
                             CID2=c(CID2_eustele$loglik, CID2_eustele$AICc),
                             CID4=c(CID4_eustele$loglik, CID4_eustele$AICc)))
models_eustele                           
colnames(models_eustele) <- c("loglik", "AICc")                             

# order models in the dataframe
models_eustele<-models_eustele[order(models_eustele$AICc), ]
models_eustele
#        loglik     AICc
#BiSSE -316.9860 644.5842
#CID2  -318.7046 650.2752
#Null  -322.1928 652.7897
#CID4  -317.7696 653.0550
#HiSSE -315.5234 653.4123

# obtaining Akaike values
aic <- models_eustele$AICc

# selecting the best model (i.e., minimum value of Akaike)
best <- min(aic)

# obtaining the difference between Akaike values
delta <- aic - best

#obtaining Akaike weights by first calculating the relative likelihood of the models using delta
sumDelta <- sum(exp(-0.5 * delta))

# The Akaike weight for a model is this value divided by the sum of these values across all models.
# http://brianomeara.info/aic.html
w_models_eustele <- round((exp(-0.5 * delta)/sumDelta), 3) %>%
  as.matrix()
colnames(w_models_eustele) <- c("AICw")
w_models_eustele

# saving into a matrix
models2_eustele <- cbind(models_eustele, w_models_eustele) %>%
  as.data.frame() %>%
  write.csv("hisse_models_eustele.csv")



##### secondary growth #####

#secondary growth data
sec_gro<-data.frame(data$taxon, data$secondary_growth, row.names = data$taxon)
head(sec_gro)

# Given that three species have unknown (i.e. "?") state, we will prune the tree to remove these species
# and continue the analyses.
sp_to_drop<-c("Tripterocalyx_carneus", "Tripterocalyx_crux_maltae", "Tripterocalyx_micranthus")
tree_prunned<-drop.tip(tree, tip = sp_to_drop)
tree_prunned # Phylogenetic tree with 101 tips

# remove the same species from the data frame
sec_grow2<-sec_gro[1:101,]
sec_grow2

#checking names
name.check(tree_prunned, sec_grow2)
#OK

# Model 1. Dull-null: turnover and extinction fraction equal for the two states
# no hidden character.
turnover <- c(1,1) #equal turnover rate
extinction.fraction<-c(1,1) #equal extinction fraction

#sampling fraction
sum(sec_grow2 == 0) # 3, so the proportion is 3/28=0.1
sum (sec_grow2 == 1) # 98, 98/569=0.17
f<-c(0.1, 0.17)

#transition matrix
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)

#running Model 1. Null
Null_sec_grow <- hisse(phy=tree_prunned, data=sec_grow2, f=f, turnover=turnover, 
                      eps=extinction.fraction, hidden.states=FALSE,
                      trans.rate=trans.rates.bisse)
Null_sec_grow

# Model 2. BiSSE: turnover different; extinction fraction equal for the two states.
turnover <- c(1,2) #different turnover rate
extinction.fraction<-c(1,1) #equal extinction fraction
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0) # different transition rates

# running Model 2
BiSSE_sec_grow <- hisse(phy=tree_prunned, data=sec_grow2, f=f, turnover=turnover,
                       eps=extinction.fraction, hidden.states=FALSE,
                       trans.rate=trans.rates.bisse)
BiSSE_sec_grow


# Model 3. HiSSE: turnover different for all the states (observed and hidden).
turnover <-c(1,2,3,4)
extinction.fraction<- c(1,1,1,1)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits = 1)

# running Model 3
HiSSE_sec_grow <- hisse (phy=tree_prunned, data=sec_grow2, f=f, turnover=turnover,
                        eps=extinction.fraction, hidden.states=TRUE,
                        trans.rate=trans.rate.hisse)
HiSSE_sec_grow

# Model 4. CID2: turnover different for the hidden states and equal for the observed states,
# that is, character-independent diversification.
turnover<- c(1,1,2,2)
extinction.fraction<-c(1,1,1,1)
trans.rate<-TransMatMakerHiSSE(hidden.traits = 1, make.null = TRUE) #Three rates
trans.rate

CID2_sec_grow<-hisse(phy=tree_prunned, data=sec_grow2, f=f, turnover=turnover, eps=extinction.fraction, 
                    hidden.states=TRUE,trans.rate=trans.rate)
CID2_sec_grow


# Model 5. CID4: turnover different for the hidden states but there are two additional
# states (C and D), equal turnover for the observed states, it has the same number of
# parameters than the Model 3-HiSSE so this evaluates the character-independent 
# model as a null model to compare the HiSSE model.

turnover<-c(1,1,2,2,3,3,4,4)
extinction.fraction<-rep(1,8)
trans.rate<-TransMatMakerHiSSE(hidden.traits = 3, make.null = TRUE)
trans.rate

CID4_sec_grow<-hisse(phy=tree_prunned, data=sec_grow2, f=f, turnover=turnover, eps=extinction.fraction, 
                    hidden.states=TRUE,trans.rate=trans.rate)
CID4_sec_grow



# Summarize the results of the models to compare the likelihoods and select models
models_sec_grow<- as.data.frame(rbind(Null=c(Null_sec_grow$loglik, Null_sec_grow$AICc),
                                     BiSSE=c(BiSSE_sec_grow$loglik, BiSSE_sec_grow$AICc),
                                     HiSSE=c(HiSSE_sec_grow$loglik, HiSSE_sec_grow$AICc),
                                     CID2=c(CID2_sec_grow$loglik, CID2_sec_grow$AICc),
                                     CID4=c(CID4_sec_grow$loglik, CID4_sec_grow$AICc)))
                      
colnames(models_sec_grow) <- c("loglik", "AICc")                             
# order models in the dataframe
models_sec_grow<-models_sec_grow[order(models_sec_grow$AICc), ]
models_sec_grow
#       loglik     AICc
#CID2  -311.4506 635.7949
#CID4  -311.0618 639.6888
#HiSSE -309.6405 641.7255
#BiSSE -318.6060 647.8436
#Null  -320.5251 649.4669

# obtaining Akaike values
aic <- models_sec_grow$AICc

# selecting the best model (i.e., minimum value of Akaike)
best <- min(aic)

# obtaining the difference between Akaike values
delta <- aic - best

#obtaining Akaike weights by first calculating the relative likelihood of the models using delta
sumDelta <- sum(exp(-0.5 * delta))

# The Akaike weight for a model is this value divided by the sum of these values across all models.
# http://brianomeara.info/aic.html
w_models_sec_grow <- round((exp(-0.5 * delta)/sumDelta), 3) %>%
  as.matrix()
colnames(w_models_sec_grow) <- c("AICw")
w_models_sec_grow

# saving into a matrix
models2_sec_grow <- cbind(models_sec_grow, w_models_sec_grow) %>%
  as.data.frame() %>%
  write.csv("hisse_models_sec_grow.csv")

#########

# plot results from the best models

# The selected models are:
# Habit - CID2 and HiSSE
# Eustele type - BiSSE
# Secondary growth - CID2 and CID4
# Conclusion: only ustele type shows an association with diversification rate.


# Reconstruction of ancestral states given the parameter estimates obtained with the hisse
# function.

# MarginReconHiSSE function returns an object of class hisse.states
# that will be the input for the model plot.hisse.states function.

recon_BiSSE_eustele<-MarginReconHiSSE(phy = tree, data = eustele, f<-c(0.07,0.21),
                                    pars = BiSSE_eustele$solution,AIC = BiSSE_eustele$AIC)

# Plot
plot.hisse.states(recon_BiSSE_eustele, rate.param="net.div",rate.colors=c("yellow","purple"),
                  state.colors=c("white", "grey50"),
                  edge.width.rate=5, edge.width.state=2,type="phylogram",
                  rate.range=NULL, show.tip.label=TRUE, fsize=0.4,
                  lims.percentage.correction=0.001,legend="tips",
                  legend.cex=0.7,legend.position=c(0, 0.2, 0, 0.2),
                  legend.kernel.rates="gaussian",
                  legend.kernel.states="auto",legend.bg="gray88")


