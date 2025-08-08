fitch_path <- "/Users/mful870/Bioinformatics_Programs/phylip-3.695/exe"

Superpose_file <- "Superpose.score"

library(ape)
library(phangorn)
library(MASS)
library(Matrix)
library(reshape2)
library(stringr)
library(Rphylip)


#Input columns are: spaces - [1]'#S#' -- [2]pdb A -- [3]pdb_B -- [4]Q-dist value -- [5]RMSD -- [6]Alnign_length
#Split on multispaces

#Create the genA, genB, Q-dist from:
# Cols 3,4,5 from original data

# Read in table with DDH info
n <- read.table(file = Superpose_file,sep="",header=FALSE,stringsAsFactors=FALSE, comment.char="")

# Rename the relevant column names
colnames(n)[c(2,3,4)] <- c("Subj","Query","Value")

# Save into two objects for independant manipulation - inefficient, but works
val_matrix <- n

#Cull unneeded columns (i.e. from Equation 1 & 3 and the unneeded cols of Equation 2)
val_matrix <- val_matrix[,-c(1,5:7)]	#Leaves cols[,3,4,5] 

#reshape into a triangle matrix
val <- as.matrix(xtabs(as.numeric(Value) ~ Subj + Query, data = val_matrix))

n_row <- nrow(val)	#How many rows have the number of Q-dist values (i.e. number of taxa)
n_col <- ncol(val)	#Ditto columns

max <- 1
dist_matrix <-(1 - val)	#Turn the DDH values into a distance
write.table(dist_matrix, file = "SuperposeDistMatrix.txt", quote = FALSE, sep = "\t")
dist_matrix <- as.dist(dist_matrix)

#RFitch.FM <- Rfitch(dist_matrix, path = fitch_path, quiet = TRUE, method = "FM", negative = FALSE)
#RFitch.LS <- Rfitch(dist_matrix, path = fitch_path, quiet = TRUE, method = "LS", negative = FALSE)
RFitch.ME <- Rfitch(dist_matrix, path = fitch_path, quiet = TRUE, method = "ME", negative = FALSE)

#write.tree(RFitch.FM, file = "FM_BestTree.tre")
#write.tree(RFitch.LS, file = "LS_BestTree.tre")
write.tree(RFitch.ME, file = "ME_BestTree.tre")

#plot(ME_best)
#plot(NJ_best)
#plot(BioNJ_best)
