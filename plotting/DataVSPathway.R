setwd("~/Uni/Semester08/KU/Harvard")

# =================================================
# Looking at the data set
# =================================================

data <- t(as.matrix(read.table("mmc4tobias.CSV",skip = 3, sep = ";")))

colnames(data) <- data[1,]
data <- data[-1,]

rownames(data) <- data[,1]
data <- data[,-1]

# only looking at cell lines
data_CL <- data[data[,2] == "cell line",]
data_CL <- data_CL[,-2]
data_CL <- data_CL[,-2]

# =================================================
# Loading Pathway activation
# =================================================

PW <- as.matrix(read.table("Pathways.csv",header = TRUE, sep = ","))
rownames(PW) <- PW[,1]
PW <- PW[,-1]

# =================================================
# Joining data and pathway activity
# =================================================

data_join <- cbind(PW, data_CL[match(rownames(PW), rownames(data_CL)),])

# =================================================
# Loading IC50-value
# =================================================

IC50 <- as.matrix(read.table("IC50.csv", skip = 4, header = TRUE, sep = ","))

colnames(IC50) <- IC50[1,]
IC50 <- IC50[-1,]
IC50 <- IC50[-1,]

rownames(IC50) <- IC50[,2]
IC50 <- IC50[,-1]
IC50 <- IC50[,-1]

IC50_reduced <- IC50[match(rownames(data_join), rownames(IC50)),]


# =================================================
# Summary on IC50
# =================================================

IC50_all <-sum(IC50_reduced != 0)
IC50_R <- sum(IC50_reduced == "R")
IC50_S <- sum(IC50_reduced == "S")
IC50_BL <- IC50_all - IC50_R - IC50_S

100*(IC50_R/IC50_all)
100*(IC50_S/IC50_all)
100*(IC50_BL/IC50_all)

# =================================================
# Finding drug target
# =================================================

setwd("~/Uni/Semester08/KU/Harvard/Drug")

data <- read.table("ScreenedCompounds.csv",
                   skip = 2,
                   header = TRUE,
                   na.strings = TRUE,
                   sep = ";")

data[,1] <- NULL
data[,3] <- NULL
data[,4] <- NULL

data < -as.data.frame(data)

tmp <- (strsplit(as.character(data$Putative.Target), "[,]"))

Putative_target <- as.data.frame(matrix(NA, ncol = 9, nrow = length(tmp)))

for(i in c(1:length(tmp))){
  for(j in 1:9){
    
    Putative_target[i,j] <- tmp[[i]][j]
    
  }
}

Putative_target_m <- as.matrix(Putative_target)
Putative_target_m <-gsub(" ", "", Putative_target_m)

a <- unclass(factor(Putative_target_m))

dim(a) <- c(length(tmp), 9)

Putative_zeros = matrix(0, ncol = 212, nrow = length(tmp))

for(i in c(1:length(tmp))){
  for(j in 1:9){
    Putative_zeros[i,a[i,j]] <- 1
  }
}

colnames(Putative_zeros) <- levels(a)
rownames(Putative_zeros) <- data$Name

drug_target <- t(Putative_zeros)

# =================================================
# Finding drug target
# =================================================

IC50_reduced_r <- IC50_reduced[,match(colnames(drug_target), colnames(IC50_reduced))]
drug_target_r <- drug_target[,match(colnames(drug_target), colnames(IC50_reduced))]

IC50_reduced_r <- IC50_reduced_r[,!is.na(colnames(drug_target_r) == colnames(IC50_reduced_r))]
drug_target_r <- drug_target_r[,!is.na(colnames(drug_target_r) == colnames(IC50_reduced_r))]

# =================================================
# 
# =================================================


# =================================================
# 
# =================================================


# =================================================
# 
# =================================================

cancer_type <- data_join[,12]
data_full <- data_join[,c(1:11,13:ncol(data_join))]

cancer_type <- as.numeric(as.factor(cancer_type))

IC50_full <- IC50_reduced_r

for(i in c(1:ncol(IC50_reduced_r))){
  
  IC50_full[,i] <- as.numeric(as.factor(IC50_reduced_r[,i]))    
  
}

drug_target_full <- drug_target_r


length(cancer_type)
ncol(data_full)
length(IC50_full)
ncol(drug_target_full)

setwd("~/Uni/Semester08/KU/Harvard/Model")

write.csv(cancer_type, "cancer_type.csv")
write.csv(data_full, "data.csv")
write.csv(drug_target_full, "drug_target.csv")
write.csv(IC50_full, "IC50.csv")



