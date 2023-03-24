# Genutzt um normalisierte Counting matrix zu generieren. Es ist möglich,
# dass für die Reprodzierbarkeit die R Version für dieses Skript in 4.2.2 geändert werden muss
# R 4.2.2 
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")

data= t(read.csv2("arab.csv", sep=","))
copy= data
design_matrix=t(data[c(1:3),])   # create a designmatrix 
counts_matrix=data[-c(1:3),]  # remove additional information to create a counting matrix
counts_matrix=apply(counts_matrix, 2 ,as.numeric)
row.names(counts_matrix)=row.names(copy[-c(1:3),])
colnames(counts_matrix)=data["X",]
row.names(design_matrix)= design_matrix[,1]
dds = DESeqDataSetFromMatrix(countData = counts_matrix, colData = design_matrix, design = ~ treatment + time + time:treatment)

norm_counts=estimateSizeFactors(dds)
norm_counts=counts(norm_counts, normalized= TRUE)

#reconstructing original dataframe
output=rbind(copy[c(1:3),],norm_counts)
write.csv2(output,file = "normed_arab.csv", sep = ",")
