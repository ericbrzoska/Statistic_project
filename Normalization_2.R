# normalization via DESeq2: estimate size factors

data= t(read.csv2("arab.csv", sep= ","))
copy=data
colnames(data)=data[1,]
data=data[-c(1,2,3),]
num_data= apply(data, 2 ,as.numeric)
colnames(num_data)=colnames(data)
row.names(num_data)=row.names(data)
num_data=num_data+1
#min_count= length(num_data[1,])*30
#num_data = num_data[which(rowSums(num_data)>min_count),]


#create an "average" sample the represents the counts that are typically sample of a dataset
pseudo_ref = sqrt(rowProds(num_data))

#we assuming, that the genes are not changing dramatically, so we calculate rations for every gene
#and define size factors as median of those ratios for each sample

test= as.double(num_data)/as.double(pseudo_ref)
dim(test)=dim(num_data)
norm_factors = colMedians(test)

num_data[,1]= num_data[,1]/norm_factors[1]
num_data[,2]= num_data[,2]/norm_factors[2]
num_data[,3]= num_data[,3]/norm_factors[3]
num_data[,4]= num_data[,4]/norm_factors[4]
num_data[,5]= num_data[,5]/norm_factors[5]
num_data[,6]= num_data[,6]/norm_factors[6]

norm_counts= num_data %/% 1 
output=t(rbind(copy[c(1:3),],norm_counts))
output=output[,-1]

output= cbind(row.names(output),output)
colnames(output)[1]="X"
row.names(output)=c("V1","V2","V3","V4","V5","V6")


write.csv2(output,file = "normed_arab_2.csv")
