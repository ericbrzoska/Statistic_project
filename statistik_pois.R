#imports
require(MASS)

#Datei einlesen
data = read.csv2("arab.csv", sep=",") 

"
verallgemeinertes lineares Modell

- init Einflussgroessen treatment, time als data frame designmatrix
- init pvalue-Vektor
"
treatment = c(rep(0,3), rep(1,3))
time = c(rep(1:3,2))

designmatrix = data.frame(treatment, time)

pvals = c(rep(0, length(data[1,])))

"
- iteriere ueber alle Gene
- erstelle data frame mit designmatrix und readcounts des aktuellen Gens
- currentglm: verallgemeinertes lineares Modell mit glm.nb fuer aktuelles Gen
"
for (i in 4:length(data[1,])){
  
  currentvals = data.frame(designmatrix, data[i])
  
  colnames(currentvals)[3] = "counts"

  currentglm = glm(counts ~ 1 + treatment + as.factor(time), data = currentvals)
  
  pvalue = summary(aov(currentglm))[[1]][2,5]
  
  pvals[i-3] = pvalue
  
}

"
- initialisiere DEgenes und NDEgenes

- iteriere über alle p-values
- sortiere alle groesser 0.05 aus
- fuer alle kleiner 0.05 adde zu neuem dataframe und berechne fold change
"
DEgenes = data.frame(data[1], data[2], data[3])
NDEgenes = data.frame(data[1], data[2], data[3])

for(i in 4:(length(data[1,]))){
  
  if(pvals[i-3] < 0.05){
    
    mean.mock = mean(cbind(data[1, i], data[2, i], data[3, i]))
    mean.hrcc = mean(cbind(data[4, i], data[5, i], data[6, i]))

    if(mean.mock/mean.hrcc > 2 | mean.mock/mean.hrcc < 0.5){
      
      DEgenes = cbind(DEgenes, data[i])
      
    }
    
  } else{
    
    NDEgenes = cbind(NDEgenes, data[i])
    
  }
  
}
