#imports
require(MASS)

#Datei einlesen
data = read.csv2("arab.csv", sep=",") 
#data = read.csv2("normed_arab.csv", sep=";") 

"
verallgemeinertes lineares Modell, p-values

- init Einflussgroessen treatment, time als data frame designmatrix
- init p-value-Vektor
"
treatment = c(rep(0,3), rep(1,3))
time = c(rep(1:3,2))

designmatrix = data.frame(treatment, time)

pvals = c(rep(0, 16529-4))

"
- iteriere ueber alle Gene
- erstelle data frame mit designmatrix und readcounts des aktuellen Gens
- currentglm: verallgemeinertes lineares Modell mit glm.nb fuer aktuelles Gen
"
for (i in 4:length(data[1,]))
{
  currentcounts = data[i]
  currentvals = data.frame(designmatrix, currentcounts)
  
  colnames(currentvals)[3] = "counts"
  
  currentglm = glm(counts ~ 1 + treatment + as.factor(time), data = currentvals, family=poisson)
  
  pvalue = summary(aov(currentglm))[[1]][2,5]
  
  pvals[i-3] = pvalue
  
}

"
- Einschub Residuenplots für ausgewählte Gene
"
currentvals = data.frame(designmatrix, data[1243])
currentvals
colnames(currentvals)[3] = "counts"

nbmodel = glm.nb(counts ~ 1 + treatment + as.factor(time), data = currentvals)

poissmodel = glm(counts ~ 1 + treatment + as.factor(time), data = currentvals, family=poisson)

p_res <- resid(poissmodel)
nb_res = resid(nbmodel)

plot(fitted(poissmodel), p_res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Poisson')

plot(fitted(nbmodel), nb_res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Negative Binomial')

"
Statistische Tests
- initialisiere DEgenes

- iteriere über alle p-values
- sortiere alle groesser 0.05 aus
- fuer alle kleiner 0.05 cbind zu neuem dataframe und berechne fold change
- falls fold change größer 2 oder kleiner 0.5, Gen ist DE Gen
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

#Histogramm für die p-values
hist(pvals,breaks=seq(from=0, to=1, by=0.05),main="Histogram of p-values",xlab="p-values",ylim=c(0,3000))

####Simulation
#Zeit ignorieren, nur noch mock oder hrcc
#Annahme: gleiche Expression von 'unwichtigen Genen', Unterschiede für DE Gene
#Einfluss Anzahl Wiederholungen auf die FDR
