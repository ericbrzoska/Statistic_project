#imports
require(MASS)

#Datei einlesen
data = read.csv2("arab.csv", sep=",") 


#verallgemeinertes lineares Modell, p-values

treatment = c(rep(0,3), rep(1,3)) #treatment Vektor
time = c(rep(1:3,2))  #time Vektor

designmatrix = data.frame(treatment, time)  #designmatrix mit Einflussgroessen

pvals = c(rep(0, length(data[1,]))) #speichert p-value fuer jedes Gen

fittedvals = c()  #speichert fitted values fuer jedes Gen
bool.de.nde = c() #speichert boolean (TRUE wenn DE Gen, false wenn NDE Gen)


"
- iteriere ueber alle Gene
- erstelle data frame mit designmatrix und readcounts des aktuellen Gens
- currentglm: verallgemeinertes lineares Modell mit glm(...,family=poisson) fuer aktuelles Gen
- berechne p-values mit aov() un speichere fitted values
"
for (i in 4:length(data[1,])){
  
  currentvals = data.frame(designmatrix, data[i])
  
  colnames(currentvals)[3] = "counts"
  
  currentglm = glm(counts ~ 1 + treatment + as.factor(time), data = currentvals, family=poisson)
  
  pvalue = summary(aov(currentglm))[[1]][2,5]
  fittedvals = cbind(fittedvals, fitted(currentglm))
  
  "
  - pruefe, ob p-value unter threshold liegt
  - wenn ja, fold change pruefen
  - wenn ja, DE Gen
  "
  if(pvalue < 0.05){
    
    mean.mock = mean(cbind(data[1, i], data[2, i], data[3, i]))
    mean.hrcc = mean(cbind(data[4, i], data[5, i], data[6, i]))
    
    if(mean.mock/mean.hrcc > 2 | mean.mock/mean.hrcc < 0.5){
      
      bool.de.nde = cbind(bool.de.nde, TRUE)
      
    }else{
      
      bool.de.nde = cbind(bool.de.nde, FALSE)
    }
  } else{
    
    bool.de.nde = cbind(bool.de.nde, FALSE)
  }
}

#Simulation der Werte

sim.values = c()  #speichert simulierte Werte fuer jedes Gen

"
- iteriere ueber Gene
- pruefe  NDE und DE Gen mittels bool.de.nde
- simuliere Daten (NDE mit selbem Mittelwert, DE mit unterschiedlichen Werten mock/hrcc)
"
for (i in 4:length(data)) {
  
  if(bool.de.nde[i-3]){
    
    DEmean.mock = mean(fittedvals[1:3,i-3])
    DEmean.hrcc = mean(fittedvals[4:6,i-3])
    
    randpois.mock = rpois(250, DEmean.mock)
    randpois.hrcc = rpois(250, DEmean.hrcc)
    
    sim.values = cbind(sim.values, c(randpois.mock, randpois.hrcc))
    
    colnames(sim.values)[i-3] = colnames(data[i])
    
  }else{
    
    NDEmean = mean(NDEfittedvals[i-3])
    
    randpois = rpois(500, NDEmean)
    
    sim.values = cbind(sim.values, randpois)
    
    colnames(sim.values)[i-3] = colnames(data[i])
  }
  
}

