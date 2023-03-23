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

pvals = c(rep(0, 16529-4))

"
- iteriere ueber alle Gene
- erstelle data frame mit designmatrix und readcounts des aktuellen Gens
- currentglm: verallgemeinertes lineares Modell mit glm.nb fuer aktuelles Gen
"
for (i in 4:16529)
{
  currentcounts = arab.data[i]
  currentvals = data.frame(designmatrix, currentcounts)
  
  colnames(currentvals)[3] = "counts"
  
  currentglm = glm.nb(counts ~ 1 + treatment + as.factor(time), data = currentvals)
  
}

#statistische Tests (p-value)

####Simulation
#Zeit ignorieren, nur noch mock oder hrcc
#Annahme: gleiche Expression von 'unwichtigen Genen', Unterschiede für DE Gene
#Einfluss Anzahl Wiederholungen auf die FDR
