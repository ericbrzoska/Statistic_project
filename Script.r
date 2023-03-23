require(MASS)


data = read.csv2("arab.csv", sep=",")


#verallgemeinertes lineares Modell
treatment = c(rep(0,3), rep(1,3))
time = c(rep(1:3,2))

designmatrix = data.frame(treatment, time)

yis = c(rep(0, 16529-4))
for (i in 4:length(16529))
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
