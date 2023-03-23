"Imports
-> if you run the code the first time, uncommend the installation lines
"
#install.packages("here")
require(MASS)
library(here)

"
Initiation
-> Shall make the import easier, so no one needs to change this path
-> pls save the Github_folder in your Download directory
"

init= getwd()
setwd(here("Downloads","Statistic_project-main"))     

"
Coding space
-> it would be great to
"

data = t(read.csv2("arab.csv", sep=","))              # Reads .csv in and transposes it by t()

"
-> Switching work directory back to normal
"
setwd(init)

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
