"Imports
-> if you run the code the first time, uncommend the installation lines
"
#install.packages("here")

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
