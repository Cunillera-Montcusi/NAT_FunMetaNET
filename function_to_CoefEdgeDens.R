

# https://www.statforbiology.com/2020/stat_nls_usefulfunctions/

install.packages("devtools")
devtools::install_github("onofriAndreaPG/aomisc")

library(drc)
library(nlme)
library(statforbiology)

X <- c(1, 3, 5, 7, 9, 11, 13, 20)
Y <- c(8.22, 14.0, 17.2, 16.9, 19.2, 19.6, 19.4, 19.6)

# nls fit
model <- nls(Y ~ NLS.asymReg(X, init, m, plateau) )

# Equits és igual a les iteracions
# LaI és igual al valor d'edge density
Equits <- c(1, 3, 5, 7, 9, 11, 13, 20)
LaI <- c(8.22, 14.0, 17.2, 16.9, 19.2, 19.6, 19.4, 19.6)
plot(Equits,LaI)

model_S <- nls(LaI~NLS.asymReg(Equits, init, m, plateau))

p<-coefficients(model_S)
p
plot(predict(model_S))

     