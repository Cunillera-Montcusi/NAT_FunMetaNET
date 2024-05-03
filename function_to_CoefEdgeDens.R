

# https://www.statforbiology.com/2020/stat_nls_usefulfunctions/

install.packages("devtools")
devtools::install_github("onofriAndreaPG/aomisc")

library(aomisc)

# Equits és igual a les iteracions
# LaI és igual al valor d'edge density

plot(Equits,LaI)

model_S <- nls(LaI~NLS.asymReg(Equits, init, m, plateau))

p<-coefficients(model_S)
p
plot(predict(model_S))

     