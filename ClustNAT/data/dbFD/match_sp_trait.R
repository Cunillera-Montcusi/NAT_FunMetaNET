
library(tidyverse)# managing data.frames and organize and edit them
library(FD)#càlcul de manera ràpida les mètriques de diversitat funcional
library(corrplot)#per fer gràfics de correlacions de matrius
library(readxl)#per carregar arxius d'excel
library(gawdis)#calcular distàncies ponderant pel pes dels traits i utilitzant fuzzy coding
library(ade4)#per fer pcoa
library(dplyr)#package per poder transformar les dades
library(mFD)
library("writexl")
library(dplyr)

#Fem la mitjana dels valors dels traits dels gèneres de la base de dades de traits que estaven repetits.
traits <- traits %>% group_by(Genus..if.description.at.this.level.) %>% 
  mutate_if(is.numeric, ~mean(.)) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)
#amb les dades resultats fem un excel per convertir a fuzzy i fer el LOG+1 (dades a fuzzy_traits.xlsx)

#Calculem la diversitat funcional: amb dades d'abundància d'sp(en matriu), dades de traits amb els variables passades a fuzzy
# necessitem fer un join entre les dades de traits_fuzzy i les dades d'sp, ho fem a continuació

# Carreguem la base de dades de traits amb les categories passades a fuzzy
# i fet el LOG+1
traits_fuz <-read_excel("C:/Users/anna.equisuany/OneDrive - Universitat de Girona/tesi/WP1 NATs/calcul de div funcional/fuzzy_traits.xlsx",sheet = "traits_fuzzy_log")
taxa_list <-read_excel("C:/Users/anna.equisuany/OneDrive - Universitat de Girona/tesi/WP1 NATs/calcul de div funcional/fuzzy_traits.xlsx",sheet = "taxa_list")
sp_abun <-read_excel("/Users/anna.equisuany/OneDrive - Universitat de Girona/tesi/WP1 NATs/calcul de div funcional/fuzzy_traits.xlsx",sheet = "sp_lakes") 
sp_abun$Year <- as.character(sp_abun$Year)

sp_pa<-sp_abun %>%  mutate_if(is.numeric, ~ifelse(.>0,1,0)) #convertim a P/A


#unim del df de traits fuzzy i taxa-list per obtenir les dades de traits per les nostres sp.
traits_ind <- taxa_list %>%
  ungroup() %>% # just in case 
  left_join(traits_fuz, by=c("genus"="Modalities"),multiple ="all")%>% #join traits with the genus of the traits database
  select(c(4:ncol(.))) %>% # select the columns with only traits
  na.omit()


 traits_num <- as.data.frame(apply(traits_ind, 2, as.numeric)) #passem totes les variables com a numèriques
 rownames(traits_num)<- c(taxa_list$id) #tornem a posar id com a rownames
#colnames(traits1)<-c(traits_cat$trait_name)


# ja tenim els tres fitxers que necessitem per calcular la div funcional:
   #1. Els traits de les espècies dels llacs en format fuzzy y LOG+1(amb els llocs i els anys) --> traits_ind
   #2. la matriu d'espècies de llacs (amb els llocs i els anys) -->sp_abun 
   #3. Un arxiu on indiquem les categories dels traits: El carreguem a continuació --> traits_cat
traits_cat <- read_excel("C:/Users/anna.equisuany/OneDrive - Universitat de Girona/tesi/WP1 NATs/calcul de div funcional/fuzzy_traits.xlsx", sheet = "traits_cat")

#Ara creem un sistema iteratiu que vagi fent el càlcul de diversitat funcional per cada any i cada llac


rownames(traits_num)<- c(traits$Modalities)
# Càlcul de la matriu de distancies funcionals
dis_traits<-gawdis(traits_num, w.type = "optimized", 
                   groups = c(1, 1, 1, 1, 1, 1, 1, 2, 2,
                              3, 3, 3, 4, 4, 4, 4,
                              5, 5, 5, 5, 5, 5, 5, 5, 
                              6, 6, 6, 6, 7, 7, 7, 7, 7, 
                              8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 
                              10, 10, 10, 10, 10, 10, 10, 10, 10,
                              11, 11, 11, 11, 11, 11, 11, 11, 
                              12, 12, 12, 12, 12, 12, 12,
                              13, 13, 13, 13, 13, 13, 13, 13,
                              14, 14, 14, 15, 15, 15, 15, 15, 
                              16, 16, 16, 16, 16, 16, 16, 16, 16,
                              17, 17, 17, 17, 18, 18, 18,
                              19, 19, 20, 20, 20, 21, 21, 21, 21, 21, 
                              22, 22, 22, 22, 22, 22
                   ), fuzzy = TRUE)

save(dis_traits, file = "dis_traits_lakes.RData")#guardem la matriu de distancies
load("dis_traits_lakes.RData") 
#Per comprovar quin triat té més pes en la mesura de distancies funcionals finals (*correls)* 
#i quin és el weight aplicat en la mitjana ponderada de les matrius de distàncies (*weights*)
attr(dis_traits,"correls")#comprovar que es correlacions son similar entre traits, per tant el weight ha funcionat
attr(dis_traits,"weights")#comproves el valor de weight que ha aplicat


# Triar nombre de dimensions de l'espai funcional--> PCoA de la matriu de distancies (dis_traits)
fspaces_quality_fruits <- mFD::quality.fspaces(
  sp_dist             = dis_traits,
  fdendro             = "average",
  maxdim_pcoa         = 20,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE)

round(fspaces_quality_fruits$"quality_fspaces", 3) #el valor minim de mad són les dimensions a agafar

#ara fem uns gràfic per comprovar el resultat. Dibuixem el resultat del dendograma promig, 
#pcoa amb 6 dimensions (com sugereix un scree plot, fet amb dudi.pco), 
#amb 7 dimensions com suggereix la comanda anterior i 
#amb 20d que seria el tope que hem posat abans i 
#on ja s'observa com empitjora la qualitat perquè els valors de qualitat tornen a augmentar.

mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fruits,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_7d","pcoa_20d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

PCOA <- dudi.pco(dis_traits, scannf = FALSE, nf = 5) #sceeplot per triar dimensions (barra en negre i gris)
scatter(PCOA, xax = 1, yax = 2, clab.row = 1, posieig = "top", sub = NULL, csub = 2) 
#script per calcular els paràmetres: loop a partir d'aquí
sp_pa <- sp_pa %>% filter(Year=="1995") #filtrem any
sp_pa <-select(sp_pa, c(3:219)) 
sp_pa <- data.frame(sp_pa, row.names = 1)#per posar el nom de les localitats 

#càlcul mètriques alfa div
resFD <- dbFD( 
  dis_traits, 
  sp_pa, #fem el LOG si treballem amb abun, si fem P/A no cal
  corr = "cailliez",
  w.abun=TRUE,
  stand.FRic = TRUE, 
  m=7) #m=número de dim

#guargar els resultats
important.indices <- cbind(resFD$nbsp, resFD$FRic, resFD$FEve, resFD$FDiv,
                           resFD$FDis, resFD$RaoQ)
colnames(important.indices) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")
important.indices<-as.data.frame(important.indices)
important.indices<-tibble::rownames_to_column(important.indices, "rownames")
save(important.indices, file = "div_func_lakes.RData")
write_xlsx(important.indices, "C:/Users/ecologia/OneDrive - Universitat de Girona/tesi/WP1 NATs/calcul de div funcional/div_func_lakes.xlsx")
