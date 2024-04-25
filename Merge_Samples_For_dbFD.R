library(FD)
library(tidyverse)
# This small LakesMergedLakes# This small script is just to put here the function to transform each one of the selected communities during the NATs 
# and prepare the "mega-communities" for the dbFD so basically we will be summing and adding all the communities that 
# we used in the NATS into a "big" sample table.

# In our example we had:
# 21 years
# A sequence of 1 11 21 31 41 51 ponds being selected
# 10 iterations 
# This makes a total of 21*6*10= 1260 combinations (AKA a table with 1200 rows)

# Each row of the LakesMergedLakes correspond to one of these combinations  

#LakesMergedLakes <- LakesMergedLakes[-1,]
# We charge the dataset and we eliminate the Lake and Year columns to make it more real to what we will have ;) 
sp_sites <- read.csv2("data/DataAnnaTest/fuzzy_traits.csv",dec = ".") #%>% mutate(Year=as.character(Year)) %>% mutate_if(is.numeric,~ifelse(.>0,1,0)) %>% mutate(Year=as.numeric(Year))
load("data/DataAnnaTest/dis_traits_lakes.RData")
# What do we want now? We want to select the ponds that are listed in each row of our table no? So we need to make a "loop"
# where for each row we will select the listed lakes, filter them from the sp_sites and sum their abundance values

# We create the output dataframe that will contain all the combinations of communities
Combis_Communities <- matrix(ncol =219,nrow = nrow(LakesMergedLakes))
for (combi in 1:nrow(LakesMergedLakes)) {
Lakes_Combin <- LakesMergedLakes[combi,which(is.na(LakesMergedLakes[combi,])==F)] # We extract The row that corresponds to a determined combination 
# Note that we select only the values that are NOT an NA

sp_sites_Combin <- sp_sites %>% filter(Year==as.numeric(LakesMergedLakes[combi,1])) %>% # Filter the sp_sites by the year 
                                filter(Lake%in%Lakes_Combin[4:length(Lakes_Combin)]) %>%  # Filter the lakes by the lakes present in the combination
                                select(-c(Year,lake_year,Lake)) # We "deselect" the three cathegorical metrics

# So, we have generated a dataframe that contains, all the samples (rows) from a year and a combination of ponds. Now we just need to 
# sum by columns and transform the result into 1/0. Finally, bind it to the output dataframe and DONE! :D

LakesMergedLakes_id<-LakesMergedLakes[combi,1:3]
#LakesMergedLakes_comb<-ifelse(apply(sp_sites_Combin,2,sum)>0,1,0) #%>% mutate_if(is.character, as.numeric)
LakesMergedLakes_comb<-apply(sp_sites_Combin,2,sum)
LakesMerged<-c(LakesMergedLakes_id,LakesMergedLakes_comb)
LakesMerged<-as.numeric(LakesMerged)
Combis_Communities[combi,] <- LakesMerged
}

setwd("C:/Users/Anna/OneDrive - Universitat de Girona/tesi/WP1 NATs/NAT_FunMetaNET")
# Each one of the rows corresponds a each one of the used combinations and can be identified with the Year, n_sites and iteration.
# Therefore we can link it with the values of the NATs
#Combis_Communities
#summary(Combis_Communities)
# In case of need we could also "merge" the three identifiers: 
sp_names <- cbind(rep("sp",(ncol(Combis_Communities)-3)),1:(ncol(Combis_Communities)-3))
colnames(Combis_Communities) <- c("Year", "n_sites", "it",paste(sp_names[,1],sp_names[,2],sep=""))
                                  
Combis_Communities<-as.data.frame(Combis_Communities) %>% mutate(Combi_ID=paste(Year,n_sites,it, sep="_"),.before = Year) %>% select(-c(Year,n_sites,it))

# PISTA!: Sempre que combinis "noms" com has fet amb el lake_year assegurat d'unir-los amb un "_" o un ".", d'aquesta 
# manera, en cas que volguessis separar-los seria fàcil perquè podries identificar per on "partir-los"
# La funció strsplit() va molt bé per fer això però necessita un identificador per on separar ;)
#apply(Combis_Communities, 2,sum)

Combis_Communities <- data.frame(Combis_Communities, row.names = 1)%>% mutate_if(is.character, as.numeric)#per posar el nom de les localitats 
#summary(Combis_Communities)
resFD <- dbFD( 
  dis_traits, 
  log(Combis_Communities+1), #fem el LOG si treballem amb abun, si fem P/A no cal
  corr = "cailliez",
  w.abun=TRUE,
  stand.FRic = TRUE, 
  m=7) #m=número de dim
important.indices <- cbind(resFD$nbsp, resFD$FRic, resFD$FEve, resFD$FDiv,
                           resFD$FDis, resFD$RaoQ)
colnames(important.indices) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")
save(important.indices,file = "res_FD_lakes_dist_def.RData")
bind_cols(LakesMergedLakes[,1:3], important.indices)%>%
  pivot_longer(cols = 4:ncol(.)) %>% 
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>%
  #filter(name=="FEve")%>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=as.numeric(n_sites),colour=as.factor(Year)), width = 0.2)+
  geom_smooth(aes(y=value, x=as.numeric(n_sites), colour=as.factor(Year)), method="loess",se=F)+
  #geom_line(aes(y=Mean_val, x=as.numeric(n_sites), colour=as.factor(Year)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()
important.indices <- rownames_to_column(as.data.frame( important.indices), "Names")
writexl::write_xlsx(as.data.frame(important.indices),path = "res_FD_PA_env.xlsx")

strsplit(df, split = "_")

#rivers####
#This small LakesMergedLakes# This small script is just to put here the function to transform each one of the selected communities during the NATs 
# and prepare the "mega-communities" for the dbFD so basically we will be summing and adding all the communities that 
# we used in the NATS into a "big" sample table.

# In our example we had:
# 21 years
# A sequence of 1 11 21 31 41 51 ponds being selected
# 10 iterations 
# This makes a total of 21*6*10= 1260 combinations (AKA a table with 1200 rows)

# Each row of the LakesMergedLakes correspond to one of these combinations  
library(FD)
library(tidyverse)
#LakesMergedLakes <- LakesMergedLakes[-1,]
# We charge the dataset and we eliminate the Lake and Year columns to make it more real to what we will have ;) 
sp_sites <- readxl::read_xlsx("data/DataAnnaTest/sp_rivers.xlsx") #%>% mutate(Year=as.character(Year)) %>% mutate_if(is.numeric,~ifelse(.>0,1,0)) %>% mutate(Year=as.numeric(Year))
load("data/DataAnnaTest/dis_traits_river.RData")
# What do we want now? We want to select the ponds that are listed in each row of our table no? So we need to make a "loop"
# where for each row we will select the listed lakes, filter them from the sp_sites and sum their abundance values

# We create the output dataframe that will contain all the combinations of communities
Combis_Communities <- matrix(ncol =381,nrow = nrow(LakesMergedLakes))
for (combi in 1:nrow(LakesMergedLakes)) {
  Lakes_Combin <- LakesMergedLakes[combi,which(is.na(LakesMergedLakes[combi,])==F)] # We extract The row that corresponds to a determined combination 
  # Note that we select only the values that are NOT an NA
  
  sp_sites_Combin <- sp_sites %>% filter(Year==as.numeric(LakesMergedLakes[combi,1])) %>% # Filter the sp_sites by the year 
    filter(Lake%in%Lakes_Combin[4:length(Lakes_Combin)]) %>%  # Filter the lakes by the lakes present in the combination
    select(-c(Year,lake_year,Lake)) # We "deselect" the three cathegorical metrics
  
  # So, we have generated a dataframe that contains, all the samples (rows) from a year and a combination of ponds. Now we just need to 
  # sum by columns and transform the result into 1/0. Finally, bind it to the output dataframe and DONE! :D
  
  LakesMergedLakes_id<-LakesMergedLakes[combi,1:3]
  #LakesMergedLakes_comb<-ifelse(apply(sp_sites_Combin,2,sum)>0,1,0) #%>% mutate_if(is.character, as.numeric)
  LakesMergedLakes_comb<-apply(sp_sites_Combin,2,sum)
  LakesMerged<-c(LakesMergedLakes_id,LakesMergedLakes_comb)
  LakesMerged<-as.numeric(LakesMerged)
  Combis_Communities[combi,] <- LakesMerged
}


# Each one of the rows corresponds a each one of the used combinations and can be identified with the Year, n_sites and iteration.
# Therefore we can link it with the values of the NATs
#Combis_Communities
#summary(Combis_Communities)
# In case of need we could also "merge" the three identifiers: 
sp_names <- cbind(rep("X",(ncol(Combis_Communities)-3)),1:(ncol(Combis_Communities)-3))
colnames(Combis_Communities) <- c("Year", "n_sites", "it",paste(sp_names[,1],sp_names[,2],sep=""))

Combis_Communities<-as.data.frame(Combis_Communities) %>% mutate(Combi_ID=paste(Year,n_sites,it, sep="_"),.before = Year) %>% select(-c(Year,n_sites,it))

# PISTA!: Sempre que combinis "noms" com has fet amb el lake_year assegurat d'unir-los amb un "_" o un ".", d'aquesta 
# manera, en cas que volguessis separar-los seria fàcil perquè podries identificar per on "partir-los"
# La funció strsplit() va molt bé per fer això però necessita un identificador per on separar ;)
#apply(Combis_Communities, 2,sum)
setwd("C:/Users/Anna/OneDrive - Universitat de Girona/tesi/WP1 NATs/NAT_FunMetaNET")
Combis_Communities <- data.frame(Combis_Communities, row.names = 1)%>% mutate_if(is.character, as.numeric)#per posar el nom de les localitats 
#summary(Combis_Communities)
resFD <- dbFD( 
  dis_traits, 
  log(Combis_Communities+1), #fem el LOG si treballem amb abun, si fem P/A no cal
  corr = "cailliez",
  w.abun=TRUE,
  stand.FRic = TRUE, 
  m=7) #m=número de dim
important.indices <- cbind(resFD$nbsp, resFD$FRic, resFD$FEve, resFD$FDiv,
                           resFD$FDis, resFD$RaoQ)
colnames(important.indices) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")
save(important.indices,file = "res_FD_rivers_DisDist_def.RData")
bind_cols(LakesMergedLakes[,1:3], important.indices)%>%
  pivot_longer(cols = 4:8) %>% 
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>%
  filter(Year==c("2008", "2010", "2017"))%>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=as.numeric(n_sites),colour=as.factor(Year)), width = 0.2)+
  geom_smooth(aes(y=value, x=as.numeric(n_sites), colour=as.factor(Year)), method="loess",se=F)+
  #geom_line(aes(y=Mean_val, x=as.numeric(n_sites), colour=as.factor(Year)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()
important.indices <- rownames_to_column(as.data.frame( important.indices), "Names")

writexl::write_xlsx(as.data.frame(important.indices),path = "res_FD_abun_env.xlsx")
data <- readxl::read_excel("res_FD_abun_env.xlsx")
data2 <- separate(data, Names, into = c("Year", "n_sites", "iter"), sep = "_")

data2 %>%
  pivot_longer(cols = 4:8) %>% 
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>%
  filter(Year==c("1995", "2000", "2010"))%>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=as.numeric(n_sites),colour=as.factor(Year)), width = 0.2)+
  geom_smooth(aes(y=value, x=as.numeric(n_sites), colour=as.factor(Year)), method="loess",se=F)+
  #geom_line(aes(y=Mean_val, x=as.numeric(n_sites), colour=as.factor(Year)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()
