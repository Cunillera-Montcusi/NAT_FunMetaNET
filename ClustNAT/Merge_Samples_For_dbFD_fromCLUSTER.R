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

# Lakes
sp_sites <- read.csv2("ClustNAT/data/dbFD/fuzzy_traits.csv",dec = ".") #%>% mutate(Year=as.character(Year)) %>% mutate_if(is.numeric,~ifelse(.>0,1,0)) %>% mutate(Year=as.numeric(Year))
load(paste(getwd(),"/data/dbFD/dis_traits_lake.RData",sep=""))

# Rivers
sp_sites <- readxl::read_excel("ClustNAT/data/dbFD/sp_rivers.xlsx")
colnames(sp_sites)[4:length(colnames(sp_sites))] <- paste("X",colnames(sp_sites)[4:length(colnames(sp_sites))],sep="")
load("ClustNAT/data/dbFD/dis_traits_river.RData")

load("ClustNAT/ClusterNATs.RData")

# Remove NAs from all the rows
LakesMergedLakes <- Final_Output$river$LakMergLak
out_check <- c()
for(Files_NA in 1:nrow(LakesMergedLakes)){out_check[Files_NA] <- length(which(is.na(LakesMergedLakes[Files_NA,])==T))}
To_Remove <- which(out_check==ncol(Final_Output$river$LakMergLak))
LakesMergedLakes <- LakesMergedLakes[-To_Remove,]

# What do we want now? We want to select the ponds that are listed in each row of our table no? So we need to make a "loop"
# where for each row we will select the listed lakes, filter them from the sp_sites and sum their abundance values

# We create the output dataframe that will contain all the combinations of communities
Combis_Communities <- matrix(ncol =ncol(sp_sites),nrow = nrow(LakesMergedLakes))
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

# Lakes is "sp" instead of "X"
sp_names <- cbind(rep("X",(ncol(Combis_Communities)-3)),1:(ncol(Combis_Communities)-3))
colnames(Combis_Communities) <- c("Year", "n_sites", "it",paste(sp_names[,1],sp_names[,2],sep=""))
                                  
Combis_Communities<-as.data.frame(Combis_Communities) %>% mutate(Combi_ID=paste(Year,n_sites,it, sep="_"),.before = Year) %>% select(-c(Year,n_sites,it))

# We identify where the sites are repeated again (Each additive process: "Random","Environment","Distance")
Equal_Ones <- which(Combis_Communities$Combi_ID==Combis_Communities$Combi_ID[1])
Combis_Communities_TYPE_Nats <- list()

Combis_Communities_TYPE_Nats[[1]] <- Combis_Communities[Equal_Ones[1]:(Equal_Ones[2]-1),]
Combis_Communities_TYPE_Nats[[2]] <- Combis_Communities[Equal_Ones[2]:(Equal_Ones[3]-1),]
Combis_Communities_TYPE_Nats[[3]] <- Combis_Communities[Equal_Ones[3]:length(Combis_Communities$Combi_ID),]
names(Combis_Communities_TYPE_Nats) <- c("Random","Environment","Distance")

DF_FunIndices <- data.frame()
for (Type_NATs in 1:length(Equal_Ones)) {
  Combis_Communities_small <- data.frame(Combis_Communities_TYPE_Nats[[Type_NATs]], row.names = 1)%>% 
                              mutate_if(is.character, as.numeric)#per posar el nom de les localitats 
  
  resFD <- dbFD( 
    dis_traits, 
    log(Combis_Communities_small+1), #fem el LOG si treballem amb abun, si fem P/A no cal
    corr = "cailliez",
    w.abun=TRUE,
    stand.FRic = TRUE, 
    m=7) #m=número de dim
  
  important.indices <- data.frame("Type_NATs"=names(Combis_Communities_TYPE_Nats)[Type_NATs],
                                  "Combi_ID"=Combis_Communities_TYPE_Nats[[1]]$Combi_ID,
                                  "NumbSpecies"=resFD$nbsp,
                                  "FRic"=resFD$FRic,"FEve"=resFD$FEve,"FDiv"=resFD$FDiv,"FDis"=resFD$FDis,"Rao"=resFD$RaoQ)
  
  #Savecopy
  save.image(file = "functionalextraction_safecopy.RData")
  #Final output storage
  DF_FunIndices <- bind_rows(DF_FunIndices,important.indices)
}
DF_FunIndices <- DF_FunIndices %>% mutate(System="river",.before = Type_NATs)

DF_FunIndices%>%
  separate(Combi_ID,c("Year","n_sites","it")) %>% 
  pivot_longer(cols = 6:ncol(.)) %>%
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>%
  #filter(name=="FEve")%>% 
  ggplot()+ 
  #geom_jitter(aes(y=value, x=as.numeric(n_sites),colour=as.factor(Year)), width = 0.2)+
  geom_smooth(aes(y=value, x=as.numeric(n_sites), colour=as.factor(Year)), method="loess",se=F)+
  #geom_line(aes(y=Mean_val, x=as.numeric(n_sites), colour=as.factor(Year)))+
  facet_grid(name~Type_NATs,scales = "free") + 
  theme_classic()

write.csv2(DF_FunIndices,file="lakes_FunctionalIndices_River.csv")

DF_FunIndices_River <- read.csv2("lakes_FunctionalIndices_River.csv")
DF_FunIndices_Lake <- read.csv2("lakes_FunctionalIndices_Lake.csv")

DF_FunIndices <- bind_rows(DF_FunIndices_Lake,DF_FunIndices_River)
write.csv2(DF_FunIndices,file="FunctionalIndices_Tot.csv")



