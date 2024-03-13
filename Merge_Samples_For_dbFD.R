

# This small script is just to put here the function to transform each one of the selected communities during the NATs 
# and prepare the "mega-communities" for the dbFD so basically we will be summing and adding all the communities that 
# we used in the NATS into a "big" sample table.

# In our example we had:
# 21 years
# A sequence of 1 11 21 31 41 51 ponds being selected
# 10 iterations 
# This makes a total of 21*6*10= 1260 combinations (AKA a table with 1200 rows)

# Each row of the LakesMergedLakes correspond to one of these combinations  
nrow(LakesMergedLakes)

# We charge the dataset and we eliminate the Lake and Year columns to make it more real to what we will have ;) 
sp_sites <- read.csv2("data/DataAnnaTest/fuzzy_traits.csv",dec = ".") %>% mutate(Year=as.character(Year)) %>% mutate_if(is.numeric,~ifelse(.>0,1,0)) %>% mutate(Year=as.numeric(Year))

# What do we want now? We want to select the ponds that are listed in each row of our table no? So we need to make a "loop"
# where for each row we will select the listed lakes, filter them from the sp_sites and sum their abundance values

# We create the output dataframe that will contain all the combinations of communities
Combis_Communities <- data.frame()
for (combi in 1:nrow(LakesMergedLakes)) {
Lakes_Combin <- LakesMergedLakes[combi,which(is.na(LakesMergedLakes[combi,])==F)] # We extract The row that corresponds to a determined combination 
# Note that we select only the values that are NOT an NA

sp_sites_Combin <- sp_sites %>% filter(Year==as.numeric(LakesMergedLakes[combi,1])) %>% # Filter the sp_sites by the year 
                                filter(Lake%in%Lakes_Combin[4:length(Lakes_Combin)]) %>%  # Filter the lakes by the lakes present in the combination
                                select(-c(Year,lake_year,Lake)) # We "deselect" the three cathegorical metrics

# So, we have generated a dataframe that contains, all the samples (rows) from a year and a combination of ponds. Now we just need to 
# sum by columns and transform the result into 1/0. Finally, bind it to the output dataframe and DONE! :D
Combis_Communities <- bind_rows(Combis_Communities,c(LakesMergedLakes[combi,1:3],ifelse(apply(sp_sites_Combin,2,sum)>0,1,0)))
}

# Each one of the rows corresponds a each one of the used combinations and can be identified with the Year, n_sites and iteration.
# Therefore we can link it with the values of the NATs
Combis_Communities

# In case of need we could also "merge" the three identifiers: 
Combis_Communities %>% mutate(Combi_ID=paste(Year,n_sites,it, sep="_"),.before = Year) %>% select(-c(Year,n_sites,it))

# PISTA!: Sempre que combinis "noms" com has fet amb el lake_year assegurat d'unir-los amb un "_" o un ".", d'aquesta 
# manera, en cas que volguessis separar-los seria fàcil perquè podries identificar per on "partir-los"
# La funció strsplit() va molt bé per fer això però necessita un identificador per on separar ;)



