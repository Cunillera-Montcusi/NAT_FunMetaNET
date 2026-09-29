#function for lakes data####
fun_to_ENV_SIM <- function(ref_year,orig_lake ){
  
  lake_pca <- readxl::read_excel("ClustNAT/data/lake_pca_coord.xlsx")
  
  lake_year <- dplyr::filter(lake_pca, year==as.character(ref_year))#filtem any
  #orig_lake <- sample(lake_year$Lake, 1)#agafem un llac a l'atzar per començar
  #loop: de mes semblant a mes diferent
  Llista_llacs <- c() #obrim un vestor on guradar els llacs per ordre
  Llista_llacs[1] <- orig_lake #posem el 1r llac al vector
  lake_year$Lake[which(lake_year$Lake==orig_lake)] <- "Orig_Lake" #retorna el nom del llac triat
  
  for (coor in 2:nrow(lake_year)) {
    
    orig_lake <- lake_year %>% filter(Lake=="Orig_Lake") %>% pull(Lake)
    
    # Primer calculem distancia
    xy <- lake_year[,2:3]  
    dist_lakes <- dist(xy,method = "euclidean")
    
    # Segon
    pos_lake <- which(lake_year$Lake==orig_lake) # Posició del Orig lake
    A <- as.matrix(dist_lakes)[,pos_lake]# Aillem la columna que correspon al llac Original
    B <-A[-pos_lake] #Eliminem el zero que està al mateix lloc que el llac original
    Dist_Neigh_lake <- min(B) # Trobem la distànica minima entre Original i altres
    C <- (which(A==Dist_Neigh_lake)) # Localitzem la posició del mínim
    Closer_Lake <- lake_year$Lake[C] # Nom del llac més proper a l'original
    
    # Tercer 
    d <- subset(lake_year, Lake%in%c(orig_lake,Closer_Lake))#filtrem els dos llacs més propers
    new_lake <- data.frame(Lake="Orig_Lake",
                           summarise(d,Dim.1=mean(Dim.1),
                                     Dim.2=mean(Dim.2)),
                           year=NA) #fem la mitjana dels dos llacs més propers
    
    # Quart 
    lake_year <- lake_year %>% filter(!Lake%in%c(orig_lake,Closer_Lake)) %>% 
      bind_rows(new_lake) #treiem de la llista els dos llacs dels quals hem fer mean, i afegim "new lake" (mean dels dos llacs)
    
    Llista_llacs[coor] <- Closer_Lake
  }
  Llista_llacs  
}# End function

#function for river data####
fun_to_ENV_SIM_riv <- function(ref_year,orig_lake ){
  
  lake_pca <- readxl::read_excel("ClustNAT/data/riv_pca_coord.xlsx")
  
  lake_year <- dplyr::filter(lake_pca, year==as.character(ref_year))#filtem any
  #orig_lake <- sample(lake_year$Lake, 1)#agafem un llac a l'atzar per començar
  #loop: de mes semblant a mes diferent
  Llista_llacs <- c() #obrim un vestor on guradar els llacs per ordre
  Llista_llacs[1] <- orig_lake #posem el 1r llac al vector
  lake_year$site[which(lake_year$site==orig_lake)] <- "Orig_Lake" #retorna el nom del llac triat
  
  for (coor in 2:nrow(lake_year)) {
    
    orig_lake <- lake_year %>% filter(site=="Orig_Lake") %>% pull(site)
    
    # Primer calculem distancia
    xy <- lake_year[,2:3]  
    dist_lakes <- dist(xy,method = "euclidean")
    
    # Segon
    pos_lake <- which(lake_year$site==orig_lake) # Posició del Orig lake
    A <- as.matrix(dist_lakes)[,pos_lake]# Aillem la columna que correspon al llac Original
    B <-A[-pos_lake] #Eliminem el zero que està al mateix lloc que el llac original
    Dist_Neigh_lake <- min(B) # Trobem la distànica minima entre Original i altres
    C <- (which(A==Dist_Neigh_lake)) # Localitzem la posició del mínim
    Closer_Lake <- lake_year$site[C] # Nom del llac més proper a l'original
    
    # Tercer 
    d <- subset(lake_year, site%in%c(orig_lake,Closer_Lake))#filtrem els dos llacs més propers
    new_lake <- data.frame(site="Orig_Lake",
                           summarise(d,Dim.1=mean(Dim.1),
                                     Dim.2=mean(Dim.2)),
                           year=NA) #fem la mitjana dels dos llacs més propers
    
    # Quart 
    lake_year <- lake_year %>% filter(!site%in%c(orig_lake,Closer_Lake)) %>% 
      bind_rows(new_lake) #treiem de la llista els dos llacs dels quals hem fer mean, i afegim "new lake" (mean dels dos llacs)
    
    Llista_llacs[coor] <- Closer_Lake
  }
  Llista_llacs  
}# End function


# lake_pca <- readxl::read_excel("data/Lakes/lake_pca_coord.xlsx")
# 
# lake_year <- dplyr::filter(lake_pca, year=="1995")
# 
# xy <- lake_year[,2:3]
# orig_lake <- sample(lake_year$Lake, 1)
# 
# dist_lakes <- dist(xy,method = "euclidean")
# 
# pos_lake <- which(lake_year$Lake==orig_lake)
# A<- as.matrix(dist_lakes)[,pos_lake]# Aillem la coumna qye cirrespon al llac Original
# B <-A[-pos_lake] #Eliminem el zero que està al mateix lloc que el llac original
# Dist_Neigh_lake <- min(B) # Trobem la distànica minima entre Original i altres
# c=(which(A==Dist_Neigh_lake))
# lake_year$Lake[c]
# d <- subset(lake_year, Lake%in%c(orig_lake,lake_year$Lake[c]))
# 
# 
# lake_year_new <- lake_year[-c(pos_lake,c),]
# 
# new_lake <- data.frame(Lake="New_Lake",
#                        summarise(d,Dim.1=mean(Dim.1),
#                                  Dim.2=mean(Dim.2)),
#                        year=NA)
# 
# lake_year_new <- lake_year %>% filter(!Lake%in%c(orig_lake,lake_year$Lake[c])) %>% 
#   bind_rows(new_lake)
# 











