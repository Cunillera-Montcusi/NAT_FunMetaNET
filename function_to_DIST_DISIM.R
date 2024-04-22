#lake_geo <- readxl::read_excel("data/Lakes/lake_coord_surface.xlsx")
#lake_geo <- lake_geo[,1:3] #seleccionem les 3 columnes (llac, x, y)
#lake_geo_m <- geosphere::distm(lake_geo[,2:3]) #fem la matriu de dist en metres

#orig_lake <- sample(lake_geo$Lake, 1) #seleccionem un llac aleatori
#pos_lake <- which(lake_geo$Lake==orig_lake) #demanem la posició del llac a la llista
#A<- as.matrix(lake_geo_m)[,pos_lake]# Aillem la columna que correspon al llac Original
#B <-A[-pos_lake] #Eliminem el zero que està al mateix lloc que el llac original
#Dist_Neigh_lake <- min(B) # Trobem la distànica minima entre Original i altres
#c=(which(A==Dist_Neigh_lake))#trobem la posicio del llac mes proper
#lake_geo$Lake[c]#trobem el nom del llac mes proper amb la seva posicio al df inicial
#d <- subset(lake_geo, Lake%in%c(orig_lake,lake_geo$Lake[c]))#filtrem els dos llacs que estem comparant

#lake_geo_new <- lake_geo[-c(pos_lake,c),]#treiem de la df inicial el llac original i el mes proper

#new_lake <- data.frame(Lake="New_Lake",
#                      summarise(d,Lon=mean(Lon),
#                                 Lat=mean(Lat)),
#                      year=NA)# fem la mitjana dels dos llacs i creem el new lake
#lake_geo_new <- lake_geo %>% filter(!Lake%in%c(orig_lake,lake_geo$Lake[c])) %>% 
#  bind_rows(new_lake)#treiem els llacs dels quals hem fet la mean i afegim el new lake

fun_to_DIST_DISIM <- function(orig_lake ){
  
  lake_geo <- readxl::read_excel("data/Lakes/lake_coord_surface.xlsx") 
  lake_geo <- lake_geo[,1:3]
  #orig_lake <- sample(lake_geo$Lake, 1)#agafem un llac a l'atzar per començar
  #loop: de mes semblant a mes diferent
  Llista_llacs <- c() #obrim un vestor on guradar els llacs per ordre
  Llista_llacs[1] <- orig_lake #posem el 1r llac al vector
  lake_geo$Lake[which(lake_geo$Lake==orig_lake)] <- "Orig_Lake" #retorna el nom del llac triat
  
  for (coor in 2:nrow(lake_geo)) {
    
    orig_lake <- lake_geo %>% filter(Lake=="Orig_Lake") %>% pull(Lake)
    
    # Primer calculem distancia
    xy <- lake_geo[,2:3]  
    dist_lakes <- geosphere::distm(xy)
    
    # Segon
    pos_lake <- which(lake_geo$Lake==orig_lake) # Posició del Orig lake
    A <- as.matrix(dist_lakes)[,pos_lake]# Aillem la columna que correspon al llac Original
    B <-A[-pos_lake] #Eliminem el zero que està al mateix lloc que el llac original
    Dist_Neigh_lake <- max(B) # Trobem la distànica minima entre Original i altres
    C <- (which(A==Dist_Neigh_lake)) # Localitzem la posició del mínim
    Closer_Lake <- lake_geo$Lake[C] # Nom del llac més proper a l'original
    
    # Tercer 
    d <- subset(lake_geo, Lake%in%c(orig_lake,Closer_Lake))#filtrem els dos llacs més propers
    new_lake <- data.frame(Lake="Orig_Lake",
                           summarise(d,Lon=mean(Lon),
                                     Lat=mean(Lat)),
                           year=NA) #fem la mitjana dels dos llacs més propers
    
    # Quart 
    lake_geo <- lake_geo %>% filter(!Lake%in%c(orig_lake,Closer_Lake)) %>% 
      bind_rows(new_lake) #treiem de la llista els dos llacs dels quals hem fer mean, i afegim "new lake" (mean dels dos llacs)
    
    Llista_llacs[coor] <- Closer_Lake
  }
  Llista_llacs  
}# End function

#Rivers####
fun_to_DIST_DISIM_riv <- function(orig_lake ){
  
  lake_geo <- readxl::read_excel("data/Rivers/swed_rivers/coord_swed_riv.xlsx")
  #orig_lake <- sample(lake_geo$Lake, 1)#agafem un llac a l'atzar per començar
  #loop: de mes semblant a mes diferent
  Llista_llacs <- c() #obrim un vestor on guradar els llacs per ordre
  Llista_llacs[1] <- orig_lake #posem el 1r llac al vector
  lake_geo$River[which(lake_geo$River==orig_lake)] <- "Orig_Lake" #retorna el nom del llac triat
  
  for (coor in 2:nrow(lake_geo)) {
    
    orig_lake <- lake_geo %>% filter(River=="Orig_Lake") %>% pull(River)
    
    # Primer calculem distancia
    xy <- lake_geo[,2:3]  
    dist_lakes <- geosphere::distm(xy)
    
    # Segon
    pos_lake <- which(lake_geo$River==orig_lake) # Posició del Orig lake
    A <- as.matrix(dist_lakes)[,pos_lake]# Aillem la columna que correspon al llac Original
    B <-A[-pos_lake] #Eliminem el zero que està al mateix lloc que el llac original
    Dist_Neigh_lake <- max(B) # Trobem la distànica minima entre Original i altres
    C <- (which(A==Dist_Neigh_lake)) # Localitzem la posició del mínim
    Closer_Lake <- lake_geo$River[C] # Nom del llac més proper a l'original
    
    # Tercer 
    d <- subset(lake_geo, River%in%c(orig_lake,Closer_Lake))#filtrem els dos llacs més propers
    new_lake <- data.frame(River="Orig_Lake",
                           summarise(d,Longitude_X=mean(Longitude_X),
                                     Latitude_Y=mean(Latitude_Y))) #fem la mitjana dels dos llacs més propers
    
    # Quart 
    lake_geo <- lake_geo %>% filter(!River%in%c(orig_lake,Closer_Lake)) %>% 
      bind_rows(new_lake) #treiem de la llista els dos llacs dels quals hem fer mean, i afegim "new lake" (mean dels dos llacs)
    
    Llista_llacs[coor] <- Closer_Lake
  }
  Llista_llacs  
}# End function
