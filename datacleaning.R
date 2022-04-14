#........................................
# MVIP data cleaning
# Hillary Topazian
# 13 April 2022
#........................................

#........................................####
# Packages & data  ####
#........................................
library(sf)
library(tidyverse)
library(gdistance)
library(abind)
library(rje)
library(malariaAtlas)
library(tableone)
library(lubridate)

# read in cases
cases <- read_csv('./Hospital Surveilance Dataset up to 30 April 2021.csv') 
# read in village coordinates
villages <- read_csv('./Sentinel District Village Cordinates.csv')

# match cases to village coordinates
# village = a combination of district, ta, village and cluster
caselink <- cases %>% left_join(villages, by=c("village_linking_key" = "key"))

nrow(caselink[!is.na(caselink$latitude),]) # 2238 with coordinates
nrow(caselink[!is.na(caselink$village),]) # 7262 with coordinates

caselink <- caselink %>% filter(!is.na(longitude)) %>% 
  st_as_sf(crs=4326, coords=c('longitude','latitude')) %>%
  mutate(hospital_f = factor(hospital, levels=c(1,2,3,4), labels=c('Ntchisi','Mchinji','Balaka','Machinga')))

table(caselink$hospital, caselink$hospital_f)


malawi <- admin0 <- readRDS("./spatial/admin0.rds") %>% filter(Country == 'Malawi') %>% st_transform(4326)
admin0 <- readRDS("./spatial/admin0.rds") %>% filter(Country %in% c('Mozambique', 'Tanzania, United Republic of', 'Zambia', 'Zimbabwe')) %>% st_transform(4326)
load('./spatial/SEAfrica.rdata')

st_crs(malawi) # wgs84
st_crs(admin0)# wgs84

# clusters
balaka <- st_read('./SENTINEL DISTRICTS SHAPEFILES/Balaka/BalakaClusters.shp') %>% st_transform(4326) %>% mutate(ID="Balaka")
machinga <- st_read('./SENTINEL DISTRICTS SHAPEFILES/Machinga/MachingaClustersUpdated_22July2021.shp') %>% st_transform(4326) %>% mutate(ID="Machinga")
ntchisi <- st_read('./SENTINEL DISTRICTS SHAPEFILES/Ntchisi/NtchisiClustersUpdated.shp') %>% st_transform(4326) %>% mutate(ID="Ntchisi")
mchinji <- st_read('./SENTINEL DISTRICTS SHAPEFILES/Mchinji/MchinjiClustersUpdated_22July2021.shp') %>% st_transform(4326) %>% mutate(ID="Mchinji")


# district hospitals
hospital <- c('Balaka', 'Machinga', 'Ntchisi', 'Mchinji')
latitude <- c(-14.98473947, -15.06293002,	-13.3637735, -13.8029785)
longitude <- c(34.94946409, 35.22546935, 33.9096052, 32.8873466)

hospitals <- as_tibble(cbind(hospital, latitude, longitude)) %>% 
  mutate(latitude = as.numeric(latitude), 
         longitude = as.numeric(longitude),
         ID ='district hospital') %>%
  st_as_sf(crs=4326, coords=c('longitude','latitude'))

# plot
ggplot() + 
  geom_sf(data=admin0, fill="cornsilk2", color="cornsilk3") +
  geom_sf(data=malawi, fill="cornsilk") + 
  geom_sf(data=lakes, fill="deepskyblue", color=NA, alpha=0.2) + 
  geom_sf(data=malawi, fill=NA, color="tan4", size=0.75) + 
  geom_sf(data=balaka, aes(fill=ID), color=NA) + 
  geom_sf(data=machinga, aes(fill=ID), color=NA) + 
  geom_sf(data=ntchisi, aes(fill=ID), color=NA) + 
  geom_sf(data=mchinji, aes(fill=ID), color=NA) + 
  geom_sf(data=hospitals, aes(color=ID)) + 
  labs(fill="", color='') + 
  theme_bw(base_size=14) + 
  scale_color_manual(values = rep('black',4)) + 
  scale_x_continuous(limits=c(32.5,36)) + 
  scale_y_continuous(limits=c(-17.2,-9.4)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank()) 


ggsave('./plots/clustermap.pdf', height=5, width=4)


# plot cases to see if they fall near district borders -------------------------
# plot
ggplot() + 
  geom_sf(data=admin0, fill="cornsilk2", color="cornsilk3") +
  geom_sf(data=malawi, fill="cornsilk") + 
  geom_sf(data=lakes, fill="deepskyblue", color=NA, alpha=0.2) + 
  geom_sf(data=malawi, fill=NA, color="tan4", size=0.75) + 
  geom_sf(data=caselink, aes(color=hospital_f), size=.1) +   
  geom_sf(data=hospitals, color='black', size=1) + 
  labs(fill="", color='') + 
  theme_bw(base_size=14) + 
  scale_x_continuous(limits=c(32.5,36)) + 
  scale_y_continuous(limits=c(-15.5,-12.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=1)))

ggsave('./plots/casemap.pdf', height=5, width=6)



# calculate distance -----------------------------------------------------------

calcdist <- function(DISTRICT, N) {
  
subset <- caselink[caselink$hospital_f==DISTRICT,]
dist <- st_distance(subset, hospitals[N,])
cbind(subset, dist)

}

dist_hosp <- map2_dfr(c('Balaka', 'Machinga', 'Ntchisi', 'Mchinji'), seq(1,4,1), calcdist)
summary(dist_hosp$dist)

ggplot(dist_hosp, aes(x=as.character(hospital_f), y=as.numeric(dist)/1000, color=hospital_f, fill=hospital_f)) +
  geom_boxplot(alpha=0.3, show.legend=F) +
  scale_y_continuous(limits=c(0, max(dist_hosp$dist/1000, na.rm=T))) +
  labs(x='District',
       y='Distance (km)') +
  theme_classic() 


ggsave('./plots/distance.pdf', height=4, width=4)



#........................................####
# MAP distance calc  ####
#........................................

# https://medium.com/@abertozz/mapping-travel-times-with-malariaatlas-and-friction-surfaces-f4960f584f08
# download shapefile
analysis.shp <- malariaAtlas::getShp(ISO = "MWI", admin_level = "admin0")
plot(analysis.shp, main="Shape for Clipping")

# download friction surface
friction <- malariaAtlas::getRaster(
  surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015",
  shp = analysis.shp)
malariaAtlas::autoplot_MAPraster(friction)

# tell function how to calculate distance and correct for globe 3D
T <- gdistance::transition(friction, function(x) 1/mean(x), 8) 
T.GC <- gdistance::geoCorrection(T)             

# read in hospital locations
point.locations <- read.csv(file = './hospitals.csv')


# create function to make surface maps
createsurface <- function(district){
  
  # subset to district of interest
  subset <- point.locations %>% filter(name==district)
  
  # keep only point coordinates within file bounds
  coordinates(subset) <- ~ X_COORD + Y_COORD
  proj4string(subset) <- proj4string(subset)
  
  points <- as.matrix(subset@coords)
  
  access.raster <- gdistance::accCost(T.GC, points)
  
  # write the resulting raster
  writeRaster(access.raster, filename=paste0('./spatial/travel_time_',district,'.tif'), overwrite=TRUE)
  
  # extract travel times around cluster points, with a 1m buffer
  travel <- raster::extract(access.raster, caselink, buffer=1, fun=mean)
  
  # plot
  p <- malariaAtlas::autoplot_MAPraster(access.raster, 
                                        shp_df=analysis.shp, printed=F)
  
  full_plot <- p[[1]] + geom_point(data=data.frame(subset@coords), 
                                   aes(x=X_COORD, y=Y_COORD)) +
    scale_fill_gradientn(colors = rev(rje::cubeHelix(gamma=1.0, 
                                                     start=1.5, 
                                                     r=-1.0, 
                                                     hue=1.5, 
                                                     n=16)), 
                         name="Minutes \n of Travel") + 
    ggtitle("Travel Time to Most Accessible Peak") +
    theme(axis.text=element_blank(),
          panel.border=element_rect(fill=NA, color="white"))
  
  print(full_plot)
  print(travel)

}

# create travel time rasters
tt_balaka <- createsurface('Balaka')
tt_machinga <- createsurface('Machinga')
tt_ntchisi <- createsurface('Ntchisi')
tt_mchinji <- createsurface('Mchinji')

# bind results to case dataframe
test <- cbind(caselink, tt_balaka,tt_machinga,tt_ntchisi,tt_mchinji) %>%
  mutate(traveltime = case_when(hospital_f=='Balaka' ~ tt_balaka,
                                hospital_f=='Machinga' ~ tt_machinga,
                                hospital_f=='Ntchisi' ~ tt_ntchisi,
                                hospital_f=='Mchinji' ~ tt_mchinji))

# plot travel time
ggplot(test, aes(x=as.character(hospital_f), y=as.numeric(traveltime), color=hospital_f, fill=hospital_f)) +
  geom_boxplot(alpha=0.3, show.legend=F) +
  scale_y_continuous(limits=c(0, max(100, na.rm=T))) +
  labs(x='District',
       y='Travel time (minutes)') +
  theme_classic() 

ggsave('./plots/travel_time.pdf', height=4, width=4)


#........................................####
# Tables & Figures  ####
#........................................
tabledat <- caselink %>% dplyr::select(
  village_linking_key, hospital, hospital_f, ta, village, cluster, admission_date, sex, 
  agem, firstadmission, admitprevno, child_have_net, 
  sleep_under_net, sprayed_walls, fevernlast7days, preillness, hiv, 
  heart_disease, sickle_cell, malnutrition, breathe_difficulty, alter_consciuos, 
  convulsions, adm_weight, adm_muac, rts1, rts2, 
  rts3, rts4, namediag, namediag1, namediag2, namediag3, outcome, outcome_date)

tabledat <- tabledat %>%
  mutate(sex = case_when(sex==1~'male',
                         sex==0~'female'),
         agem = case_when(agem==99~NA_real_,
                          TRUE~agem),
         preillness = case_when(preillness==2~NA_real_,
                                TRUE~preillness),
         outcome = case_when(outcome==1~'discharged/alive',
                             outcome==2~'dead',
                             outcome==3~'referred',
                             outcome==4~'absconded'),
         namediag = case_when(namediag==1~'Asthma',
                              namediag==2~'Gastroentritis',
                              namediag==3~'Malaria',
                              namediag==4~'Pneumonia',
                              namediag==5~'Severe Malaria',
                              namediag==6~'Septicemia',
                              namediag==7~'Sepsis',
                              namediag==8~'Severe Pneumonia',
                              namediag==9~'Severe anaemia secondary to malaria',
                              namediag==10~'Sickle cell crisis',
                              namediag==11~'Cerebral malaria',
                              namediag==12~'Bacterial Meningitis',
                              namediag==13~'Viral Meningitis',
                              namediag==14~'Cryptococcal Meningitis',
                              namediag==15~'TB Meningitis',
                              namediag==16~'Meningitis Unspecified',
                              namediag==95~'Other'),
         hospstay = interval(lubridate::dmy(admission_date), lubridate::dmy_hm(outcome_date)) %/% days(1)
         
         )

vars <- c('sex', 'preillness', 'outcome', 'namediag')
tab1a <- tableone::CreateCatTable(vars = vars, strata = c('hospital_f'), 
                         data = as.data.frame(tabledat), includeNA = F, test = TRUE, testExact = 'fischer.test')
print(tab1a, showAllLevels = TRUE, formatOptions = list(big.mark = ","))

vars <- c('agem', 'hospstay')
tab1b <- tableone::CreateContTable(vars = vars, strata = c('hospital_f'), 
                                 data = as.data.frame(tabledat), test = TRUE)
print(tab1b, showAllLevels = TRUE, formatOptions = list(big.mark = ","))
