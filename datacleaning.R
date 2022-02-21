library(excel.link)
library(sf)
library(tidyverse)

# pass <- read.table('./pass.txt') %>% as.character()
#   
# data <- excel.link::xl.read.file('C:/Users/htopazia/Documents/MVIP/Hospital Surveilance Dataset up to 30 April 2021.xlsx', password = pass)

# read in cases
cases <- read_csv('C:/Users/htopazia/Documents/MVIP/Copy of Hospital Surveilance Dataset up to 30 April 2021.csv') 
# read in village coordinates
villages <- read_csv('C:/Users/htopazia/Documents/MVIP/Sentinel District Village Cordinates.csv')

# match cases to village coordinates
# Note that a village for us is defined as a combination of district, ta, village and cluster (as combined in the key column). In your file merge exercise please keep this in mind because some village names repeat in several clusters/districts and can only be distinguished through this key combination.

nrow(cases[!is.na(cases$village),]) # 7262 with village names
nrow(cases[is.na(cases$village),]) 

caselink <- cases %>% left_join(villages, by=c("village" = "Village"))

nrow(caselink[!is.na(caselink$latitude),]) # 5223 with coordinates

caselink <- caselink %>% filter(!is.na(longitude)) %>% 
  st_as_sf(crs=4326, coords=c('longitude','latitude')) %>%
  mutate(hospital_f = factor(hospital, levels=c(1,2,3,4), labels=c('Ntchisi','Mchinji','Balaka','Machinga')))

table(caselink$hospital, caselink$hospital_f)


# recovering by hamming's distance
barcodematedist <- stringdist::stringdistmatrix(
  a=cases$village,
  b=villages$Village,
  method=c("hamming"),
  useNames = c("strings"))
barcodematedist_tbl <- cbind.data.frame(rownames(barcodematedist), barcodematedist)
colnames(barcodematedist_tbl)[1] <- "HIVrecode_barcode"
barcodematedist_tbl_oneoffs <- barcodematedist_tbl %>% 
  as.tibble(.) %>% 
  tidyr::gather(., key="qPCR_barcode", value = "stringdist", 2:ncol(.)) %>% 
  dplyr::filter(stringdist == 1) %>% 
  #dplyr::filter( ! c( HIVrecode_barcode == "C3S3Z" & qPCR_barcode == "C2S3Z" ) ) %>%  # manual fix, see below 
  dplyr::select(-c(stringdist))


malawi <- admin0 <- readRDS("./spatial/admin0.rds") %>% filter(Country == 'Malawi') %>% st_transform(4326)
admin0 <- readRDS("./spatial/admin0.rds") %>% filter(Country %in% c('Mozambique', 'Tanzania, United Republic of', 'Zambia', 'Zimbabwe')) %>% st_transform(4326)
load('./spatial/SEAfrica.rdata')

st_crs(malawi) # wgs84
st_crs(admin0)# wgs84

# clusters
balaka <- st_read('./SENTINEL DISTRICTS SHAPEFILES/Balaka/BalakaClusters.shp') %>% st_transform(4326) %>% mutate(ID="Balaka")
machinga <- st_read('./SENTINEL DISTRICTS SHAPEFILES/Machinga/MachingaClusters.shp') %>% st_transform(4326) %>% mutate(ID="Machinga")
ntchisi <- st_read('./SENTINEL DISTRICTS SHAPEFILES/Ntchisi/NtchisiClusters.shp') %>% st_transform(4326) %>% mutate(ID="Ntchisi")
mchinji <- st_read('./SENTINEL DISTRICTS SHAPEFILES/Mchinji_New/MchinjiClusters.shp') %>% st_transform(4326) %>% mutate(ID="Mchinji")

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



# tESTS_---------------

district <- st_read("./GADM/gadm36_MWI_1.shp",stringsAsFactors=F) # ggplot() + geom_sf(data=district)


ggplot() + 
  geom_sf(data=balaka) + 
  geom_sf(data=machinga) + 
  geom_sf(data=ntchisi) +
  geom_sf(data=mchinji)
