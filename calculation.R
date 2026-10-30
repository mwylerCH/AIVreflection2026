# Results of AIV clustering

library(data.table)
library(tidyverse)
library(Rtsne)
library(ggbiplot)
library(umap) 


## Metadata Tab --------------------

### genotype name
metaTab <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/temp_Metadata.txt',
                 data.table = F, header = F, fill = T, sep = "\t")
colnames(metaTab) <- c('ID', 'Isolate')


# subset for only analyzed (IUPA+8 segments only, single genome)
fatti <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/analizedSamples.txt',
               data.table = F, header = F, fill = T, sep = "\t")
metaTab <- metaTab %>% 
  filter(ID %in% fatti$V1)

### Year extraction
metaTab$Year <- metaTab$Isolate

# extract year
metaTab$Year <- gsub("Mixed$", '', metaTab$Year)
metaTab$Year <- gsub("(.*)H[[:digit:]]N[[:digit:]]$", "\\1", metaTab$Year) # if HxNx at the end
metaTab$Year <- gsub("(.*)_+$", "\\1", metaTab$Year) # if HxNx at the end
metaTab$Year <- gsub("(.*)__H[[:digit:]]N[[:digit:]]$", "\\1", metaTab$Year) # if HxNx at the end
metaTab$Year <- gsub(".*/([[:digit:]]{4}$)", "\\1", metaTab$Year)
metaTab$Year <- gsub(".*/([[:digit:]]{4})-[[:digit:]]{2}-[[:digit:]]{2}$", "\\1", metaTab$Year) # if detailed date
metaTab$Year <- gsub(".*-([[:digit:]]{4}$)", "\\1", metaTab$Year)

# remove remaining
metaTab$Year[metaTab$Year < 1900] <- NA
metaTab$Year[grepl('[[:alpha:]]', metaTab$Year)] <- NA


### Geo vector
library(r2country)

# prep package data
continentReference <- continentOf %>% unlist %>% as.data.frame()
continentReference$country <- rownames(continentReference)
colnames(continentReference)[1] <- 'continent'
continentReference$country <- gsub(' ','_', continentReference$country)

# more then 200 samples from russian regions, but also many cities
russianProvince <- grep("russia,",names(timeIn), value = T)
russianProvince <- gsub('russia, ', '', russianProvince)
russianProvince <- c(russianProvince, 'rostov', 'ossetia', 'krasnodar', 'leningrad', 'saratov', 'kostroma',
                     'voronezh', 'tyumen', 'tatarstan', 'kursk', 'kurgan', 'penza', 'tyva')

# same for china
chineseProvince <- grep("china,",names(timeIn), value = T)
chineseProvince <- gsub('china, ', '', chineseProvince)
chineseProvince <- gsub(', ', '\\|', chineseProvince)


# us states
USstates <- grep("usa,",names(timeIn), value = T)
USstates <- gsub('usa, ', '', USstates)
USstates <- gsub(', .*', '', USstates)
USstates <- gsub(' ', '_', USstates)

# get ID of samples (lower capital)
metaTab$Country <- tolower(metaTab$Isolate)

# UK vs single states
metaTab$Country  <- gsub('england', 'united_kingdom', metaTab$Country )
metaTab$Country  <- gsub('scotland', 'united_kingdom', metaTab$Country )
metaTab$Country  <- gsub('wales', 'united_kingdom', metaTab$Country )
metaTab$Country  <- gsub('nl-', 'netherlands', metaTab$Country )
metaTab$Country  <- ifelse(grepl(paste0(russianProvince, collapse = '|'), metaTab$Country ), 'russia', metaTab$Country )
metaTab$Country  <- ifelse(grepl(paste0(chineseProvince, collapse = '|'), metaTab$Country ), 'china', metaTab$Country )
metaTab$Country  <- ifelse(grepl(paste0(USstates, collapse = '|'), metaTab$Country ), 'united_states', metaTab$Country )

# convert to continent
continentReference <- continentReference[continentReference$country != "turkey",] # exclude turkey too dangerous
for (ROW in 1:nrow(continentReference)){
  continente <- continentReference$continent[ROW]
  paese <- continentReference$country[ROW]
  # pick IDs of interest
  metaTab$Country[grepl(paese, metaTab$Country)] <- continente
  
}
table(metaTab$Country)

# manually correct 
metaTab$Country  <- ifelse(grepl('korea', metaTab$Country), 'Asia', metaTab$Country)
metaTab$Country  <- ifelse(grepl('chiba', metaTab$Country), 'Asia', metaTab$Country)
metaTab$Country  <- ifelse(grepl('czech_reublic', metaTab$Country), 'Europe', metaTab$Country)
metaTab$Country  <- ifelse(grepl('victoria|tasmania|northern_territory|queensland"', metaTab$Country), 'Oceania', metaTab$Country)


# make some vectors of remaining
asianRegions <- c(
  "aichi", "akita", "ang_thong", "angthong", "anhui", "aomori", "anyang", "assam", 
  "ba_ria-vung_tau", "ba_ria_vung_tau", "ba_ria_vungtau", "ba_rja_vung_tau", 
  "bac_giang", "bac_gieng", "bac_lieu", "bac_ninh", "bali", "bandar_lampung", 
  "bangkok", "bangli_bali", "banten", "bantul", "ben_tre", "binh_dinh", "binh_duong", 
  "binh_thuan", "boyolali", "ca_mau", "camau", "can_tho", "cao_bang", "changhua", 
  "chany", "chany_lake", "chiayi", "chiping", "da_nang", "dairi", "dak_lak", 
  "deli_serdang", "dong_nai", "dong_thap", "dongguan", "dongnai", "dovolnoe", 
  "dubai", "east_java", "fukuoka", "fukushima", "ganzhou", "garut", "gifu", 
  "gui_yang", "gunung_kidal", "ha_nam", "ha_tay", "ha_tinh", "hai_duong", 
  "hai_phong_", "hainan", "hamanaka", "hanoi", "hau_giang", "hiroshima", "hk", 
  "hn", "hokkaido", "huaian", "huadong", "hulu_sungai_utara", "hunnan", "hyogo", 
  "ibaraki", "indramayu", "iwate", "jalgaon", "jembrana", "jiang_xi", "kagawa", 
  "kagoshima", "kalasin", "kampot", "kanagawa", "kedaton", "khabarovsk", "khanh_hoa", 
  "khunt_lake", "khuvsgul", "kien_giang", "kochi", "kon_tum", "krasnoozerskoe", 
  "kulon_progo", "kumamoto", "kyaing_tong", "kyoto", "lahore", "laiwu", "lang_son", 
  "lao", "long_an", "longan", "magetan", "malang", "mangarai-ntt", "mie", "miyazaki", 
  "mong_cai", "nagano", "nakhon_pathom", "nakhonsawan", "nam_dinh", "nandurbar", 
  "nara", "ngawi", "nghe_an", "nha_trang", "niigata", "ninh_binh", "ninh_thuan", 
  "ningxia", "north_sumatra", "oita", "okayama", "pangkalpinang", "parepare", 
  "pavlodar", "pekalongan", "ph", "phu_tho", "pingtung", "pulau_rampang", "purwakarta", 
  "purworejo", "qinghai", "quang", "quang_nam", "quang_ngai", "quang_ninh", 
  "quang_tri", "quanh_ninh", "saga", "sagaing", "saitama", "sakha", "salatiga", 
  "sanjiang", "sanmenxia", "sartlan", "shiga", "shimane", "shizuoka", 
  "sidenreng_rappang", "sikkim", "simalanggang", "soc_trang", "son_la", 
  "south_kalimantan", "sragen", "st", "suphanburi", "suzdalka", "taichung", 
  "taishun", "taizhou", "tainan", "takeo", "taoyuan", "tarutung", "tasikmalaya", 
  "tay_ninh", "tebing_tinggi", "thai_binh", "thai_ninh", "thanatpin", "thanh_hoa", 
  "tien_giang", "tochigi", "tokushima", "tokyo", "tonghai", "tottori", "toyama", 
  "tripura", "tsukuba", "tumen", "tuva", "ugii_lake", "viet_nam", "vinh_long", 
  "wajo", "wakayama", "way_kanan", "wenzhou", "west_bengal", "west_java", 
  "westjava", "wuxi", "xianghai", "xuzhou", "yamaguchi", "yangon", "yangzhou", 
  "yilan", "yn", "yogjakarta", "yuhang", "zhenjiang", "gaza", "taipei", "kaohsiung"
)
northAmericaRegions <- c(
  "ab", "aguascalientes", "aguascallientes", "alb", "alberta", "bc", 
  "british_columbia", "ca", "chiapas", "coahuila", "colima", "de", "durango", 
  "guanajuato", "hidalgo", "jalisco", "ma", "massachussetts", "mb", "md", "me", 
  "mi", "michoacan", "mn", "morelos", "new_brunswick", "newyork", "nj", "ns", 
  "nv", "on", "ontario", "pa", "pei", "pensylvania", "puebla", "qc", "queretaro", 
  "queretaroa", "queretero", "san_luis_potosi", "saskatchewan", "sk", "sonora", 
  "tabasco", "tx", "va", "veracruz", "wi", "yucatan"
)
europeRegions <- c("adygea", "altai", "bih", "castillalamancha", "crimea", "czech_reublic", 
                   "djankoy", "domodedovo", "jersey", "kalmykia", "kamchatka", "lisbon", "nl", 
                   "poitou_charentes", "potsdam", "primorje", "pt", "reshoty", "sergiyev_posad", 
                   "uvs-nuur_lake", "berlin")

africaRegions <- c(
  "al-sharqia", "albehra", "beni_sueif", "cairo", "cote_d'ivoire", "dakahlia", 
  "el_fayoum", "giza", "ivory_coast", "kafr_elsheikh", "menia", "menofia", 
  "qalubia", "sharkia", "sharkeya"
)
southAmericaRegions <- c("antofagasta", "araucania", "arica_y_parinacota", "atacama", "bolivar", 
                         "choco", "coquimbo", "cordoba", "magdalena", "maule", "nuble", "ohiggins", 
                         "tarapaca", "valparaiso")

metaTab$Country  <- ifelse(grepl(paste0("/", northAmericaRegions, "/", collapse = '|'), metaTab$Country ), 'North America', metaTab$Country )
metaTab$Country  <- ifelse(grepl(paste0(africaRegions, collapse = '|'), metaTab$Country ), 'Africa', metaTab$Country )
metaTab$Country  <- ifelse(grepl(paste0(asianRegions, collapse = '|'), metaTab$Country ), 'Asia', metaTab$Country )
metaTab$Country  <- ifelse(grepl(paste0(europeRegions, collapse = '|'), metaTab$Country ), 'Europe', metaTab$Country )
metaTab$Country  <- ifelse(grepl(paste0(southAmericaRegions, collapse = '|'), metaTab$Country ), 'South America', metaTab$Country )

# remove bad samples with no geo position
metaTab[!metaTab$Country %in% continentReference$continent, 'Country'] <- NA


# Eu classification name
EUTab <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/EU_classification.txt',
               data.table = F, header = F, fill = T, sep = "\t")
EUTab$V1 <- gsub('EPIISL', 'EPI_ISL_', EUTab$V1)
colnames(EUTab) <- c('ID', 'EU_class')
EUTab$EU_class <- gsub('\\[unassigned\\]', NA, EUTab$EU_class)


# US classification
USTab <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/US_classification.txt',
               data.table = F, header = T, fill = T, sep = "\t")
colnames(USTab) <- c('ID', 'US_class')
USTab$ID <- gsub('EPIISL', 'EPI_ISL_', USTab$ID)
USTab$ID <- gsub('.fa', '', USTab$ID)
USTab$US_class <- ifelse(grepl('Not assigned', USTab$US_class), NA, USTab$US_class)


### Merge tables
metaTab <- left_join(metaTab,EUTab) %>% 
  left_join(., USTab)


colnames(PCAbackup)[duplicated(colnames(PCAbackup))]
# AVOR Distances -------------------------------------

# prep empty tab
PCAtable <- metaTab[,1:2]
PCAtable$Isolate <- NULL

# representative sequences
nomi <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/ReferenceNames.txt',
              header = F, sep = '\t')
nomi[is.na(nomi$V2),2] <- 'NA'

# read in segment distance
for (NR in 1:8){
  nameSeg <- c('PB2', 'PB1', 'PA', 'HA','NP', 'NA','MP','NS')
  nome <- paste0('//wsl.localhost/Ubuntu-22.04/home/mwyler/temp_AIVpca/Segment_', nameSeg[NR], '_2.DistMatrix')
  pezzo <- fread(nome,
                 data.table = F)
  
  # filter using representative sequences
  keepNames <- c('V1', as.vector(unlist(nomi[nomi$V2 == nameSeg[NR], 1])))
  pezzo <- pezzo[, colnames(pezzo) %in% keepNames]
  
  # clean names rows
  pezzo$V1 <- gsub('_.*$', '', pezzo$V1)
  pezzo$V1 <- gsub('EPIISL', 'EPI_ISL_', pezzo$V1)
  
  
  # order based on classification file
  ordine <- pezzo[match(metaTab$ID, pezzo$V1),]
  
  # remove duplicated columns
  ordine <- ordine[, !duplicated(colnames(ordine))]
  
  # remove ID
  ordine$V1 <- NULL
  
  PCAtable <- cbind(PCAtable, ordine)
  print(dim(PCAtable))
}

PCAbackup <- PCAtable
# colnames(PCAbackup)[duplicated(colnames(PCAbackup))]

# some are NAs (removed because of trimming)
removeFromMeta <- PCAtable %>% filter(if_any(everything(), is.na)) %>% pull(ID)
PCAtable <- PCAtable %>% na.omit()
metaTab <- metaTab %>% 
  filter(!ID %in% removeFromMeta)

# duplicate IDs table
PCAtable %>%
  group_by(across(2:ncol(PCAtable))) %>%
  filter(n() > 1) %>%
  ungroup()

# make table with duplicated IDs
noId <- PCAtable[, 2:ncol(PCAtable)]
key <- do.call(paste, c(noId, sep = "\r"))
matches <- duplicated(key)
duplicatedIDs <- data.frame(duplicated_row = which(matches),
                            same_as_row = match(key[matches], key))
for (NR in duplicatedIDs$same_as_row){
  duplicatedIDs[duplicatedIDs$same_as_row == NR, 2] <- PCAtable$ID[NR]
}
for (NR in duplicatedIDs$duplicated_row){
  duplicatedIDs[duplicatedIDs$duplicated_row == NR, 1] <- PCAtable$ID[NR]
}

# remove duplicated rows
duplicatedRows <- duplicated(noId)
PCAtable <- PCAtable[!duplicatedRows,]

# find identical columns
PCAtable <- PCAtable[!duplicated(as.list(PCAtable))]


## RTsne -----------------------------------------

# prepare matrix
matrice <- as.matrix(PCAtable[, 2:ncol(PCAtable)])
X <- normalize_input(matrice)

#set.seed(12345)
#tsneRes <- Rtsne(X,
# perplexity should been between 5 and 50
#                 perplexity = 50,
#                 num_threads = 20)

set.seed(1234)
tsneRes <- Rtsne(X,
                 # perplexity should been between 5 and 50
                 perplexity = 50,
                 early_exaggeration = 20,
                 pca = F,
                 num_threads = 20)

# alternative UMAP
#umapRes = umap(X, n_components = 2, random_state = 15) 
#layout <- umapRes[["layout"]] 
#layout <- data.frame(layout) 
#final <- cbind(layout, PCAtable$V3) 
#colnames(final) <- c('x','y','col')

# AVOR plot Tab
tsne_plot <- data.frame(x = tsneRes$Y[,1], y = tsneRes$Y[,2], ID = PCAtable$ID)
# add metadata
tsne_plot <- tsne_plot %>% 
  left_join(., metaTab)

# plot -----------------------------------------

### Colors ----------------------

# select random colors for an ad hoc palette
colori <- colors()
colori <- colori[!grepl('grey|gray', colori)]

# EU colors
EUcolors <- sample(colori, length(unique(tsne_plot$EU_class)), replace = F)

# US colors
UScolors <- sample(colori, length(unique(tsne_plot$US_class)), replace = F)

# inHouse
#only after DB scan

### Continents ------------------
ContinentPlot <- tsne_plot %>% 
  ggplot() + 
  # Note: Make sure to color by Continent/Region if you only want a short legend here
  geom_point(aes(x = x, y = y, color = Country), alpha = 0.7) + 
  scale_color_discrete(na.value = "grey90") +
  theme_bw() +
  labs(title = 'Geographic Distribution', x = "t-SNE 1", y = "t-SNE 2") +
  theme(
    # 1. Place the legend INSIDE the plot area (X=15% from left, Y=85% from bottom)
    legend.position = c(0.15, 0.85), 
    # 2. Make the legend box background transparent so it doesn't block data points
    legend.background = element_rect(fill = alpha("white", 0.6), color = "grey80"),
    # 3. Shrink text and spacing to make it very compact
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(3, 5, 3, 5),
    
    plot.margin = margin(15, 15, 15, 15)
  )
ContinentPlot
ggsave('MWy/forPaper_AIV/tsne_AIV_country.pdf', width = 10, height = 10)


Continentlegend <- tsne_plot %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=Country))+
  scale_color_discrete(na.value = "grey90")+
  theme_bw()+
  labs(title = 'Continents')+
  theme(legend.position = 'bottom', legend.title = element_blank()) 
cowplot::plot_grid(cowplot::get_legend(Continentlegend))
ggsave('MWy/forPaper_AIV/Countrylegend.pdf', width = 10, height = 10)


### European ------------------


EUplot <- tsne_plot %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=EU_class))+
  scale_color_manual(values = EUcolors, na.value = "grey90")+
  theme_bw()+
  labs(title = 'EURL Classification', x = "t-SNE 1", y = "t-SNE 2")+
  theme(legend.position = '',
        plot.margin = margin(15, 15, 15, 15))
EUplot
ggsave('MWy/forPaper_AIV/tsne_AIV_European.pdf', width = 10, height = 10)

EUlegend <- tsne_plot %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=EU_class))+
  scale_color_manual(values = EUcolors, na.value = "grey90")+
  theme_bw()+
  labs(title = 'European')+
  theme(legend.position = 'bottom', legend.title = element_blank()) 
cowplot::plot_grid(cowplot::get_legend(EUlegend))
ggsave('MWy/forPaper_AIV/EUlegend.pdf', width = 10, height = 10)


### US ------------------

USplot <- tsne_plot %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=US_class))+
  scale_color_manual(values = UScolors, na.value = "grey90")+
  theme_bw()+
  labs(title = 'USDA Classification', x = "t-SNE 1", y = "t-SNE 2")+
  theme(legend.position = '',
        plot.margin = margin(15, 15, 15, 15))
USplot
ggsave('MWy/forPaper_AIV/tsne_AIV_American.pdf', width = 10, height = 10)


USlegend <- tsne_plot %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=US_class))+
  scale_color_manual(values = UScolors, na.value = "grey90")+
  theme_bw()+
  labs(title = 'us')+
  theme(legend.position = 'bottom', legend.title = element_blank()) 
cowplot::plot_grid(cowplot::get_legend(USlegend))
ggsave('MWy/forPaper_AIV/USlegend.pdf', width = 10, height = 10)


### DBscan ------------------


library(dbscan)

# get ellbow (for minPTS, k is given by dimensionality)
tsne2Dim <- tsne_plot[, 1:2]
kNNdistplot(tsne2Dim, k = 10) # ln(rows)
abline(h=1, col = "red", lty=2)

# cluster
dbscan_res <- dbscan(tsne2Dim, eps = 1.5, minPts = 10)
dbscan_res
dbscan_res$cluster

# merge table
tsne_plot$cluster <- as.factor(dbscan_res$cluster)
inHousecolors <- sample(c(colori,colori,colori), length(unique(tsne_plot$cluster))-1, replace = F)

DBplot <- tsne_plot %>% 
  mutate(cluster = ifelse(cluster == "0", NA, cluster)) %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  ggplot() + 
  geom_point(aes(x=x, y=y, color=cluster)) +
  scale_color_manual(values = inHousecolors, na.value = "grey90")+
  theme_bw()+
  labs(title = 'Statistical Genotyping', x = "t-SNE 1", y = "t-SNE 2")+
  theme(legend.position = '',
        plot.margin = margin(15, 15, 15, 15))
DBplot
ggsave('MWy/forPaper_AIV/tsne_AIV_DBclust.pdf', width = 10, height = 10)


### Panel ------------------

cowplot::plot_grid(ContinentPlot, EUplot,
                   USplot, DBplot,
                   labels = c('A', 'B', 'C', 'D'))
ggsave('MWy/forPaper_AIV/panel.pdf', width = 15, height = 15)
ggsave('MWy/forPaper_AIV/panel.svg', width = 15, height = 15)
ggsave('MWy/forPaper_AIV/panel.tiff', width = 15, height = 15)
ggsave('MWy/forPaper_AIV/panel.png', width = 15, height = 15)


# cluster centroids ------------------

# data for the metadata table


# calculate centroids for all clusters
centroids <- tsne_plot %>%
  filter(cluster != 0) %>% 
  group_by(cluster) %>%
  summarize(centroid_x = mean(x, na.rm = TRUE),
            centroid_y = mean(y, na.rm = TRUE))

# get closest entry (euclidean)
closest_samples <- tsne_plot %>%
  filter(cluster != 0) %>%
  left_join(centroids, by = "cluster") %>%
  mutate(distance = sqrt((x - centroid_x)^2 + (y - centroid_y)^2)) %>%
  group_by(cluster) %>%
  slice_min(distance, n = 1, with_ties = FALSE) %>%
  ungroup()

# use ID for cluster Name
closest_samples$clusterName <- paste0('cluster_', closest_samples$ID)

# reformat
clusterNames <- closest_samples %>% 
  dplyr::select(cluster, clusterName)

# add to PCA table
clusterTab <- tsne_plot %>% 
  left_join(., clusterNames) %>% 
  dplyr::select(ID, Isolate, clusterName)

# update metatab 
supplementaryTab <- left_join(metaTab, clusterTab)

# add duplicated cluster information
for (NR in 1:nrow(duplicatedIDs)){
  # get single (with cluster info)
  originale <- duplicatedIDs$same_as_row[NR]
  # get other
  copia <- duplicatedIDs$duplicated_row[NR]
  # supplementary
  supplementaryTab[supplementaryTab$ID == copia, 'clusterName'] <- supplementaryTab[supplementaryTab$ID == originale, 'clusterName']
}


# Table values -----------------------------------------

# get values for a classification overview: Tab 1

# table for each continent
supplementaryTab %>% 
  mutate(Country = ifelse(Country == 'Europe & Asia', 'others', Country)) %>% 
  mutate(Country = ifelse(Country == 'South America', 'others', Country)) %>% 
  mutate(Country = ifelse(Country == 'Oceania', 'others', Country)) %>% 
  mutate(Country = ifelse(is.na(Country), 'others', Country)) %>% 
  group_by(Country) %>% 
  summarise(n = n(),
            perCent_classified_genoflu = round(sum(!is.na(US_class))/n*100, digits = 2),
            perCent_classified_genin2 = round(sum(!is.na(EU_class))/n*100, digits = 2),
            perCent_classified_inHouse = round(sum(!is.na(clusterName))/n*100, digits = 2),
            count_genoflu = length(na.omit(US_class)),
            count_genin2 = length(na.omit(EU_class)),
            count_inHouse = length(na.omit(clusterName))) %>% 
  fwrite(., file = 'MWy/forPaper_AIV/Table1_country.csv')


# last row for total
supplementaryTab %>% 
  summarise(n = n(),
            perCent_classified_genoflu = round(sum(!is.na(US_class))/n*100, digits = 2),
            perCent_classified_genin2 = round(sum(!is.na(EU_class))/n*100, digits = 2),
            perCent_classified_inHouse = round(sum(!is.na(clusterName))/n*100, digits = 2)) %>% 
  fwrite(., file = 'MWy/forPaper_AIV/Table1_total.csv')


# number of identified cluster in house
supplementaryTab$clusterName %>% unique() %>% na.omit() %>% length()
# number of identified cluster US
supplementaryTab$US_class %>% unique() %>% na.omit() %>% length()
# number of identified cluster EU
supplementaryTab$EU_class %>% unique() %>% na.omit() %>% length()

# number of new clusters not identified by any other
newIdentifiedClusters <- c()
for (CL in na.omit(unique(supplementaryTab$clusterName))) {
  miniTab <- supplementaryTab %>% 
    filter(clusterName == CL)
  classifiedRow <- miniTab %>% 
    filter(is.na(US_class) & is.na(EU_class)) %>% nrow()
  if (nrow(miniTab) == classifiedRow){
    newIdentifiedClusters <- c(newIdentifiedClusters, CL)
  }
  
}
length(newIdentifiedClusters)

onlyNew <- unique(tsne_plot[is.na(tsne_plot$col) & is.na(tsne_plot$USA), 'newCluster'])
for (CL in onlyNew){
  # table subset
  subTable <- tsne_plot[tsne_plot$newCluster == CL, ]
  NRsubTable <- nrow(subTable)
  # count unclassified
  NRunclassified <- nrow(unique(tsne_plot[is.na(tsne_plot$col) & is.na(tsne_plot$USA) & tsne_plot$newCluster == CL, ]))
  # remove if less then 5% not covered
  amountUnclassified <- NRunclassified/NRsubTable*100
  if (amountUnclassified < 5){
    onlyNew <- onlyNew[onlyNew != CL]
  }
}
length(onlyNew )