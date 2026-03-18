# header ----
grupa="Papilionoidea"
Grassland_index = TRUE  # Mainīt uz FALSE, ja nevajag apakšmapi
sugaEGV="COEPAM"
suga="COEPAM"
suga_versija=paste0(suga,"_v0b1")
home_range=500
piepules_radiuss_m=home_range*3
LULCfiltrs_cell=99.99
LULCfiltrs_hr=99.99
gada_filtrs=2016
grupas_mainigajiem="./SpeciesModels/00_EGVtables/EGVtable_AllSpecies_v3_01032026.xlsx"
noverojumu_tabula="./SpeciesModels/00_Observations/Observations_v11_AAvj_tikls.csv"
apaksgala_limitacija=0.1
vpi_slieksnis=1
extraeffort="./SpeciesModels/00_Observations/ObsEffort_Papilionoidea.csv"

suppressPackageStartupMessages(library(tidyverse))
svarosanai=readr::read_csv("./SpeciesModels/00_FilteringWeighting/NoverojumuFiltresanai_Novirzem_20260301.csv")
svarosana_suga=svarosanai %>% 
  filter(CODE==suga)
svarosanas_sugas=svarosanai %>% 
  filter(WeightingGroup == svarosana_suga$WeightingGroup)

biotopu_svars=svarosana_suga$HabitatWeight
sugu_svars=1-biotopu_svars

papildpiepulei=read_csv(extraeffort)
papildpiepulei=papildpiepulei %>% 
  mutate(AccumDW_hr=paste0("AccumDW_",home_range),
         AccumTCL_hr=paste0("AccumTCL_",home_range))

#papilddati=readr::read_csv("./SpeciesModels/00_Observations/extra_AmphibiaReptilia.csv")
#papildsuga=papilddati %>% 
#  filter(Group==grupa) %>% 
#  filter(CODE==suga) %>% 
#  filter(Year>=gada_filtrs)
#papildgrupa=papilddati %>% 
#  filter(Group==grupa) %>% 
#  filter(Year>=gada_filtrs)

#izsledzamie_egv=c("egv_093", # papildinajums
#                  "egv_150",
#                  "egv_423",
#                  "egv_424",
#                  "egv_425",
#                  "egv_426",
#                  "egv_427",
#                  "egv_255",
#                  "egv_256",
#                  "egv_257",
#                  "egv_258",
#                  "egv_259")



# JAUNS: Ceļu loģika

grassland_folder = "00GrasslandIndex" # Šeit ieraksti vēlamo mapes nosaukumu

if (Grassland_index) {
  base_path = paste0("./TestingScripts/JekaterinaButkevica/", grupa, "/", grassland_folder, "/", suga_versija, "/")
} else {
  base_path = paste0("./TestingScripts/JekaterinaButkevica/", grupa, "/", suga_versija, "/")
}




fs::dir_create(base_path)

# libs -----
suppressPackageStartupMessages(library(plotROC))
suppressPackageStartupMessages(library(ecospat))
suppressPackageStartupMessages(library(maxnet))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(sfarrow))
suppressPackageStartupMessages(library(usdm))
suppressPackageStartupMessages(library(maps))
suppressPackageStartupMessages(library(rasterVis))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(SDMtune))
suppressPackageStartupMessages(library(ENMeval))
suppressPackageStartupMessages(library(zeallot))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(blockCV))
suppressPackageStartupMessages(library(overlapping))

# vif -----

print("Sāku VIF: ")
print(Sys.time())

grupas_mainigie=read_excel(grupas_mainigajiem)
grupas_mainigie=grupas_mainigie #%>% 
#  filter(!(egv_layername %in% izsledzamie_egv))

mainigajiem2=grupas_mainigie %>% 
  pivot_longer(cols=9:length(names(grupas_mainigie)),
               names_to="sugas_kodi",
               values_to="atzimes") %>% 
  filter(sugas_kodi==sugaEGV) %>% 
  filter(!is.na(atzimes))

grupai=paste0("./RasterGrids_100m/2024/Scaled/",mainigajiem2$egv_filename)

mainigie=terra::rast(grupai)
set.seed(1);solots_vif=usdm::vifstep(mainigie,th=10,size=10000)

write_rds(solots_vif,
          paste0(base_path,
                 "VIFstep10k_",suga,".RDS"))
vifrez=as.data.frame(solots_vif@results)
mainigajiem3=left_join(mainigajiem2,vifrez,by=c("egv_layername"="Variables"))
openxlsx::write.xlsx(mainigajiem3,
                     paste0(base_path,
                            "VIFs_",suga,".xlsx"))

# spatial autocorrelation -----

print("Sāku autokorelācijas analīzi: ")
print(Sys.time())




mazmainigajiem=mainigajiem3 %>% 
  filter(!is.na(VIF))
mazmainigo_celi=paste0("./RasterGrids_100m/2024/Scaled/",
                       mazmainigajiem$egv_filename)
mazmainigie=rast(mazmainigo_celi)

set.seed(1);blokiem=blockCV::cv_spatial_autocor(r=mazmainigie,num_sample=30000)
write_rds(blokiem,paste0(base_path,
                         "blockCV_",suga,".RDS"))

# effort bias -----

print("Sāku piepūles raksturošanu: ")
print(Sys.time())


visinoverojumi=read_csv(noverojumu_tabula,
                        guess_max=225000)
merksuga=visinoverojumi %>% 
  filter(CODE==suga) %>% 
  filter(Year>=gada_filtrs) %>%
  mutate(DoY=lubridate::yday(Date)) %>% 
  filter(!is.na(DoY)) %>% 
  dplyr::select(Group,CODE,Year,DoY,LKS_X,LKS_Y,id100)
#merksuga=rbind(merksuga,papildsuga)
visagrupa=visinoverojumi %>% 
  filter(Year>=gada_filtrs) %>%
  filter(Group==grupa) %>% 
  filter(CODE %in% svarosanas_sugas$CODE) %>% 
  mutate(DoY=lubridate::yday(Date)) %>% 
  filter(!is.na(DoY)) %>% 
  group_by(CODE) %>% 
  mutate(skaits=n()) %>% 
  ungroup() %>% 
  filter(skaits>=2) %>% 
  dplyr::select(Group,CODE,Year,DoY,LKS_X,LKS_Y,id100)
papildgrupa=papildpiepulei %>% 
  filter(Year>=gada_filtrs) %>%
  filter(Group==grupa) %>% 
#  filter(CODE %in% svarosanas_sugas$CODE) %>% 
  mutate(DoY=lubridate::yday(Date)) %>% 
  filter(!is.na(DoY)) %>% 
  group_by(CODE) %>% 
  mutate(skaits=n()) %>% 
  ungroup() %>% 
  filter(skaits>=2) %>% 
  dplyr::select(Group,CODE,Year,DoY,LKS_X,LKS_Y,id100)
visagrupa=rbind(visagrupa,papildgrupa)

table(visagrupa$CODE)

kodi=levels(factor(visagrupa$CODE))
parklajumiem=data.frame(sugas=kodi,
                        parklajumi=NA)

for(i in seq_along(kodi)){
  kods=kodi[i]
  vektors=visagrupa %>% 
    filter(CODE==kods)
  saraksts=list(pirmais=merksuga$DoY,
                otrais=vektors$DoY)
  parklajums=overlap(saraksts, 
                     plot=FALSE,
                     pairsOverlap = TRUE,
                     type="2")
  parklajumiem$parklajumi[i]=parklajums$OV
}

parklajumiem=parklajumiem %>% 
  mutate(svarots=parklajumi/sum(parklajumi))

write_csv(parklajumiem,paste0(base_path,
                              "ActivityWeights_",suga,".csv"))


visagrupa2=visinoverojumi %>% 
  filter(Year>=gada_filtrs) %>%
  filter(Group==grupa) %>% 
  mutate(DoY=lubridate::yday(Date)) %>%
  filter(Mislocation_Sea==0) %>% 
  filter(Mislocation_CLC==0) %>% 
  filter(AccumDW_cell<=LULCfiltrs_cell) %>% 
  filter(AccumTCL_cell<=LULCfiltrs_cell) %>% 
  filter(AccumDW_hr<=LULCfiltrs_hr) %>% 
  filter(AccumTCL_hr<=LULCfiltrs_hr) %>% 
  filter(CODE %in% kodi) %>% 
  dplyr::select(Group,CODE,Year,DoY,LKS_X,LKS_Y,id100)
papildgrupa=papildpiepulei %>% 
  filter(Year>=gada_filtrs) %>%
  filter(Group==grupa) %>% 
  mutate(DoY=lubridate::yday(Date)) %>%
  filter(Mislocation_Sea==0) %>% 
  filter(Mislocation_CLC==0) %>% 
  filter(AccumDW_cell<=LULCfiltrs_cell) %>% 
  filter(AccumTCL_cell<=LULCfiltrs_cell) %>% 
  filter(AccumDW_hr<=LULCfiltrs_hr) %>% 
  filter(AccumTCL_hr<=LULCfiltrs_hr) %>% 
  filter(CODE %in% kodi) %>% 
  dplyr::select(Group,CODE,Year,DoY,LKS_X,LKS_Y,id100)
visagrupa2=rbind(visagrupa2,papildgrupa)

kodi=levels(factor(visagrupa2$CODE))
template100=terra::rast("./Templates/TemplateRasters/LV100m_10km.tif")

sasvaroti_noverojumi=visagrupa2 %>% 
  left_join(parklajumiem,by=c("CODE"="sugas")) %>% 
  group_by(id100) %>% 
  summarise(svars=sum(svarots)) %>% 
  ungroup()

punkti100=sfarrow::st_read_parquet("./Templates/TemplateGridPoints/pts100_sauzeme.parquet")
punkti100=st_transform(punkti100,crs=3059)

noverojumi100=punkti100 %>% 
  right_join(sasvaroti_noverojumi,by=c("id"="id100"))

source("./RScripts_final/sdmhelpers_sfKDEterraff.R")
sugu_piepule=sfKDEterraff(x=noverojumi100,
                                weight_field = "svars",
                                sigma=piepules_radiuss_m,
                                ref = template100,
                                normalize="pdf",
                                mask=TRUE)
sugu_piepule=sugu_piepule*sugu_svars

# biotopi

#biotopi=sfarrow::st_read_parquet("./SpeciesModels/00_effortEUhabitats/EUHabitats_points.parquet")

#atlasiti_biotopi=biotopi %>% 
#  mutate(grupas_svars=case_when(EUhabitat=="atsegumi"~svarosana_suga$Atsegumi,
#                                EUhabitat=="kapas"~svarosana_suga$Kapas,
#                                EUhabitat=="koki"~svarosana_suga$Koki,
#                                EUhabitat=="purvi"~svarosana_suga$Purvi,
#                                EUhabitat=="udeni"~svarosana_suga$Udeni,
#                                EUhabitat=="virsaji"~svarosana_suga$Virsaji,
#                                EUhabitat=="zalaji"~svarosana_suga$Zalaji,
#                                .default=NA)) %>% 
#  filter(grupas_svars>0)


#merksuga=visinoverojumi %>% 
#  filter(CODE==suga) %>% 
#  filter(Year>=gada_filtrs) %>%
#  mutate(DoY=lubridate::yday(Date)) %>% 
#  filter(!is.na(DoY))


#source("./RScripts_final/sdmhelpers_DoYKDEweights.R")
#svari_biotopiem=doy_kde_weights(A_dates = merksuga$Date,
#                                B_dates=atlasiti_biotopi$OBS_DATE,
#                                scale="minmax_A")
#atlasiti_biotopi$sezonas_svars=svari_biotopiem$weights

#svaroti_biotopi=data.frame(atlasiti_biotopi) %>% 
#  mutate(kopejais_svars=grupas_svars*sezonas_svars) %>% 
#  filter(kopejais_svars>0) %>% 
#  group_by(id,x,y) %>% 
#  summarise(beigu_svars=sum(kopejais_svars)) %>% 
#  ungroup()
#svaroti_biotopi_sf=st_as_sf(svaroti_biotopi,coords=c("x","y"),remove = FALSE,crs=3059)

#biotopu_piepule=sfKDEterraff(x=svaroti_biotopi_sf,
#                             weight_field = "beigu_svars",
#                             sigma=piepules_radiuss_m,
#                             ref = template100,
#                             normalize="pdf",
#                             mask=TRUE)
#biotopu_piepule=biotopu_piepule*biotopu_svars

#svarota_piepule=biotopu_piepule+sugu_piepule
svarota_piepule=sugu_piepule
writeRaster(svarota_piepule,
            filename=paste0(base_path,
                            "EffortBiasLayer_",suga,".tif"),
            overwrite=TRUE)

svarota_piepule=terra::rast(paste0(base_path,
                                   "EffortBiasLayer_",suga,".tif"))





# TrainValTest -----

print("Sāku veidot datu kopas: ")
print(Sys.time())



klatbutnes=visinoverojumi %>% 
  filter(CODE==suga) %>% 
  filter(Year>=gada_filtrs) %>%
  filter(Mislocation_Sea==0) %>% 
  filter(Mislocation_CLC==0) %>% 
  filter(AccumDW_cell<=LULCfiltrs_cell) %>% 
  filter(AccumTCL_cell<=LULCfiltrs_cell) %>% 
  filter(AccumDW_hr<=LULCfiltrs_hr) %>% 
  filter(AccumTCL_hr<=LULCfiltrs_hr) %>% 
  filter(!duplicated(id100)) %>% 
  dplyr::select(CODE,LKS_X,LKS_Y)


videjais=as.numeric(terra::global(svarota_piepule,"mean",na.rm=TRUE))

cita_piepule3=ifel(svarota_piepule<=videjais*apaksgala_limitacija,videjais*apaksgala_limitacija,svarota_piepule)
set.seed(1);fonavietas3=terra::spatSample(cita_piepule3, 
                              size = 40000, 
                              na.rm = TRUE, 
                              values = TRUE, 
                              xy = TRUE, 
                              method="weights") |> as.data.frame()
#ggplot(fonavietas3,aes(x,y))+geom_point(alpha=0.25)+labs(title="fonavietas3")

fonavietas3=fonavietas3 %>% 
  mutate(CODE=suga,
         LKS_X=x,
         LKS_Y=y) %>% 
  dplyr::select(CODE,LKS_X,LKS_Y)


kvadrati=blokiem$plots$map_plot$data

klatbutnes_sf=st_as_sf(klatbutnes,coords=c("LKS_X","LKS_Y"),remove=FALSE,crs=3059)
klatbutnes_sf$id=rownames(klatbutnes_sf)
klatbutnes_tikla=st_join(klatbutnes_sf,kvadrati)
klatbutnes_tikla=klatbutnes_tikla %>% 
  filter(!duplicated(id))
klatbutnes_tikla3=klatbutnes_tikla


tikls100=sfarrow::st_read_parquet("./Templates/TemplateGrids/tikls100_sauzeme_3059.parquet")

fonavietas3_sf=st_as_sf(fonavietas3,coords=c("LKS_X","LKS_Y"),remove=FALSE,crs=3059)
fonavietas3_sf$id=rownames(fonavietas3_sf)
fonavietas3_tikla=st_join(fonavietas3_sf,kvadrati)
fonavietas3_tikla=fonavietas3_tikla %>% 
  filter(!duplicated(id)) 

fonavietas3_tikla2=st_join(fonavietas3_tikla,tikls100)
fonavietas3_tikla2=fonavietas3_tikla2 %>% 
  filter(!duplicated(id.y))
fonavietas3_tikla3=fonavietas3_tikla2

source("./RScripts_final/sdmhelpers_FoldSelectionFunction.R")
funcija1=select_joint_folds(pres_folds = klatbutnes_tikla3$folds,
                            bg_folds = fonavietas3_tikla3$folds,
                            target_prop     = 0.25,
                            pres_bounds     = c(0.20, 0.30),
                            bg_bounds     = c(0.20, 0.30),
                            max_tries       = 10000,
                            bg_match_weight = 1,
                            include_na      = FALSE,
                            seed            = 1,
                            print_report    = TRUE)

suga_testpres=klatbutnes_tikla3 %>% 
  filter(folds %in% funcija1$selected_fold_ids)
suga_testpres_df=data.frame(suga_testpres) %>% 
  mutate(x=LKS_X,
         y=LKS_Y) %>% 
  dplyr::select(x,y)
sfarrow::st_write_parquet(suga_testpres,
                          paste0(base_path,
                                 "PointsPresTest_",suga,".parquet"))
suga_testbg=fonavietas3_tikla3 %>% 
  filter(folds %in% funcija1$selected_fold_ids)
suga_testbg_df=data.frame(suga_testbg) %>% 
  mutate(x=LKS_X,
         y=LKS_Y) %>% 
  dplyr::select(x,y)
sfarrow::st_write_parquet(suga_testbg,
                          paste0(base_path,
                                 "PointsBgTest_",suga,".parquet"))


suga_trainpres=klatbutnes_tikla3 %>% 
  filter(!(folds %in% funcija1$selected_fold_ids))
suga_trainpres_df=data.frame(suga_trainpres) %>% 
  mutate(x=LKS_X,
         y=LKS_Y) %>% 
  dplyr::select(x,y)
sfarrow::st_write_parquet(suga_trainpres,
                          paste0(base_path,
                                 "PointsPresTrain_",suga,".parquet"))

suga_trainfons=fonavietas3_tikla3 %>% 
  filter(!(folds %in% funcija1$selected_fold_ids))
suga_trainfons_df=data.frame(suga_trainfons) %>% 
  mutate(x=LKS_X,
         y=LKS_Y) %>% 
  dplyr::select(x,y)
sfarrow::st_write_parquet(suga_trainfons,
                          paste0(base_path,
                                 "PointsBgTrain_",suga,".parquet"))



block_folds <- get.block(occ = suga_trainpres_df, 
                         bg = suga_trainfons_df)

write_rds(block_folds,
          paste0(base_path,
                 "blocks_",suga,".RDS"))


write_csv(suga_testpres_df,
          paste0(base_path,
                 "PointsPresTest_",suga,".csv"))
write_csv(suga_testbg_df,
          paste0(base_path,
                 "PointsBgTest_",suga,".csv"))
write_csv(suga_trainpres_df,
          paste0(base_path,
                 "PointsPresTrain_",suga,".csv"))
write_csv(suga_trainfons_df,
          paste0(base_path,
                 "PointsBgTrain_",suga,".csv"))

# parameterisation -----

print("Sāku parametrizāciju: ")
print(Sys.time())


objekti=data.frame(objekti=ls())
objekti2=objekti %>% 
  filter(!(objekti %in% c("grupa","suga","suga_versija","vpi_slieksnis", "base_path")))

rm(list=objekti2$objekti)
rm(objekti2)

videi=read_excel(paste0(base_path,
                        "VIFs_",suga,".xlsx"))
videi=videi %>% 
  filter(!is.na(VIF))
vide=terra::rast(paste0("./RasterGrids_100m/2024/Scaled/",
                        videi$egv_filename))

train_pres=read_csv(paste0(base_path,
                           "PointsPresTrain_",suga,".csv"))
train_bg=read_csv(paste0(base_path,
                         "PointsBgTrain_",suga,".csv"))
test_pres=read_csv(paste0(base_path,
                          "PointsPresTest_",suga,".csv"))
test_bg=read_csv(paste0(base_path,
                        "PointsBgTest_",suga,".csv"))

trenin_dati <- prepareSWD(species = suga,
                          p = train_pres,
                          a = train_bg,
                          env = vide)
trenin_dati=addSamplesToBg(trenin_dati)
block_folds <- get.block(occ = trenin_dati@coords[trenin_dati@pa == 1, ], 
                         bg = trenin_dati@coords[trenin_dati@pa == 0, ])



testa_dati=prepareSWD(species=suga,
                      p = test_pres,
                      a = test_bg,
                      env = vide)
testa_dati=addSamplesToBg(testa_dati)

klatbutnes=cbind(trenin_dati@coords[trenin_dati@pa == 1, ],
                 trenin_dati@data[trenin_dati@pa == 1, ])
fons=cbind(trenin_dati@coords[trenin_dati@pa == 0, ],
           trenin_dati@data[trenin_dati@pa == 0, ])


# bezvariabilitātes pazīmju izmešana
kraa=data.frame(standartnovirzes=apply(trenin_dati@data,2,sd))
kraa$nosaukums=rownames(kraa)
sliktie=which(kraa$standartnovirzes==0)
sliktie_egv=kraa$nosaukums[sliktie]
videi=videi %>% filter(!(egv_layername %in% sliktie_egv))

fons_kopa=rbind(testa_dati@coords[testa_dati@pa == 0, ],
                trenin_dati@coords[trenin_dati@pa == 0, ])
fonsX=data.frame(fons_kopa) %>% 
  mutate(x=X,
         y=Y) %>% 
  dplyr::select(x,y)
fons_SWD <- prepareSWD(species = suga, 
                       a = fonsX, 
                       env = vide)
kraa=data.frame(standartnovirzes=apply(fons_SWD@data,2,sd))
kraa$nosaukums=rownames(kraa)
sliktie=which(kraa$standartnovirzes==0)
sliktie_egv=kraa$nosaukums[sliktie]
videi=videi %>% filter(!(egv_layername %in% sliktie_egv))
vide=terra::rast(paste0("./RasterGrids_100m/2024/Scaled/",
                        videi$egv_filename))

set.seed(1);bg <- terra::spatSample(vide,
                                    size = 10000,
                                    method = "random",
                                    na.rm = TRUE,
                                    xy = TRUE,
                                    values = FALSE)
set.seed(1);bg <- prepareSWD(species = "Bgs",
                             a=bg,
                             env = vide)
kraa=data.frame(standartnovirzes=apply(bg@data,2,sd))
kraa$nosaukums=rownames(kraa)
sliktie=which(kraa$standartnovirzes==0)
sliktie_egv=kraa$nosaukums[sliktie]
videi=videi %>% filter(!(egv_layername %in% sliktie_egv))
vide=terra::rast(paste0("./RasterGrids_100m/2024/Scaled/",
                        videi$egv_filename))


trenin_dati <- prepareSWD(species = suga,
                          p = train_pres,
                          a = train_bg,
                          env = vide)
trenin_dati=addSamplesToBg(trenin_dati)



# no šejienes ----

pres_data <- trenin_dati@data[trenin_dati@pa == 1, , drop = FALSE]
pres_sd   <- apply(pres_data, 2, sd)

sliktie_egv_pres <- names(pres_sd)[pres_sd == 0]
sliktie_egv_pres

videi <- videi %>%
  filter(!(egv_layername %in% sliktie_egv_pres))

vide <- terra::rast(paste0(
  "./RasterGrids_100m/2024/Scaled/",
  videi$egv_filename
))

trenin_dati <- prepareSWD(
  species = suga,
  p = train_pres,
  a = train_bg,
  env = vide
)
trenin_dati <- addSamplesToBg(trenin_dati)

pres_df <- as.data.frame(trenin_dati@data[trenin_dati@pa == 1, , drop = FALSE])

sliktie_foldos <- c()

for (k in sort(unique(block_folds$occs.grp))) {
  train_idx <- block_folds$occs.grp != k
  this_sd <- apply(pres_df[train_idx, , drop = FALSE], 2, sd)
  bad <- names(this_sd)[this_sd == 0]
  sliktie_foldos <- union(sliktie_foldos, bad)
}

sliktie_foldos

videi <- videi %>%
  filter(!(egv_layername %in% sliktie_foldos))

vide <- terra::rast(paste0(
  "./RasterGrids_100m/2024/Scaled/",
  videi$egv_filename
))

trenin_dati <- prepareSWD(
  species = suga,
  p = train_pres,
  a = train_bg,
  env = vide
)
trenin_dati <- addSamplesToBg(trenin_dati)

# līdz šejienei ----



block_folds <- get.block(occ = trenin_dati@coords[trenin_dati@pa == 1, ], 
                         bg = trenin_dati@coords[trenin_dati@pa == 0, ])

testa_dati=prepareSWD(species=suga,
                      p = test_pres,
                      a = test_bg,
                      env = vide)
testa_dati=addSamplesToBg(testa_dati)

klatbutnes=cbind(trenin_dati@coords[trenin_dati@pa == 1, ],
                 trenin_dati@data[trenin_dati@pa == 1, ])
fons=cbind(trenin_dati@coords[trenin_dati@pa == 0, ],
           trenin_dati@data[trenin_dati@pa == 0, ])


#
source("./RScripts_final/sdmhelpers_TSS.R")

pirmais <- ENMevaluate(occs = klatbutnes, 
                       bg = fons,
                       algorithm = 'maxnet', 
                       partitions = 'user', 
                       user.grp=list(occs.grp=block_folds$occs.grp,
                                     bg.grp=block_folds$bg.grp),
                       tune.args = list(fc = c("L","Q","LQ"), 
                                        rm = c(0.2,0.25,1/3,0.5,1,2,3,4,5,10)),
                       other.settings=list(validation.bg="partition"),
                       user.eval = tss_user_eval)

labakajam=pirmais@results %>% 
  filter(cbi.val.avg==max(cbi.val.avg,na.rm=TRUE)) %>% 
  filter(tss.val.avg==max(tss.val.avg,na.rm=TRUE)) %>% 
  filter(auc.val.avg==max(auc.val.avg,na.rm=TRUE)) %>% 
  mutate(cbi.diff=abs(cbi.train-cbi.val.avg),
         tss.diff=abs(tss.train.avg-tss.val.avg)) %>% 
  filter(cbi.diff==min(cbi.diff,na.rm=TRUE)) %>% 
  filter(tss.diff==min(tss.diff,na.rm=TRUE)) %>% 
  filter(auc.diff.avg==min(auc.diff.avg,na.rm=TRUE)) %>% 
  mutate(fc_score=ifelse(fc=="L",1,
                         ifelse(fc=="Q",2,3))) %>% 
  filter(fc_score==min(fc_score,na.rm=TRUE))

fc_mazais=tolower(labakajam$fc)
fc_lielais=toupper(labakajam$fc)
rm=as.numeric(as.character(labakajam$rm))

write_rds(pirmais,
          paste0(base_path,
                 "ModelSelection_",suga,".RDS"))

# simplification correlations -----


print("Sāku vienkāršošanu sugas fona vietās: ")
print(Sys.time())

pirmais_modelis=train(method = "Maxnet",
                      data=trenin_dati,
                      folds=block_folds,
                      progress = TRUE,
                      fc=fc_mazais,
                      reg=rm)



fons_kopa=rbind(testa_dati@coords[testa_dati@pa == 0, ],
                trenin_dati@coords[trenin_dati@pa == 0, ])
fonsX=data.frame(fons_kopa) %>% 
  mutate(x=X,
         y=Y) %>% 
  dplyr::select(x,y)
fons_SWD <- prepareSWD(species = suga, 
                 a = fonsX, 
                 env = vide)



set.seed(1);korelejosais=SDMtune::varSel(pirmais_modelis,
                                         metric="tss",
                                         bg4cor = fons_SWD,
                                         cor_th = 0.8,
                                         permut=9,
                                         verbose=TRUE,
                                         progress=TRUE,
                                         interactive=FALSE)
write_rds(korelejosais,
          paste0(base_path,
                 "ModelVarSel_SpeciesBg_",suga,".RDS"))


print("Sāku vienkāršošanu nejaušās fona vietās: ")
print(Sys.time())


set.seed(1);bg <- terra::spatSample(vide,
                                    size = 10000,
                                    method = "random",
                                    na.rm = TRUE,
                                    xy = TRUE,
                                    values = FALSE)
set.seed(1);bg <- prepareSWD(species = "Bgs",
                             a=bg,
                             env = vide)
set.seed(1);korelejosais2=SDMtune::varSel(korelejosais,
                                         metric="tss",
                                         bg4cor = bg,
                                         cor_th = 0.8,
                                         permut=9,
                                         verbose=TRUE,
                                         progress=TRUE,
                                         interactive=FALSE)
write_rds(korelejosais2,
          paste0(base_path,
                 "ModelVarSel_BgBg_",suga,".RDS"))

# simplifications importance ------

print("Sāku vienkāršošanu pēc permutāciju ietekmes: ")
print(Sys.time())

set.seed(1);reducets=reduceVar(korelejosais2,
                               th=vpi_slieksnis,
                               metric="tss",
                               permut=9,
                               verbose=TRUE,
                               interactive=FALSE,
                               use_jk = TRUE)
write_rds(reducets,
          paste0(base_path,
                 "ModelVarReduce_",suga,".RDS"))
#reducets=read_rds(paste0("./SpeciesModels/",
#                         grupa,"/",
#                         suga_versija,"/",
#                         "ModelVarReduce_",suga,".RDS"))

# validation nulls ----

print("Sāku veidot vienkāršoto ENM: ")
print(Sys.time())

mazajai_videi=videi %>% 
  filter(egv_layername %in% names(reducets@data@data))
maza_vide=terra::rast(paste0("./RasterGrids_100m/2024/Scaled/",
                             mazajai_videi$egv_filename))

set.seed(1);mazas_vifvertibas=usdm::vif(maza_vide,
                            size=30000)
esosas_vifvertibas=read_excel(paste0(base_path,
                                     "VIFs_",suga,".xlsx"))
mazas_vifvertibas=mazas_vifvertibas %>% 
  mutate(VIFfinal=VIF) %>% 
  dplyr::select(-VIF)

esosas_vifvertibas=esosas_vifvertibas %>% 
  left_join(mazas_vifvertibas,by=c("egv_layername"="Variables"))

openxlsx::write.xlsx(esosas_vifvertibas,
                     paste0(base_path,
                            "VIFs_",suga,".xlsx"))


trenin_dati2 <- prepareSWD(species = suga,
                          p = train_pres,
                          a = train_bg,
                          env = maza_vide)
trenin_dati2=addSamplesToBg(trenin_dati2)
block_folds2 <- get.block(occ = trenin_dati2@coords[trenin_dati2@pa == 1, ], 
                         bg = trenin_dati2@coords[trenin_dati2@pa == 0, ])


testa_dati2=prepareSWD(species=suga,
                      p = test_pres,
                      a = test_bg,
                      env = maza_vide)
testa_dati2=addSamplesToBg(testa_dati2)




klatbutnes2=cbind(trenin_dati2@coords[trenin_dati2@pa == 1, ],
                  trenin_dati2@data[trenin_dati2@pa == 1, ])
fons2=cbind(trenin_dati2@coords[trenin_dati2@pa == 0, ],
            trenin_dati2@data[trenin_dati2@pa == 0, ])

isais_modelis <- ENMevaluate(occs = klatbutnes2[,1:2], 
                       bg = fons2[,1:2],
                       envs= maza_vide,
                       algorithm = 'maxnet', 
                       partitions = 'user', 
                       user.grp=list(occs.grp=block_folds2$occs.grp,
                                     bg.grp=block_folds2$bg.grp),
                       tune.args = list(fc = fc_lielais, 
                                        rm = rm),
                       other.settings=list(validation.bg="partition"),
                       user.eval = tss_user_eval)

write_rds(isais_modelis,
          paste0(base_path,
                 "Model_cvENM_",suga,".RDS"))


print("Sāku veidot validācijas nulles: ")
print(Sys.time())


nulles_modeli_val <- ENMnulls(isais_modelis, 
                              mod.settings = list(fc = fc_lielais, 
                                                  rm = rm), 
                              user.eval.type="kspatial",
                              no.iter = 100)

write_rds(nulles_modeli_val,
          paste0(base_path,
                 "Model_NullsCV_",suga,".RDS"))


# independent testing nulls -----


print("Sāku veidot testēšanas ENM: ")
print(Sys.time())


testesanas2=ENMevaluate(occs = klatbutnes2[,1:2], 
                        bg = fons2[,1:2],
                        envs= maza_vide,
                        occs.testing = test_pres,
                        algorithm = 'maxnet', 
                        partitions = 'testing', 
                        tune.args = list(fc = fc_lielais, 
                                         rm = rm))
write_rds(testesanas2,
          paste0(base_path,
                 "Model_testENM_",suga,".RDS"))



print("Sāku veidot testēšanas nulles: ")
print(Sys.time())




nulles_modeli_test <- ENMnulls(testesanas2, 
                              mod.settings = list(fc = fc_lielais, 
                                                  rm = rm), 
                              user.eval.type="testing",
                              no.iter = 100)
write_rds(nulles_modeli_test,
          paste0(base_path,
                 "Model_NullsTest_",suga,".RDS"))

#testesanas_tabula_apkopots=nulles_modeli_test@null.emp.results
#testesanas_tabula_nulles=nulles_modeli_test@null.results




# thresholds -----


print("Sāku veidot reducēto SDMtune: ")
print(Sys.time())


partaisits_reducetais=train(method = "Maxnet",
                            data=trenin_dati2,
                            folds=block_folds2,
                            progress = TRUE,
                            fc=fc_mazais,
                            reg=rm)

write_rds(partaisits_reducetais,
          paste0(base_path,
                 "Model_cvSDMtune_",suga,".RDS"))
#names(partaisits_reducetais@data@data)

#partaisits_reducetais=read_rds(paste0("./SpeciesModels/",
#                                      grupa,"/",
#                                      suga_versija,"/",
#                                      "Model_cvSDMtune_",suga,".RDS"))


print("Sāku veidot kombinēto SDMtune: ")
print(Sys.time())

labakais_comb=combineCV(partaisits_reducetais)
write_rds(labakais_comb,
          paste0(base_path,
                 "Model_combSDMtune_",suga,".RDS"))


print("Sāku sliekšņus: ")
print(Sys.time())

ths <- SDMtune::thresholds(labakais_comb, 
                           type = "cloglog",
                           test=testa_dati2)
#ths
write_csv(ths,
          paste0(base_path,
                 "thresholds_",suga,".csv"))




# projection map ------

print("Sāku HS projekciju: ")
print(Sys.time())


#map_best <- predict(labakais_comb,
#                    data = maza_vide,
#                    type = "cloglog",
#                    filename=paste0("./SpeciesModels/",
#                           grupa,"/",
#                           suga_versija,"/",
#                           "HSmap_",suga,".tif"),
#                    overwrite=TRUE)
#
tks_rezgis <- sfarrow::st_read_parquet("./Templates/TemplateGrids/tks93_50km.parquet")
tks_rezgis=st_transform(tks_rezgis,crs=3059)
fs::dir_create(paste0(base_path,
                      "HSmap_tiles"))
lapas <- levels(factor(tks_rezgis$NUMURS))
for(i in seq_along(lapas)){
  print(i)
  sakums=Sys.time()
  lapai=lapas[i]
  lapa=tks_rezgis %>%
    filter(NUMURS == lapai)
  egv_mazs=terra::crop(maza_vide,lapa)
  projekcijai=predict(labakais_comb,
                      data = egv_mazs,
                      type = "cloglog",
                      filename=paste0(base_path,
                                      "HSmap_tiles/",
                                      "HSmap_",suga,"_",lapai,".tif"),
                      overwrite=TRUE)
  #plot(projekcijai)
  beigas=Sys.time()
  ilgums=beigas-sakums
  print(ilgums)
}
slani=list.files(paste0(base_path,
                        "HSmap_tiles/"),full.names=TRUE)
virt_slanis=terra::vrt(slani)
slana_nosaukums=paste0("HSmap_",suga_versija)
names(virt_slanis)=slana_nosaukums
writeRaster(virt_slanis,
            paste0(base_path,
                   "HSmap_",suga,".tif"),
                   overwrite = TRUE)
map_best=rast(paste0(base_path,
                     "HSmap_",suga,".tif"))
#plot(projekcija)
#

unlink(paste0(base_path,
              "HSmap_tiles"),
       recursive=TRUE)

print("Sāku HS attēlu: ")
print(Sys.time())



slieksni=ths
slieksnis=slieksni[3,2]
slieksnis_vert=as.numeric(slieksnis)
apaksdala=mean(c(0,slieksnis_vert))
augsdala=mean(c(1,slieksnis_vert))

#map_best=terra::rast(paste0("./SpeciesModels/",
#                            grupa,"/",
#                            suga_versija,"/",
#                            "HSmap_",suga,".tif"))
slanis=map_best
names(slanis)="lyr1"
slanis_df=terra::as.data.frame(slanis,xy=TRUE)
slanis_df_augstie=slanis_df[slanis_df$lyr1>=slieksnis_vert,]
slanis_df_zemie=slanis_df[slanis_df$lyr1<slieksnis_vert,]

krasas <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
parejas <- c(0, apaksdala, slieksnis_vert, augsdala, 1)

karte=ggplot(slanis_df,aes(x=x,y=y,fill=lyr1))+
  geom_raster() +
  coord_fixed(ratio=1)+
  scale_fill_gradientn("",
                       colors = krasas,
                       values = parejas,
                       breaks=c(0,round(slieksnis_vert,3),1),
                       limits=c(0,1))+
  ggthemes::theme_map()+
  theme(legend.position = "inside",
        legend.position.inside=c(0,0.6),
        plot.background = element_rect(fill="white",color="white"))+
  labs(title=paste0(grupa,": ",suga))
karte


ggsave(karte,
       filename=paste0(base_path,
                       "PicHSmap_",
                       suga,
                       ".png"),
       width=900,height=550,units="px",dpi=100)


# stretch no 0 -> slieksnis -> 1 uz 0 -> 0.5 -> 1


print("Sāku HS stretch: ")
print(Sys.time())



objekts=terra::as.data.frame(slanis, cell=TRUE)
objekts2=objekts %>% 
  mutate(vert2=ifelse(lyr1>=slieksnis,
                      scales::rescale(objekts$lyr1,to=c(0.5,1),from=c(slieksnis,1)),
                      scales::rescale(objekts$lyr1,to=c(0,0.5),from=c(0,slieksnis))))

karte2=slanis
objekts3=objekts2[,c(1,3)]
karte2[objekts3$cell] <- objekts3[,-1]
#par(mfrow=c(1,2))
#terra::plot(slanis)
#terra::plot(karte2)
#par(mfrow=c(1,1))
terra::writeRaster(karte2,
                   filename=paste0(base_path,
                                   "HSstretch_",grupa,"_",suga_versija,".tif"),
                   overwrite=TRUE)



# marginal responses -----

print("Sāku marginālās atbildes: ")
print(Sys.time())


modelis_CV=partaisits_reducetais

#esosas_vifvertibas=read_excel(paste0("./SpeciesModels/",
#                                     grupa,"/",
#                                     suga_versija,"/",
#                                     "VIFs_",suga,".xlsx"))
#mazajai_videi=esosas_vifvertibas %>% 
#  filter(!is.na(PermImp_avg))
mainigo_tabula=mazajai_videi

mainigie=mainigo_tabula

augstums=ceiling(length(mainigie$longname_english)/7)*400

a=plotResponse(modelis_CV, 
               var = mainigie$egv_layername[1], 
               type = "cloglog",
               only_presence = TRUE,
               marginal = TRUE, 
               fun = median,
               rug = TRUE,
               col="black")
b=ggplot2::ggplot_build(a)

saistibas_funkcija=b$plot$data
saistibas_funkcija$nosaukums=mainigie$longname_english[1]
saistibas_funkcija$egv=mainigie$egv_layername[1]
vietas_presence=b$data[[3]]
vietas_presence$y=1
vietas_presence$nosaukums=mainigie$longname_english[1]
vietas_presence$egv=mainigie$egv_layername[1]
vietas_absence=b$data[[4]]
vietas_absence$y=-0.03
vietas_absence$nosaukums=mainigie$longname_english[1]
vietas_absence$egv=mainigie$egv_layername[1]

for(i in 2:length(mainigie$egv_layername)){
  a=plotResponse(modelis_CV, 
                 var = mainigie$egv_layername[i], 
                 type = "cloglog",
                 only_presence = TRUE,
                 marginal = TRUE, 
                 fun = median,
                 rug = TRUE,
                 col="black")
  b=ggplot2::ggplot_build(a)
  
  saistibas_funkcija_i=b$plot$data
  saistibas_funkcija_i$nosaukums=mainigie$longname_english[i]
  saistibas_funkcija_i$egv=mainigie$egv_layername[i]
  vietas_presence_i=b$data[[3]]
  vietas_presence_i$y=1
  vietas_presence_i$nosaukums=mainigie$longname_english[i]
  vietas_presence_i$egv=mainigie$egv_layername[i]
  vietas_absence_i=b$data[[4]]
  vietas_absence_i$y=-0.03
  vietas_absence_i$nosaukums=mainigie$longname_english[i]
  vietas_absence_i$egv=mainigie$egv_layername[i]
  
  saistibas_funkcija=bind_rows(saistibas_funkcija,saistibas_funkcija_i)
  vietas_presence=bind_rows(vietas_presence,vietas_presence_i)
  vietas_absence=bind_rows(vietas_absence,vietas_absence_i)
}


facet_levels <- saistibas_funkcija %>%
  distinct(nosaukums, egv) %>%
  mutate(egv_id = as.integer(sub("^egv_", "", egv))) %>%
  arrange(egv_id) %>%
  pull(nosaukums)
saistibas_funkcija <- saistibas_funkcija %>%
  mutate(nosaukums = factor(nosaukums, levels = facet_levels))
vietas_presence <- vietas_presence %>%
  mutate(nosaukums = factor(nosaukums, levels = facet_levels))
vietas_absence <- vietas_absence %>%
  mutate(nosaukums = factor(nosaukums, levels = facet_levels))

attels=ggplot(saistibas_funkcija)+
  geom_ribbon(data=saistibas_funkcija,aes(x=x,y=y,ymin=y_min,ymax=y_max),alpha=0.5)+
  geom_line(data=saistibas_funkcija,aes(x=x,y=y))+
  geom_point(data=vietas_presence,aes(x=x,y=y),size=0.5,alpha=0.5)+
  geom_point(data=vietas_absence,aes(x=x,y=y),size=0.5,alpha=0.5)+
  coord_cartesian(ylim=c(0,1))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank())+
  labs(y="Marginal response (cloglog)")+
  facet_wrap(~nosaukums,scales = "free_x",ncol=7,
             labeller = label_wrap_gen(width=25,multi_line = TRUE))
attels

ggsave(attels,
       filename=paste0(base_path,
                       "PicMargResp_",suga,".png"),
       width=2000,height=augstums,units="px",dpi=120)

# TSS ----

kopeja_tss=SDMtune::tss(labakais_comb)

tss_tabulai=data.frame(Group=grupa,
                       CODE=suga,
                       Version=suga_versija,
                       combTSS=kopeja_tss)
write_excel_csv(tss_tabulai,
                paste0(base_path,
                       "combTSS_",suga,".csv"))

# AUC ----


print("Sāku ROC: ")
print(Sys.time())



labakais_proc=SDMtune::plotROC(labakais_comb,test = testa_dati2)
labakais_proc

write_rds(labakais_proc,
          paste0(base_path,
                 "ROCcurves_",suga,".RDS"))
#labakais_proc=read_rds(paste0("./SpeciesModels/",
#                              grupa,"/",
#                              suga_versija,"/",
#                              "ROCcurves_",suga,".RDS"))
pic_roc=labakais_proc+
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=11),
        plot.subtitle = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=11),
        legend.position="inside",
        legend.position.inside = c(0.8,0.15))

ggsave(pic_roc,
       filename=paste0(base_path,
                       "PicROC_",suga,".png"),
       width=550,height=550,units="px",dpi=120)

# permutation importance -----


print("Sāku permutāciju ietekmi: ")
print(Sys.time())



vi_tss=varImp(partaisits_reducetais,
              permut=99)


mainigo_tabula=mazajai_videi %>% 
  left_join(mazas_vifvertibas,by=c("egv_layername"="Variables")) %>% 
  left_join(vi_tss,by=c("egv_layername"="Variable"))
mainigie=mainigo_tabula %>% 
  mutate(PermImp_avg=Permutation_importance,
         PermImp_sd=sd)

augstums=ceiling(length(mainigie$longname_latvian)/10)*300

pic_varimp=ggplot(mainigie,aes(x=reorder(longname_english,PermImp_avg),y=PermImp_avg))+
  geom_col()+
  geom_pointrange(data=mainigie,aes(x=reorder(longname_english,PermImp_avg),
                                    y=PermImp_avg,
                                    ymin=PermImp_avg-PermImp_sd,
                                    ymax=PermImp_avg+PermImp_sd))+
  scale_y_continuous("Permutation importance (%)",breaks=seq(0,100,10))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 75))+
  geom_text(data=mainigie,aes(x=reorder(longname_english,PermImp_avg),
                              y=103,
                              label=round(VIFfinal,3)),hjust=0,size=3)+
  coord_flip(ylim=c(0,105))+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        panel.grid.major.x = element_line(colour="grey"))
pic_varimp


ggsave(pic_varimp,
       filename=paste0(base_path,
                       "PicVarImpVIFs_",suga,".png"),
       width=1250,height=augstums,units="px",dpi=120)


mainigie_beigam=mainigie %>% 
  mutate(PermImp_avg=Permutation_importance,
         PermImp_sd=sd) %>% 
  dplyr::select(egv_layername,PermImp_avg,PermImp_sd)

pamata_tabula=read_excel(paste0(base_path,
                              "VIFs_",suga,".xlsx"))

beigu_tabula_mainigie=left_join(pamata_tabula,
                                mainigie_beigam,
                                by="egv_layername")

openxlsx::write.xlsx(beigu_tabula_mainigie,
                     paste0(base_path,
                            "VIFs_",suga,".xlsx"))

# parametrizacijas izveles attels -----


print("Sāku parametrizācijas attēlu: ")
print(Sys.time())



parametrizacijas_tabula=pirmais@results
openxlsx::write.xlsx(parametrizacijas_tabula,
                     paste0(base_path,
                            "ParametrizacijasTabula_",suga,".xlsx"))

parametrizacijas_tabula2=parametrizacijas_tabula %>% 
  dplyr::select(fc,rm,
                auc.train,auc.val.avg,
                cbi.train,cbi.val.avg,
                tss.train.avg,tss.val.avg,
                or.10p.avg,
                or.mtp.avg) %>% 
  pivot_longer(cols = auc.train:or.mtp.avg) %>% 
  mutate(nosaukums=case_when(name=="auc.train" ~ "AUC \n (training)",
                             name=="auc.val.avg" ~ "AUC \n (validation; \n average)",
                             name=="cbi.train" ~ "CBI \n (training)",
                             name=="cbi.val.avg" ~ "CBI \n (validation; \n average)",
                             name=="tss.train.avg" ~ "TSS \n (training)",
                             name=="tss.val.avg" ~ "TSS \n (validation; \n average)",
                             name=="or.10p.avg" ~ "Omission rate \n 10 % \n training",
                             name=="or.mtp.avg" ~ "Omission rate \n minimum \n training",
                             .default=NA))


rm_skaitlis=rm
parametrizacijas_tabula_labakajam=parametrizacijas_tabula2 %>% 
  filter(fc==fc_lielais) %>% 
  mutate(rm2=as.numeric(as.character(rm))) %>% 
  filter(rm2==rm_skaitlis)

selection_figure=ggplot()+
  theme_classic()+
  geom_point(data=parametrizacijas_tabula2,aes(nosaukums,value),
             col="grey",
             alpha=0.5,
             position = position_jitterdodge(jitter.width = 0.5,
                                             jitter.height = 0))+
  geom_point(data=parametrizacijas_tabula_labakajam,aes(nosaukums,value),
             size=2,
             position = position_jitterdodge(jitter.width = 0.5,
                                             jitter.height = 0))+
  scale_y_continuous("Value",
                     breaks=seq(0,1,by=0.1),
                     labels=scales::label_percent())+
  coord_cartesian(ylim=c(0,1))+
  theme(axis.title.x = element_blank())+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=11),
        plot.subtitle = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=11),
        legend.position="inside",
        legend.position.inside = c(0.8,0.15))

ggsave(selection_figure,
       filename=paste0(base_path,
                       "PicSelection_",suga,".png"),
       width=1250,height=550,units="px",dpi=120)


# nulles modelu atteli -----


## validacijai -----


print("Sāku validācijas nuļļu attēlu: ")
print(Sys.time())



nullem_val=read_rds(paste0(base_path,
                           "Model_NullsCV_",suga,".RDS"))
nullu_tabula_apkopots_val=nullem_val@null.emp.results
nullu_tabula_nulles_val=nullem_val@null.results

nullu_tabulai_val=evalplot.nulls(nullem_val, 
                             stats = c("auc.val","cbi.val","or.mtp","or.10p"), 
                             plot.type = "violin",
                             return.tbl = TRUE)
nulles_val=nullu_tabulai_val$null.avgs %>% 
  mutate(nosaukumi=case_when(metric=="auc.val"~"AUC\n(validation)",
                             metric=="cbi.val"~"CBI\n(validation)",
                             metric=="or.10p"~"Omission rate \n 10 % training\n(validation)",
                             metric=="or.mtp"~"Omission rate \n minimum training\n(validation)")) %>% 
  mutate(secibai=case_when(metric=="auc.val"~1,
                           metric=="cbi.val"~2,
                           metric=="or.10p"~3,
                           metric=="or.mtp"~4))
empiriskie_val=nullu_tabulai_val$empirical.results %>% 
  mutate(nosaukumi=case_when(metric=="auc.val"~"AUC\n(validation)",
                             metric=="cbi.val"~"CBI\n(validation)",
                             metric=="or.10p"~"Omission rate \n 10 % training\n(validation)",
                             metric=="or.mtp"~"Omission rate \n minimum training\n(validation)")) %>% 
  mutate(secibai=case_when(metric=="auc.val"~1,
                           metric=="cbi.val"~2,
                           metric=="or.10p"~3,
                           metric=="or.mtp"~4))

pic_nulles_val=ggplot()+
  geom_violin(data=nulles_val,
              aes(reorder(nosaukumi,secibai),avg),
              col="grey",
              fill="grey",
              alpha=0.2)+
  geom_jitter(data=nulles_val,
              aes(reorder(nosaukumi,secibai),avg),
              col="grey",
              width=0.2,
              shape=3)+
  geom_point(data=empiriskie_val,
             aes(reorder(nosaukumi,secibai),avg),
             col="black",
             size=3)+
  coord_cartesian(ylim=c(0,1))+
  theme_classic()+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     labels=scales::label_percent())+
  labs(y="Value")+
  theme(axis.title.x = element_blank())

pic_nulles_val


ggsave(pic_nulles_val,
       filename=paste0(base_path,
                       "PicNullsVal_",suga,".png"),
       width=700,height=350,units="px",dpi=120)

## neatkarigajai testesanai -----


print("Sāku testēšanas nuļļu attēlu: ")
print(Sys.time())



nullem_test=read_rds(paste0(base_path,
                            "Model_NullsTest_",suga,".RDS"))
nullu_tabula_apkopots_test=nullem_test@null.emp.results
nullu_tabula_nulles_test=nullem_test@null.results

nullu_tabulai_test=evalplot.nulls(nullem_test, 
                                 stats = c("auc.val","cbi.val","or.mtp","or.10p"), 
                                 plot.type = "violin",
                                 return.tbl = TRUE)
nulles_test=nullu_tabulai_test$null.avgs %>% 
  mutate(nosaukumi=case_when(metric=="auc.val"~"AUC\n(indep. testing)",
                             metric=="cbi.val"~"CBI\n(indep. testing)",
                             metric=="or.10p"~"Omission rate \n 10 % training\n(indep. testing)",
                             metric=="or.mtp"~"Omission rate \n minimum training\n(indep. testing)")) %>% 
  mutate(secibai=case_when(metric=="auc.val"~1,
                           metric=="cbi.val"~2,
                           metric=="or.10p"~3,
                           metric=="or.mtp"~4))
empiriskie_test=nullu_tabulai_test$empirical.results %>% 
  mutate(nosaukumi=case_when(metric=="auc.val"~"AUC\n(indep. testing)",
                             metric=="cbi.val"~"CBI\n(indep. testing)",
                             metric=="or.10p"~"Omission rate \n 10 % training\n(indep. testing)",
                             metric=="or.mtp"~"Omission rate \n minimum training\n(indep. testing)")) %>% 
  mutate(secibai=case_when(metric=="auc.val"~1,
                           metric=="cbi.val"~2,
                           metric=="or.10p"~3,
                           metric=="or.mtp"~4))

pic_nulles_test=ggplot()+
  geom_violin(data=nulles_test,
              aes(reorder(nosaukumi,secibai),avg),
              col="grey",
              fill="grey",
              alpha=0.2)+
  geom_jitter(data=nulles_test,
              aes(reorder(nosaukumi,secibai),avg),
              col="grey",
              width=0.2,
              shape=3)+
  geom_point(data=empiriskie_test,
             aes(reorder(nosaukumi,secibai),avg),
             col="black",
             size=3)+
  coord_cartesian(ylim=c(0,1))+
  theme_classic()+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     labels=scales::label_percent())+
  labs(y="Value")+
  theme(axis.title.x = element_blank())

pic_nulles_test


ggsave(pic_nulles_test,
       filename=paste0(base_path,
                       "PicNullsTest_",suga,".png"),
       width=700,height=350,units="px",dpi=120)
