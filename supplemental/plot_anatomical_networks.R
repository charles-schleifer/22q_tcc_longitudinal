# C. Schleifer 8/4/2023
# plot Huang et al 2021 networks derived from FreeSurfer regions in single example subject

# clear environment
rm(list = ls(all.names = TRUE))

# list of packages to load
packages <- c("optparse","ciftiTools", "dplyr", "tidyr", "magrittr", "DescTools","parallel")

# Install packages not yet installed
# Note: ciftiTools install fails if R is started without enough memory
installed_packages <- packages %in% rownames(installed.packages())
# comment out line below to skip install
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# set up workbench
wbpath <- "/Applications/workbench/bin_macosx64/"
#wbpath <- "/u/project/cbearden/data/scripts/tools/workbench/bin_rh_linux64/wb_command"
ciftiTools.setOption("wb_path", wbpath)

# set up hoffman path
#hoffman <- "/u/project/cbearden/data/"
hoffman <- "~/Desktop/hoffman_mount/"


sesh="Q_0001_09242012"
sessions_dir <- file.path(hoffman,"22q/qunex_studyfolder/sessions")
bold_name_use="resting"
after_dir="/images/functional/"
file_end="_Atlas_s_hpss_res-mVWMWB1d_lpss.dtseries.nii"

## load parcellation data
#ji_path <- file.path(hoffman,"/22q/qunex_studyfolder/analysis/fcMRI/roi/ColeAnticevicNetPartition-master")
#ji_key <- read.table(file.path(ji_path,"/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"),header=T)
#ji_net_keys <- ji_key[,c("NETWORKKEY","NETWORK")] %>% distinct %>% arrange(NETWORKKEY)
#print(ji_net_keys)
#xii_Ji_network <- read_cifti(file.path(ji_path,"/data/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dscalar.nii"), brainstructures = "all")

# read subject-level thalamic atlas cifti
#xii_subcort_structs <- read_cifti(file.path(hoffman,"22q/qunex_studyfolder/analysis/fcMRI/roi/structures.dtseries.nii"), brainstructures = "all")
#sesh="Q_0001_09242012"
thal_path <- file.path(sessions_dir,sesh,"hcp",sesh,"T1w",sesh,"mri/ThalamicNuclei_Atlas_2mm.dscalar.nii")
xii_thal_orig <- read_cifti(thal_path, brainstructures = "all")
#view_xifti_volume(xii_thal_orig)

# read subject-level cortical DK atlas cifti (32K mesh)
dk_name <- paste0(sesh,".aparc.32k_fs_LR.dlabel.nii")
#dk_path <- file.path(hoffman,"22q/qunex_studyfolder/sessions",sesh,"hcp",sesh,"MNINonLinear/fsaverage_LR32k",dk_name)
dk_path <- file.path(sessions_dir,sesh,"hcp",sesh,"MNINonLinear/fsaverage_LR32k",dk_name)
xii_dk_orig <- read_cifti(dk_path, brainstructures = "all")
#view_xifti_surface(xii_dk_orig)

## group surface regions
# first get key with labels and values for DK dlabel
dk_key <- xii_dk_orig$meta$cifti$labels %>% as.data.frame
colnames(dk_key) <- c("Key","Red","Green","Blue","Alpha")
dk_key$label <- rownames(dk_key)

# match region names to Huang et al 2021 groupings
dk_groups <- list(prefrontal=c("superiorfrontal",
                           "caudalanteriorcingulate",
                           "rostralanteriorcingulate",
                           "medialorbitofrontal",
                           "lateralorbitofrontal",
                           "rostralmiddlefrontal",
                           "parsopercularis",
                           "parsorbitalis",
                           "parstriangularis"), 
              parietal=c("inferiorparietal",
                         "superiorparietal",
                         "precuneus",
                         "isthmuscingulate",
                         "posteriorcingulate",
                         "supramarginal"),
              temporal=c("superiortemporal",
                         "transversetemporal",
                         "middletemporal",
                         "fusiform",
                         "inferiortemporal",
                         "parahippocampal",
                         "entorhinal"),
              motor=c("caudalmiddlefrontal",
                      "paracentral",
                      "precentral"),
              somatosensory="postcentral",
              visual=c("pericalcarine",
                       "lateraloccipital",
                       "lingual",
                       "cuneus"))


# get groupings for left and right hemispheres
dk_left <- lapply(dk_groups, function(x) paste0("L_",x))
names(dk_left) <- paste0("L_",names(dk_left))
dk_right <- lapply(dk_groups, function(x) paste0("R_",x))
names(dk_right) <- paste0("R_",names(dk_right))
dk_all <- c(dk_left, dk_right)

# all region names
dk_all_regions <- unlist(dk_all) %>% as.vector %>% sort

# merge names to manually check matching spelling and missing regions
dk_check <- merge(x=dk_key, y=data.frame(label=dk_all_regions, name=dk_all_regions, exists="x"), by="label", all.x=TRUE, all.y=TRUE)

# get dscalar atlas from dlabel
xii_dk_scalar <- convert_xifti(xii_dk_orig, to="dscalar")

# get scalar values corresponding to each group of regions
dk_group_nums_l <- lapply(dk_left, function(x) dk_key[x,"Key"])
dk_group_nums_r <- lapply(dk_right, function(x) dk_key[x,"Key"])

# get cifti vertices corresponding to each group of regions 
dk_vertex_groups_l <- lapply(dk_group_nums_l, function(x)which(xii_dk_scalar$data$cortex_left %in% x))
dk_vertex_groups_r <- lapply(dk_group_nums_r, function(x)which(xii_dk_scalar$data$cortex_right %in% x))

# create new blank cifti with only zeros or NaNs
xii_blank <- xii_dk_scalar/xii_dk_scalar-1

# create cifti relabeled by region groups
xii_group <- xii_blank
ngroup <- length(dk_vertex_groups_l)
new_roi_ids <- NULL
new_roi_names <- NULL
for (i in 1:ngroup){
  # get different numbers for left and right IDs
  n <- plyr::round_any(ngroup+1, accuracy=10, f=ceiling)
  i_l <- i
  i_r <- n+i
  # update roi key left
  new_roi_ids <- c(new_roi_ids,i_l)
  new_roi_names <- c(new_roi_names,names(dk_vertex_groups_l[i]))
  # update roi key right
  new_roi_ids <- c(new_roi_ids,i_r)
  new_roi_names <- c(new_roi_names,names(dk_vertex_groups_r[i]))
  # get cifti row indices to change
  rows_l <- dk_vertex_groups_l[[i]]
  rows_r <- dk_vertex_groups_r[[i]]
  # edit cifti
  xii_group$data$cortex_left[rows_l,] <- i_l
  xii_group$data$cortex_right[rows_r,] <- i_r
}

# make new roi key
new_roi_key <- data.frame(id=new_roi_ids, name=new_roi_names)
out_path_key <- file.path(sessions_dir, sesh, after_dir, "anatomical_roi_key.csv")
#write.table(new_roi_key, file=out_path_key, col.names=T, row.names=F, quote=F, na="NA", sep=",", eol = "\n")

# set unused regions to NA to finalize new atlas
l_unused <- which(!xii_group$data$cortex_left %in% new_roi_key$id)
xii_group$data$cortex_left[l_unused,] <- NA 
r_unused <- which(!xii_group$data$cortex_right %in% new_roi_key$id)
xii_group$data$cortex_right[r_unused,] <- NA 

#view_xifti_surface(xii_group)

# match thal IDs to Huang et al 2021 groupings
# These networks were (1) prefrontal–mediodorsal (2) motor–ventral lateral (3) somatosensory–ventral posterolateral (4) temporal–medial geniculate (5) parietal–pulvinar (6) occipital–lateral geniculate (7) hippocampus–anterior nuclear groups
thal_groups <- list(prefrontal=c(12,13),
                  parietal=c(20,21,22,23),
                  temporal=09,
                  motor=c(26,27,28,29,30),
                  somatosensory=33,
                  visual=15)

thal_left <- lapply(thal_groups, function(x) 8100+x)
names(thal_left) <- paste0("L_", names(thal_left))
thal_right <- lapply(thal_groups, function(x) 8200+x)
names(thal_right) <- paste0("R_", names(thal_right))

# get cifti vertices corresponding to each group of regions 
thal_vertex_groups_l <- lapply(thal_left, function(x)which(xii_thal_orig$data$subcort %in% x))
thal_vertex_groups_r <- lapply(thal_right, function(x)which(xii_thal_orig$data$subcort %in% x))

# create new blank cifti with only zeros or NaNs
xii_blank_thal <- xii_thal_orig/xii_thal_orig-1

# create cifti relabeled by region groups
xii_thal_new <- xii_blank_thal
ngroup_thal <- length(thal_vertex_groups_l)
new_roi_ids_thal <- NULL
new_roi_names_thal <- NULL
for (i in 1:ngroup_thal){
  # get different numbers for left and right IDs
  n <- plyr::round_any(ngroup_thal+1, accuracy=10, f=ceiling)
  i_l <- i
  i_r <- n+i
  # update roi key left
  new_roi_ids_thal <- c(new_roi_ids_thal,i_l)
  new_roi_names_thal <- c(new_roi_names_thal,names(thal_vertex_groups_l[i]))
  # update roi key right
  new_roi_ids_thal <- c(new_roi_ids_thal,i_r)
  new_roi_names_thal <- c(new_roi_names_thal,names(thal_vertex_groups_r[i]))
  # get cifti row indices to change
  rows_l <- thal_vertex_groups_l[[i]]
  rows_r <- thal_vertex_groups_r[[i]]
  # edit cifti
  xii_thal_new$data$subcort[rows_l,] <- i_l
  xii_thal_new$data$subcort[rows_r,] <- i_r
}

# make new roi key
new_roi_key_thal <- data.frame(id=new_roi_ids_thal, name=new_roi_names_thal)

# set unused regions to NA to finalize new atlas
thal_unused <- which(!xii_thal_new$data$subcort %in% new_roi_key_thal$id)
xii_thal_new$data$subcort[thal_unused,] <- NA 
xii_thal_new <- remove_xifti(xii_thal_new, remove=c("cortex_left","cortex_right"))

#view_xifti_volume(xii_thal_new, slices=seq(30,46, by=2))

# combine surface and subcort ciftis
xii_new_atlas <- combine_xifti(xii_group, xii_thal_new)
xii_new_atlas$data$subcort <- xii_thal_new$data$subcort

## plot
# make dlabel with color key
view_xifti_volume(xii_new_atlas, fname="~/Dropbox/github/22q_tcc_longitudinal/figures/anatomical/thal_anatomical_axial_34_42.pdf", bg="white",crop=FALSE,colors="turbo", slices=c(34,42),color_mode="qualitative")

view_xifti_surface(xii_new_atlas,fname="~/Dropbox/github/22q_tcc_longitudinal/figures/anatomical/surface_anatomical.pdf", colors="turbo", color_mode = "qualitative")

  
