
### Ppsbm result -------
data_folder = "../Results/Rdata/RDA/ppsbm/"
subj_name_vec = list.files(data_folder, full.names = FALSE, recursive = FALSE)


for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_ppsbm.Rmd", 
                    params = list(subj_name=subj_name),
                    output_file = paste0("../../Results/Html/ppsbm/",subj_name,'.html'))
}

### Our method -------
res_name_vec = c("CDF_v8_rmv1+3_keeptop_Nrestart1_totaltime336")
# res_name_vec = c("CDF_v1_selectNclus_Nrestart1_totaltime336")
# res_name_vec = list.files("../Results/Rdata/RDA/", pattern="CDF_v3_*",
#                           full.names = FALSE, recursive = FALSE)[-(1:16)]
for (res_name in res_name_vec) {
  res_folder = paste0("../Results/Rdata/RDA_v3/",res_name,"/")
  subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)
  subj_name_vec = c("func_20150417")
  dir.create(paste0("../Results/Html_v2/",res_name,"/"), recursive = TRUE)
  for (subj_name in subj_name_vec) {
    rmarkdown::render("./RDA/Visualize_res_our_v2.Rmd", 
                      params = list(subj_name=subj_name,
                                    rmvtop=FALSE,
                                    # rmvtop=TRUE,
                                    data_folder=paste0('../',res_folder)),
                      output_file = paste0("../../Results/Html_v2/",res_name,"/",
                                           subj_name,'.html'))
  }
  
}


### Our method -------
res_name_vec = c("CDF_v9_rmv2+3_keeptop_Nrestart1_totaltime336")
# res_name_vec = list.files("../Results/Rdata/RDA_v3/", pattern="CDF_v9_*",
#                           full.names = FALSE, recursive = FALSE)
for (res_name in res_name_vec) {
  res_folder = paste0("../Results/Rdata/RDA_v3/",res_name,"/")
  subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)
  # subj_name_vec = c("func_20150417")
  dir.create(paste0("../Results/Html/",res_name,"/"), recursive = TRUE)
  for (subj_name in subj_name_vec) {
    rmarkdown::render("./RDA/Visualize_res_our.Rmd", 
                      params = list(subj_name=subj_name,
                                    data_folder=paste0('../',res_folder)),
                      output_file = paste0("../../Results/Html/",res_name,"/",
                                           subj_name,'.html'))
  }
}
  
