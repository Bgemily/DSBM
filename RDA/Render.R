
### Ppsbm result -------
data_folder = "../Results/Rdata/RDA/ppsbm/"
subj_name_vec = list.files(data_folder, full.names = FALSE, recursive = FALSE)


for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_ppsbm.Rmd", 
                    params = list(subj_name=subj_name),
                    output_file = paste0("../../Results/Html/ppsbm/",subj_name,'.html'))
}

### Our method -------
res_folder = "../Results/Rdata/RDA/CDF_freqtrun4/"
subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)

dir.create("../Results/Html/CDF_freqtrun4/", recursive = TRUE)
for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_our.Rmd", 
                    params = list(subj_name=subj_name,
                                  data_folder=paste0('../',res_folder)),
                    output_file = paste0("../../Results/Html/CDF_freqtrun4/",subj_name,'.html'))
}

### Our method -------
res_name_vec = c("CDF_v14_keeptop_Nrestart20_freqtrun3_totaltime200")
# res_name_vec = list.files("../Results/Rdata/RDA/", pattern="CDF_v3_*",
#                           full.names = FALSE, recursive = FALSE)[-(1:16)]
for (res_name in res_name_vec) {
  res_folder = paste0("../Results/Rdata/RDA/",res_name,"/")
  subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)
  subj_name_vec = c("func_20150417")
  dir.create(paste0("../Results/Html_v2/",res_name,"/"), recursive = TRUE)
  for (subj_name in subj_name_vec) {
    rmarkdown::render("./RDA/Visualize_res_our_v2.Rmd", 
                      params = list(subj_name=subj_name,
                                    data_folder=paste0('../',res_folder)),
                      output_file = paste0("../../Results/Html_v2/",res_name,"/",
                                           subj_name,'.html'))
  }
  
}


### Our method -------
# res_name_vec = c("CDF_v12_rmvtop_Nrestart10_freqtrun3_totaltime200")
res_name_vec = list.files("../Results/Rdata/RDA/", pattern="CDF_v12_rmvtop*",
                          full.names = FALSE, recursive = FALSE)
for (res_name in res_name_vec) {
  res_folder = paste0("../Results/Rdata/RDA/",res_name,"/")
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
  
### Our method (restart) -------
  res_name_vec = c("CDF_v13_rmvtop_Nrestart20_freqtrun3_totaltime200",
                   "CDF_v14_keeptop_Nrestart20_freqtrun3_totaltime200")
  # res_name_vec = list.files("../Results/Rdata/RDA/", pattern="CDF_v3_*",
  #                           full.names = FALSE, recursive = FALSE)[-(1:16)]
  for (res_name in res_name_vec) {
    res_folder = paste0("../Results/Rdata/RDA/",res_name,"/")
    subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)
    # subj_name_vec = c("func_20150417")
    dir.create(paste0("../Results/Html/",res_name,"/"), recursive = TRUE)
    for (subj_name in subj_name_vec) {
      rmarkdown::render("./RDA/Visualize_res_our_debug.Rmd", 
                        params = list(subj_name=subj_name,
                                      data_folder=paste0('../',res_folder)),
                        output_file = paste0("../../Results/Html/",res_name,"/",
                                             subj_name,'.html'))
    }
    
}




