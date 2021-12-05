
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
res_name_vec = c("CDF_freqtrun4_totaltime220",
                 "CDF_freqtrun4_totaltime240",
                 "CDF_freqtrun4_totaltime260",
                 "CDF_freqtrun4_totaltime280",
                 "CDF_freqtrun4_totaltime300",
                 "CDF_freqtrun5_totaltime220",
                 "CDF_freqtrun5_totaltime240",
                 "CDF_freqtrun5_totaltime260",
                 "CDF_freqtrun5_totaltime280",
                 "CDF_freqtrun5_totaltime300"
                 # "CDF_freqtrun6_totaltime220",
                 # "CDF_freqtrun6_totaltime240",
                 # "CDF_freqtrun6_totaltime260",
                 # "CDF_freqtrun6_totaltime280",
                 # "CDF_freqtrun6_totaltime300"
                 )
for (res_name in res_name_vec) {
  res_folder = paste0("../Results/Rdata/RDA/",res_name,"/")
  subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)
  
  dir.create(paste0("../Results/Html/",res_name,"/"), recursive = TRUE)
  for (subj_name in subj_name_vec) {
    rmarkdown::render("./RDA/Visualize_res_our.Rmd", 
                      params = list(subj_name=subj_name,
                                    data_folder=paste0('../',res_folder)),
                      output_file = paste0("../../Results/Html/",res_name,"/",
                                           subj_name,'.html'))
  }
  
}



