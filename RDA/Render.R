
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
res_folder = "../Results/Rdata/RDA/CDF_freqtrun4_L&R/"
subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)

dir.create("../Results/Html/CDF_freqtrun4_L&R/", recursive = TRUE)
for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_our_L&R.Rmd", 
                    params = list(subj_name="func_20150410",
                                  # subj_name = subj_name,
                                  data_folder="../../Results/Rdata/RDA/CDF_freqtrun4_L&R/"
                                  # data_folder=paste0('../',res_folder)
                                  ),
                    output_file = paste0("./",
                                         # "../../Results/Html/CDF_freqtrun4_L&R/",
                                         subj_name,'.html')
                    )
}


