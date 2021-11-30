
### Ppsbm result -------
data_folder = "../Results/Rdata/RDA/ppsbm/"
subj_name_vec = list.files(data_folder, full.names = FALSE, recursive = FALSE)


for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_ppsbm.Rmd", 
                    params = list(subj_name=subj_name),
                    output_file = paste0("../../Results/Html/ppsbm/",subj_name,'.html'))
}

### Our method -------
res_folder = "../Results/Rdata/RDA/PDF+pairwise_v2/"
subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)

dir.create("../Results/Html/PDF+pairwise_v2/", recursive = TRUE)
for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_our.Rmd", 
                    params = list(subj_name=subj_name,
                                  data_folder=paste0('../',res_folder)),
                    output_file = paste0("../../Results/Html/PDF+pairwise_v2/",subj_name,'.html'))
}


### Our method -------
res_folder = "../Results/Rdata/RDA/CDF_Nclus3_randinit40/"
subj_name_vec = list.files(res_folder, full.names = FALSE, recursive = FALSE)

dir.create("../Results/Html/CDF_Nclus3_randinit40/", recursive = TRUE)
for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_our.Rmd", 
                    params = list(subj_name=subj_name,
                                  data_folder=paste0('../',res_folder)),
                    output_file = paste0("../../Results/Html/CDF_Nclus3_randinit40/",subj_name,'.html'))
}


