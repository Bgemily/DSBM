
### Ppsbm result -------
data_folder = "../Results/Rdata/RDA/ppsbm/"
subj_name_vec = list.files(data_folder, full.names = FALSE, recursive = FALSE)


for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_ppsbm.Rmd", 
                    params = list(subj_name=subj_name),
                    output_file = paste0("../../Results/Html/ppsbm/",subj_name,'.html'))
}

### Our method -------
data_folder = "../Results/Rdata/RDA/PDF+pairwise/"
subj_name_vec = list.files(data_folder, full.names = FALSE, recursive = FALSE)

dir.create("../Results/Html/PDF+pairwise/", recursive = TRUE)
for (subj_name in subj_name_vec) {
  rmarkdown::render("./RDA/Visualize_res_our.Rmd", 
                    params = list(subj_name=subj_name),
                    output_file = paste0("../../Results/Html/PDF+pairwise/",subj_name,'.html'))
}


