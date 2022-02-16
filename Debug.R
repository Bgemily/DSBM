
# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


# -------------------------------------------------------------------------

res = 3
L2_vec = c()

generated_network = do.call(generate_network2_v3, results[[res]]$network_param)
tmp=sum(apply((generated_network$pdf_true_array)^2,c(1,2),sum))


L2_vec = c(L2_vec, tmp)


plot(L2_vec)
