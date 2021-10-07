#!/usr/bin/env Rscript

# Case 3 ------------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)



library("optparse")

option_list = list(
  make_option(c("-n", "--NSim"), type="integer", default=1, 
              help="number of repeated trials")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


NSim = opt$NSim


# pp = FALSE
# NSim = 1000
Ncores = 20

# step_size = 0.02
SEED_vec = seq(139,8397,length.out=NSim)



library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results3 <- foreach(i = 1:NSim) %dopar% {
  SEED = SEED_vec[i]
  main(case=3, SEED=SEED, N_clus=3, N_overclus=5, MaxIter = 2)
}


now = format(Sys.time(), "%Y%m%d_%H%M")
save.image(paste0('case3_NSim', NSim, '_', now, '.Rdata'))




