#!/usr/bin/env Rscript

# Case 2 ------------------------------------------------------------------

rm(list=ls())
file_path = "./functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)


library("optparse")

option_list = list(
  make_option(c("-n", "--NSim"), type="integer", default=1, help="number of repeated trials"),
  make_option(c("-k", "--N_overclus"), type="integer", default=5)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

NSim = opt$NSim
N_overclus = opt$N_overclus


Ncores = 20

# step_size = 0.02
SEED_vec = seq(189,110765,length.out=NSim)


library(foreach)
library(doParallel)
registerDoParallel(cores=Ncores)

results2 <- foreach(i = 1:NSim) %dopar%
{
  SEED = SEED_vec[i]
  main(case=2, SEED=SEED, N_clus=3, N_overclus=N_overclus, MaxIter = 10)
}


now = format(Sys.time(), "%Y%m%d_%H%M")
save.image(paste0('case2_NSim', NSim, '_Noverclus', N_overclus, '_', now, '.Rdata'))

# debug
SEED_vec = seq(189,110765,length.out=100)
i=2
SEED = SEED_vec[i]
main(case=2, SEED=SEED, N_clus=3, N_overclus=3, MaxIter = 5, bw=1)->r
tmp$clus_result$clusters
