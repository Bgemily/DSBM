#### Read data ####
# Use R.matlab to read the .mat file directly 
library(R.matlab)

# Read data: 


# We will calculate the dF/F following the method in Wan et al. 
# (Page 25, Step 6)
# dF/F = (F_cell - F_bsl)/(F_bsl-F_bkgd)
# F_bsl: 20th percentile in the nearby 61 time points
# F_bkgd is set to be zero 

calculate.dFF<-function(F.trace.full,bsl.prob,length.window){
  dFF.full=F.trace.full;
  if(is.nan(F.trace.full[1])){
    
  }else{
    
    for(i in 1:length(F.trace.full)){
      F.trace.ind= max(1,i-(length.window-1)/2):min(length(F.trace.full),i+(length.window-1)/2); 
      F.trace<-F.trace.full[F.trace.ind];
      F.bsl=quantile(F.trace,probs=bsl.prob);  
      dFF.full[i]= (F.trace.full[i]-F.bsl)/F.bsl;  
    }
    
  }
  return(dFF.full)
}

path.list=list.files('../Zebrafish_spinal_cord_development/AblationData/');

for(k in 1:length(path.list)){
  path=path.list[[k]]
  dat<- readMat(paste('../Zebrafish_spinal_cord_development/AblationData/',path,'/profile.mat',sep=''));
  
  dat.activity=dat$profile.all;
  n.neurons = dim(dat.activity)[1]; 
  n.timeframes=dim(dat.activity)[2];
  
  # Resolution: We then seamlessly transitioned to longitudinal functional imaging from 17.5 to 22 hpf of all post-mitotic neurons across 9-10 segments of the spinal cord at 4 Hz
  
  length.windows=61;
  
  # Obtain dFF for each trace
  dat.dFF<-dat.activity;
  for(j in 1:n.neurons){
    dat.dFF[j,]=calculate.dFF(dat.activity[j,],bsl.prob=0.2,length.window=61);
  }
  
  write.csv(dat.dFF,paste('../processed_AblationData/',path,'/dFF.csv',sep=''),col.names=F)
}




