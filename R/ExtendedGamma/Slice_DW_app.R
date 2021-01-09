

file_path="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes"
setwd("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes")


for( dataid in 2:20){
  
    theta=rvmf(n=1000, mu=c(1,0) , k = 1)
    Mc=Dir_data_inf_DW_slice( theta, 50000 )
  
    MC_VONV=list(MC=Mc, data=theta, mu_true=c(1,0), kappa_true=1, method="Damian Walker Slice Sampler ")
    save(MC_VONV, file=paste0(file_path,"/WorkSpace_DW_slice/kappa_01/MC_kappa1_data_",dataid,".Rdata"))
print(paste0("data", dataid, "completed for kappa 1")) 
}



for( dataid in 1:20){
  
  theta=rvmf(n=1000, mu=c(1,0) , k = 5)
  Mc=Dir_data_inf_DW_slice( theta, 50000 )
  
  MC_VONV=list(MC=Mc, data=theta, mu_true=c(1,0), kappa_true=5, method="Damian Walker Slice Sampler ")
  save(MC_VONV, file=paste0(file_path,"/WorkSpace_DW_slice/kappa_5/MC_kappa5_data_",dataid,".Rdata"))
  print(paste0("data", dataid, "completed for kappa 5"))
}



for( dataid in 1:20){
  
  theta=rvmf(n=1000, mu=c(1,0) , k = 10)
  Mc=Dir_data_inf_DW_slice( theta, 50000 )
  
  MC_VONV=list(MC=Mc, data=theta, mu_true=c(1,0), kappa_true=10, method="Damian Walker Slice Sampler ")
  save(MC_VONV, file=paste0(file_path,"/WorkSpace_DW_slice/kappa_10/MC_kappa10_data_",dataid,".Rdata"))
  print(paste0("data", dataid, "completed for kappa 10"))
}




for( dataid in 1:20){
  
  theta=rvmf(n=1000, mu=c(1,0) , k = 15)
  Mc=Dir_data_inf_DW_slice( theta, 50000 )
  
  MC_VONV=list(MC=Mc, data=theta, mu_true=c(1,0), kappa_true=15, method="Damian Walker Slice Sampler ")
  save(MC_VONV, file=paste0(file_path,"/WorkSpace_DW_slice/kappa_15/MC_kappa15_data_",dataid,".Rdata"))
  print(paste0("data", dataid, "completed for kappa 15"))
}




