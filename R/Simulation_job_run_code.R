
sim_job_run<-function(){
library(Rfast)
library(benchmarkme)
library(Bessel)
library(gsl)





dir_loc="/Users/subhadippal/Desktop/BCUP/Sim Circular/"
for(ii in 29:40){
  data=rvmf(n =1000, mu=c(1,0),k = 10)
MC_OBject=VONF_pD_MPG_INVCDF_DA_POSTERIOR(Y=data, MCSamplerSize=5000)
file_loc=paste0(dir_loc, "MC_OBject_data_kappa_10_id_", ii,".RData")
save(MC_OBject, file=file_loc)
print(paste0("data_", ii, "Saved."))
}
}







