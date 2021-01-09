
plot(density(AllList_3D lst_MPG_pD[[1]][1:100]),main='MPG Augmentation')




acf_value=cbind( acf= acf(AllList_3D[[1]]$lst_GEO_3D$McSample_kappa,lag=10)$acf,lag=as.factor(1:11))
for(i in 2:20){
acf_value= rbind(acf_value,cbind( acf= acf(AllList_3D[[i]]$lst_GEO_3D$McSample_kappa,lag=10)$acf,lag=as.factor(1:11)))

}



plot(density(acf_value[3,]))

acf_value=data.frame(acf_value)
library(vioplot)
x1 <- acf_value$acf[acf_value$lag==2]
x2 <- acf_value$acf[acf_value$lag==3]
x3 <- acf_value$acf[acf_value$lag==4]
x4 <- acf_value$acf[acf_value$lag==5]
x5 <- acf_value$acf[acf_value$lag==6]
x6 <- acf_value$acf[acf_value$lag==7]
x7 <- acf_value$acf[acf_value$lag==8]
x8 <- acf_value$acf[acf_value$lag==9]
x9 <- acf_value$acf[acf_value$lag==10]
x10 <- acf_value$acf[acf_value$lag==10]

vioplot(x1, x2, x3,x4,x5,x6,x7,x8,x9,x10,  col="magenta")
title("Violin Plots of Miles Per Gallon")

library(ggplot2)
# Basic violin plot
p <- ggplot(acf_df, aes(factor(acf_df$lag),acf_df$acf)) + geom_violin()
p
# Rotate the violin plot
p + coord_flip()
# Set trim argument to FALSE
ggplot(ToothGrowth, aes(x=dose, y=len)) +    geom_violin(trim=FALSE)

ToothGrowth
