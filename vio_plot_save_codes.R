

file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/ACF_Sim1_3D_MPG_DA_vioplot.pdf"
pdf(file=file, width = 6.5, height =3.5 )

data("ACF_Sim1_3D_MPG_DA_vioplot_data")
plot_Acf_violiin(Vio_ACF_MPG, size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()




file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/ACF_Sim1_3D_GEO_DA_vioplot.pdf"
pdf(file=file, width = 6.5, height =3.5 )

data(ACF_Sim1_3D_GEO_DA_vioplot_data)
plot_Acf_violiin(Vio_ACF_GEO, size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()






file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/ACF_Sim2_pD_MPG_DA_vioplot.pdf"
pdf(file=file, width = 6.5, height =3.5 )
data(ACF_sim2_pD_MPG_DA_vioplot_data)
plot_Acf_violiin(Vio_ACF_MPG_pd,size=.2, ifGray=TRUE, GridLines="Vertical")
dev.off()




file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/SLICE_ACF_vioplot.pdf"
#pdf(file=file, width = 6.5, height =3.5 )
data("ACF_Slice")
plot_Acf_violiin(ACF_Slice)
plt=plot_Acf_violiin(ACF_Slice,size=.2, ifGray=TRUE, GridLines="Vertical", Ylimit = c(.9, 1))
Sys.sleep(5)
ggsave(filename = file, plot=plt,width = 6.5, height =3.5)
#Sys.sleep(5)
#dev.off()



file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/Circular_MPG_ACF_vioplot1.pdf"

dd=circular_data_MPG_acf(25)
plt=plot_Acf_violiin(dd,size=.2, ifGray=TRUE, GridLines="Vertical", Ylimit = c(0, 1))
Sys.sleep(5)
ggsave(filename = file, plot=plt,width = 6.5, height =3.5)

tab=CalculateTable_MPG_circular(dd)
kappa_str=paste0("$\\kappa= ", c(1,5,10,15), "$")
tab1=cbind(kappa_str, tab$TableFormat1)
aa=apply(tab1,MARGIN = 1 ,FUN = function(x){ paste0( paste0(x, collapse = "&"), "\\\\ \n \\hline")})
cat(aa)





