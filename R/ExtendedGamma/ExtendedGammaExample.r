

rExtendedGamma(100, 10000000,1000000)
#rExtendedGamma(10000000, 90000000,9000000)# time consuming as alpha=100000 takes time to allocate a vector of that many size. We can make it efficient by writing the code in a better way

rExtendedGamma(10000000000, 10000000,-1000000)

alpha=1;a=400;b=+100
samlpe_g= rExtendedGamma(alpha, a,b)

for(iii in 1:100000){
  samlpe_g[iii]=rExtendedGamma(alpha,a,b)
}

par(mfrow=c(1,2))

plot(density(samlpe_g))
plot(function(x){x^(alpha-1)*exp(-a*x^2+b*x)},xlim=c(0,.4))



