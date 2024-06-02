#Bifactor ESEM in R
library(lavaan)
library(semTools)
library(GPArotation)
library(psych)

#Simulate data (Lavaan)
#One general factor and 2 specific factors
set.seed(123456)
sim.model <- 'g=~.4*i1+.4*i2+.4*i3+.4*i4+.4*i5+.4*i6+.4*i7+.4*i8+.4*i9+.4*i10 + .4*i11+.4*i12+.4*i13+.4*i14+.4*i15

f1=~.5*i1+.5*i2+.5*i3+.5*i4+.5*i5
f2=~.6*i6+.6*i7+.6*i8+.6*i9+.6*i10
f3=~.6*i11+.6*i12+.6*i13+.6*i14+.6*i15

#orthogonal factors
g~~0*f1+0*f2 +0*f3
f1~~0*f2
f1~~0*f3
f2~~0*f3'
d <- simulateData(model=sim.model, model.type = "cfa", sample.nobs = 200)
datasim=d
datasim = ifelse(d>=0.5,1,0)
