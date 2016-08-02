context("Testing lmp and aovp")
library(lmPerm)

## Fixture
set.seed(485)
f1 = factor(round(runif(20,0,5)+0.5),1:5,labels=paste("F1",1:5))
f2 = factor(round(runif(20,0,5)+0.5),1:5,labels=paste("F2",1:5))
out = runif(20,1,20)
df = data.frame(f1,f2,out)


coeffs = structure(c(8.04386605673159, -3.46741918787981, 7.16350775379688, 
                     1.95094546703622, -6.8150014689813, -1.34869177614649, 
                     -1.75041950954279, -9.61050823992895, 9.62426867865337), 
                   .Names = c("(Intercept)", "f11", "f12", "f13", "f14", "f21", 
                              "f22", "f23", "f24"))

test_that("Basic lmp",{
    expect_output(mdl <- lmp(out ~ f1 + f2,data = df),"Settings:  unique SS ")
    expect_equal(length(mdl$coefficients),9)
    expect_equal(mdl$coefficients,coeffs)
})


test_that("Basic aovp",{
  
  expect_silent(mdl <- aovp(out ~ f1 + f2,data = df,settings=FALSE))
  expect_equal(length(mdl$coefficients),9)
  expect_equal(mdl$coefficients,coeffs)
})


