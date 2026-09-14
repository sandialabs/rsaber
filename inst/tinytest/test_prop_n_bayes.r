



if ( requireNamespace("tinytest", quietly=TRUE) ){
  #graphics.off()
  tinytest::test_package("rsaber")
}


library("rsaber")
library("tinytest")
#graphics.off()


#
#
# Testing prop_n_bayes: finite n ------------------------------------------------------
#
#

vec_N <- c(20, 50, 100)
vec_p <- c(0.05, 0.10, 0.20)
vec_a <- c(0.01, 0.05, 0.10 )

N0 <- vec_N[2]
p0 <- vec_p[2]
a0 <- vec_a[3]


## Compare result to constant ---------------------------------

pnb01 <- rsaber::prop_n_bayes( x=0, N=vec_N, preq=p0, alpha=a0, prior=c(1,1) )$n
pnb02 <- rsaber::prop_n_bayes( x=0, N=N0, preq=vec_p, alpha=a0, prior=c(1,1) )$n
pnb03 <- rsaber::prop_n_bayes( x=0, N=N0, preq=p0, alpha=vec_a, prior=c(1,1) )$n

pnb04 <- rsaber::prop_n_bayes( x=1, N=vec_N, preq=p0, alpha=a0, prior=c(1,1) )$n


expect_equal(
  pnb01,
  c(10, 15, 18)
)

expect_equal(
  pnb02,
  c(26, 15, 8)
)

expect_equal(
  pnb03,
  c(25, 19, 15)
)

expect_equal(
  pnb04,
  c(16, 25, 30)
)





## Compare result to interval ---------------------------------



pib01a <- rsaber::prop_interval_bayes( x=0, n=pnb01, N=vec_N, alpha=a0, mesh="recycle" )
pib01b <- rsaber::prop_interval_bayes( x=0, n=pnb01-1, N=vec_N, alpha=a0, mesh="recycle" )

pib02a <- rsaber::prop_interval_bayes( x=0, n=pnb02, N=N0, alpha=a0, mesh="recycle" )
pib02b <- rsaber::prop_interval_bayes( x=0, n=pnb02-1, N=N0, alpha=a0, mesh="recycle" )

pib03a <- rsaber::prop_interval_bayes( x=0, n=pnb03, N=N0, alpha=vec_a, mesh="recycle" )
pib03b <- rsaber::prop_interval_bayes( x=0, n=pnb03-1, N=N0, alpha=vec_a, mesh="recycle" )

pib04a <- rsaber::prop_interval_bayes( x=1, n=pnb04, N=vec_N, alpha=a0, mesh="recycle" )
pib04b <- rsaber::prop_interval_bayes( x=1, n=pnb04-1, N=vec_N, alpha=a0, mesh="recycle" )



expect_true(  all(pib01a[["pucl"]] <= rep(p0, 3) ) )
expect_false( all(pib01b[["pucl"]] <= rep(p0, 3) ) )

expect_true(  all(pib02a[["pucl"]] <= vec_p ) )
expect_false( all(pib02b[["pucl"]] <= vec_p ) )

expect_true(  all(pib03a[["pucl"]] <= rep(a0, 3) ) )
expect_false( all(pib03b[["pucl"]] <= rep(a0, 3) ) )

expect_true(  all(pib04a[["pucl"]] <= rep(p0, 3) ) )
expect_false( all(pib04b[["pucl"]] <= rep(p0, 3) ) )




#
#
# Testing prop_n_bayes: infinite n ------------------------------------------------------
#
#


























