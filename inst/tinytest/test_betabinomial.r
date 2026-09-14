



if ( requireNamespace("tinytest", quietly=TRUE) ){
  #graphics.off()
  tinytest::test_package("rsaber")
}


library("rsaber")
library("tinytest")
#graphics.off()



#
#
# Testing dbetabinomial ------------------------------------------------------
#
#

## Compare to rmutil ------------------

NN <- 100
aa <- 1
bb <- 1

t01_q01a <- rsaber::dbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=bb )
t01_q01b <- rmutil::dbetabinom(   seq(0,NN) , size=NN, m=aa/(aa+bb), s=aa+bb )


aa <- 1
bb <- 5
t01_q02a <- rsaber::dbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=bb )
t01_q02b <- rmutil::dbetabinom(   seq(0,NN) , size=NN, m=aa/(aa+bb), s=aa+bb )


aa <- 1
bb <- 50
t01_q03a <- rsaber::dbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=bb )
t01_q03b <- rmutil::dbetabinom(   seq(0,NN) , size=NN, m=aa/(aa+bb), s=aa+bb )


expect_equal( round(t01_q01a,6), round(t01_q01b,6) )
expect_equal( round(t01_q02a,6), round(t01_q02b,6) )
expect_equal( round(t01_q03a,6), round(t01_q03b,6) )


## Compare to constants ------------------

NN <- 20
aa <- 1
bb <- 2

# ## For generating the original
# sol_gen <- rsaber::dbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=bb )
# paste( round(sol_gen,5), collapse=", " )

dreference_20_1_1 <- c( 
  0.04762, 0.04762, 0.04762, 0.04762, 0.04762, 0.04762, 0.04762,
  0.04762, 0.04762, 0.04762, 0.04762, 0.04762, 0.04762, 0.04762, 
  0.04762, 0.04762, 0.04762, 0.04762, 0.04762, 0.04762, 0.04762 
  )

dreference_20_1_2 <- c(
  0.09091, 0.08658, 0.08225, 0.07792, 0.07359, 0.06926, 0.06494, 
  0.06061, 0.05628, 0.05195, 0.04762, 0.04329, 0.03896, 0.03463, 
  0.0303, 0.02597, 0.02165, 0.01732, 0.01299, 0.00866, 0.00433
  )

dreference_20_1_50 <- c( 
  0.71429, 0.20704, 0.05785, 0.01554, 0.004, 0.00099, 0.00023, 
  5e-05, 1e-05, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0 
  )


t01_q04a <- rsaber::dbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=1 )
t01_q04b <- rsaber::dbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=2 )
t01_q04c <- rsaber::dbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=50 )


expect_true( 
  sum(abs(
    round(t01_q04a,5) - dreference_20_1_1
  )) <= 1e-5
)

expect_true( 
  sum(abs(
    round(t01_q04b,5) - dreference_20_1_2
  )) <= 1e-5
)

expect_true( 
  sum(abs(
    round(t01_q04c,5) - dreference_20_1_50
  )) <= 1e-5
)


#
#
# Testing pbetabinomial ------------------------------------------------------
#
#

## Compare to rmutil ------------------

NN <- 100
aa <- 1
bb <- 1
t02_q01a <- rsaber::pbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=bb )
t02_q01b <- rmutil::pbetabinom(   seq(0,NN) , size=NN, m=aa/(aa+bb), s=aa+bb )


aa <- 1
bb <- 5
t02_q02a <- rsaber::pbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=bb )
t02_q02b <- rmutil::pbetabinom(   seq(0,NN) , size=NN, m=aa/(aa+bb), s=aa+bb )


aa <- 1
bb <- 50
t02_q03a <- rsaber::pbetabinomial( seq(0,NN) , Nn=NN, a=aa, b=bb )
t02_q03b <- rmutil::pbetabinom(   seq(0,NN) , size=NN, m=aa/(aa+bb), s=aa+bb )



expect_equal( round(t02_q01a,6), round(t02_q01b,6) )
expect_equal( round(t02_q02a,6), round(t02_q02b,6) )
expect_equal( round(t02_q03a,6), round(t02_q03b,6) )


## Compare to constants ------------------





#
#
# Testing qbetabinomial -------------------------------------------------------
#
#

## Compare to rmutil ------------------

NN <- 100
aa <- 1
bb <- 1

t03_q01a <- rsaber::qbetabinomial( seq(0,1,0.10) , Nn=NN, a=aa, b=bb )
t03_q01b <- rmutil::qbetabinom(   seq(0,1,0.10) , size=NN, m=aa/(aa+bb), s=aa+bb )


NN <- 100
aa <- 1
bb <- 5

t03_q02a <- rsaber::qbetabinomial( seq(0,1,0.10) , Nn=NN, a=aa, b=bb )
t03_q02b <- rmutil::qbetabinom(   seq(0,1,0.10) , size=NN, m=aa/(aa+bb), s=aa+bb )

expect_equal( round(t03_q01a,6), round(t03_q01b,6) )
expect_equal( round(t03_q02a,6), round(t03_q02b,6) )


## Compare to constants ------------------
















