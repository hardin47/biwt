test_that("function runs", {
  expect_equal({set.seed(4747)
               samp.data <- MASS::mvrnorm(30, mu=c(0,0,0),
                                          Sigma=matrix(c(1,.75,-.75,.75,1,-.75,-.75,-.75,1),
                                                       ncol=3))

               r<-0.2 # breakdown

               biwt_cor(samp.data[,1:2], r=.2)$biwt_cor}, 0.6605525,
               tolerance = 0.01)
})
