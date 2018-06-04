context("Old unit tests")

parallel <- TRUE

# prevents funcitons from spamming the the console with unrequestted output
silent <- function(expr) {
  null <- capture.output(result <- expr)
  return(result)
}

# when unit test plot functions, we just check for the presence of a file 
# and delete it immediately
expect_plot <- function(..., file = "Rplots.pdf", delete = TRUE) 
{
  pdf(file = file)
  x <- expect_silent(...)
  dev.off()
  
  exists <- file.exists("Rplots.pdf")
  y <- expect_true(exists)
  if (delete & exists) file.remove(file)
  
  return(invisible(x))
}

methods <- "ALL"
k <- 5
set.seed(3804)
ds <- sampleFuncy(obsNr = 100, timeNr = 20, reg = TRUE, k = k, sd = 0.4)

test_that("Plotting regular data sets", {
  expect_plot(plotFuncy(ds@data))
})


test_that("All Cluster algorithms work", {
  
  res1 <- silent(funcit(
    methods = methods, data = ds@data, k = k,
    clusters = ds@clusters, seed = 2405,
    save.data = TRUE, parallel = parallel
  ))
  
  expect_plot(plot(res1))
  
  expect_plot(plot(res1, select = 1, type = "all"))
  expect_plot(plot(res1, select = c(1, 2), type = "centers"))
  expect_plot(plot(res1, select = 1, type = "dist2centers"))
  
  ## special plot functions for fitclust-object
  expect_plot(plot(res1, select = "fitfclust", type = "discrim"))
  expect_plot(plot(res1, select = "fitfclust", type = "conf"))
  
  ## special plot functions for FSCM-object
  expect_plot(plot(res1, select = "fscm", type = "overview"))
  expect_plot(plot(res1, select = "fscm", type = "deviations"))
  expect_plot(plot(res1, select = "fscm", type = "locations"))
})




set.seed(3805)
ds <- sampleFuncy(reg = FALSE, timeNrMin = 5, timeNrMax = 10, k = k, sd = 0.3)

test_that("Plotting irregular data set ", {
  expect_plot(plotFuncy(ds))
})

test_that("Making an irregular data set regular", {
  data <- expect_silent(regFuncy(ds@data, timeNr = 10, nbasis = 5,
                                 method = "interpolate"))
  
  
  data <- expect_silent(regFuncy(ds@data, baseType = "splines", timeNr = 10, 
                                 nbasis = 10, method = "project")$data)
  
  # to save time, maybe don't cluster again
  # Test all cluster methods again",
  res2 <- silent(funcit(methods = methods, seed = 2506, data = data, k = 4))    
  cl <- lapply(res2@models, class)
  
  expected <- structure(
    c("funcyOutMbc-fitfclust", "funcyOut", "funcyOut-iterSubspace", 
      "funcyOutMbc", "funcyOutMbc", "funcyOutMbc-fscm", "funcyOutMbc"
    ), .Names = c("fitfclust", "distclust", "iterSubspace", "funclust", 
                  "funHDDC", "fscm", "waveclust"))
  
  expect_equal(object = unlist(cl), expected = expected)
})


test_that("Test control arguments", {
  k <- 4
  set.seed(3806)
  ds <- sampleFuncy(timeNrMin = 5, timeNrMax = 10, reg = FALSE, k = k, sd = 0.3)
  
  a <- list(coeffsCalc = "estimate", average = TRUE)
  fpcCtrl <- as(a, "fpcCtrl")
  
  b <- list(maxit = 5, baseType = "eigenbasis", flexDim = TRUE)
  funcyCtrlMbc <- as(b, "funcyCtrlMbc")
  
  res3 <- silent(funcit(
    methods = c(1, 2, 3), data = ds@data, k = k,
    clusters = ds@clusters,
    fpcCtrl = fpcCtrl,
    funcyCtrl = funcyCtrlMbc,
    save.data = TRUE,
    parallel = parallel
  ))
  
  expect_plot(plot(res3, type = "all"))
  expect_plot(plot(res3, type = "accordance"))
  
  library(scatterplot3d)
  expect_plot(plot(res3, select = 3, type = "fpc"))
})

