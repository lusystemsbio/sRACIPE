
data("demoCircuit")
racipe <- RacipeSE()
sracipeCircuit(racipe) <- demoCircuit
expect_equal(annotation(racipe),"demoCircuit")

tmp1 <- sracipeCircuit(racipe)
storage.mode(demoCircuit$Type) <- "character"
sracipeCircuit(racipe) <- demoCircuit
tmp2 <- sracipeCircuit(racipe)
expect_equal(tmp1,tmp2)
