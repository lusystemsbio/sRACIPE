getInteractions <- function(value) {
    interaction = 1
    nInteractions <- sum(value > 0)
    circuitInteractions <-
        data.frame(matrix(nrow = nInteractions, ncol = 3))
    geneNames <- rownames(value)
    for (i in seq_len(nrow(value))) {
        for (j in seq_len(nrow(value))) {
            if (as.integer(value[i, j]) > 0) {
                circuitInteractions[interaction, 2] <- geneNames[i]
                circuitInteractions[interaction, 1] <- geneNames[j]
                circuitInteractions[interaction, 3] <-
                    as.integer(value[i, j])
                interaction <- interaction + 1
            }
        }
    }
    names(circuitInteractions) <- c("Source", "Target", "Type")
    return(circuitInteractions)
}

#' @export
#' @title Generate parameter names for a circuit
#' @param circuit RacipeSE object or topology as data.frame or filename
#' @examples
#' rSet <- RacipeSE()
#' data("demoCircuit")
#' sracipeCircuit(rSet) <- demoCircuit
#' paramNames <- sRACIPE::sracipeGenParamNames(rSet)
#'
#' @return list
sracipeGenParamNames <- function(circuit = "inputs/test.tpo") {
    if (methods::is(circuit, "RacipeSE")) {
        rSet <- circuit
    }  else {
        rSet <- RacipeSE()
        sracipeCircuit(rSet) <- circuit
    }
    geneInteraction <- as.matrix(rowData(rSet))
    geneNames <- names(rSet)

    paramList <- list()
    tmp <- lapply(geneNames,  function(x)
        paste("G_", x, sep = ""))
    paramList <- append(paramList, tmp)
    tmp <- lapply(geneNames,  function(x)
        paste("K_", x, sep = ""))
    paramList <- append(paramList, tmp)
    tmp <- list()
    tmp2 <- list()
    tmp3 <- list()
    for (gene1 in seq_along(geneNames))
    {
        for (gene2 in seq_along(geneNames))
        {
            if (geneInteraction[gene1, gene2] > 0) {
                tmp <- append(tmp,
                              paste("TH_", geneNames[[gene2]], "_",
                                    geneNames[[gene1]], sep = ""))
                tmp2 <- append(tmp2,
                               paste("N_", geneNames[[gene2]], "_",
                                     geneNames[[gene1]], sep = ""))
                tmp3 <- append(tmp3,
                               paste("FC_", geneNames[[gene2]], "_",
                                     geneNames[[gene1]], sep = ""))
            }
        }
    }
    paramList <- do.call(c, list(paramList, tmp, tmp2, tmp3))

    return(as.character(paramList))
}
