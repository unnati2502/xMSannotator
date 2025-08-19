# Helper: load and standardize Stage 3 input
.load_stage3_data <- function(outloc) {
    # Ensure output directory exists
    if (!dir.exists(outloc)) {
        suppressWarnings(dir.create(outloc, recursive = TRUE))
    }
    setwd(outloc)

    # Validate existence of Stage3.csv before attempting to read
    if (!file.exists("Stage3.csv")) {
        stop("Stage3.csv not found in outloc: ", normalizePath(outloc, winslash = "/", mustWork = FALSE))
    }

    chemscoremat_highconf <- read.csv("Stage3.csv")
    chemscoremat_highconf <- as.data.frame(chemscoremat_highconf)
    chemscoremat_highconf$mz <- as.numeric(chemscoremat_highconf$mz)

    cnames <- colnames(chemscoremat_highconf)
    cnames <- gsub(cnames, pattern = ".x", replacement = "")
    colnames(chemscoremat_highconf) <- cnames

    return(chemscoremat_highconf)
}

# Helper: get unique chemical IDs
.get_unique_chemids <- function(chemscoremat_highconf) {
    chemids <- unique(chemscoremat_highconf$chemical_ID)
    return(chemids)
}

# Helper: initialize adduct-related settings
.prepare_adduct_settings <- function(adduct_weights) {
    data(adduct_table)

    if (is.na(adduct_weights)[1] == TRUE) {
        data(adduct_weights)
        adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
        adduct_weights1[1,] <- c("M+H", 1)
        adduct_weights1[2,] <- c("M-H", 1)
        adduct_weights <- as.data.frame(adduct_weights1)
        colnames(adduct_weights) <- c("Adduct", "Weight")
    }

    adduct_table <- adduct_table[order(adduct_table$Adduct),]

    return(list(adduct_table = adduct_table, adduct_weights = adduct_weights))
}

# Helper: compute confidence levels per chemical ID (potentially parallel)
.compute_confidence_levels <- function(chemscoremat_highconf,
                                       chemids,
                                       max_diff_rt,
                                       adduct_weights,
                                       filter.by,
                                       max_isp,
                                       min_ions_perchem,
                                       adduct_table,
                                       num_nodes) {
    # winsize retained for compatibility with existing logic (not actively used)
    winsize <- 500

    # Initialize in case parallelization fails
    chemscoremat_conf_levels <- data.frame(Confidence = numeric(0), chemical_ID = character(0))

    {
        cl <- makeSOCKcluster(num_nodes)

        # Register foreach backend and guarantee cleanup/reset
        doSNOW::registerDoSNOW(cl)
        on.exit({
            try(suppressWarnings(stopCluster(cl)), silent = TRUE)
            try(foreach::registerDoSEQ(), silent = TRUE)
        }, add = TRUE)

        # Keep exports/eval calls to preserve original behavior and environment propagation
        clusterEvalQ(cl, "get_confidence_stage4")
        clusterEvalQ(cl, "check_element")
        clusterEvalQ(cl, "group_by_rt")
        clusterEvalQ(cl, "library(doSNOW)")

        clusterExport(cl, "multilevelannotationstep2")

        clusterEvalQ(cl, "library(Rdisop)")
        clusterEvalQ(cl, "library(plyr)")
        # clusterEvalQ(cl, 'library(pryr)')
        # clusterEvalQ(cl, 'library(profmem)')
        # clusterEvalQ(cl, 'library(gdata)')

        clusterExport(cl, "getMolecule", envir = environment())
        clusterExport(cl, "ldply", envir = environment())
        clusterExport(cl, "max_isp", envir = environment())
        clusterExport(cl, "get_confidence_stage4", envir = environment())

        clusterExport(cl, "min_ions_perchem", envir = environment())
        clusterExport(cl, "check_element", envir = environment())
        clusterExport(cl, "group_by_rt", envir = environment())
        clusterExport(cl, "adduct_table", envir = environment())
        clusterExport(cl, "adduct_weights", envir = environment())
        clusterExport(cl, "filter.by", envir = environment())
        clusterExport(cl, "max_diff_rt", envir = environment())
        clusterExport(cl, "chemscoremat_highconf", envir = environment())
        clusterExport(cl, "chemids", envir = environment())

        # NOTE: foreach backend registered above. Errors per task are handled to avoid full-stop.
        chemscoremat_conf_levels <- foreach(c = 1:length(chemids), .combine = rbind) %dopar% {
            tryCatch({
                cur_chemid <- chemids[c]
                curdata <- chemscoremat_highconf[chemscoremat_highconf$chemical_ID == cur_chemid, ]
                curdata <- curdata[order(curdata$Adduct), ]
                bool_check <- 1

                if (!is.na(filter.by)[1] && bool_check == 1) {
                    check_adduct <- which(curdata$Adduct %in% filter.by)
                    bool_check <- if (length(check_adduct) > 0) 1 else 0
                }

                if (bool_check == 1) {
                    final_res <- get_confidence_stage4(
                        curdata,
                        max_diff_rt,
                        adduct_weights = adduct_weights,
                        filter.by = filter.by,
                        max_isp = max_isp,
                        min_ions_perchem = min_ions_perchem
                    )

                    if (!identical(final_res, "None") && !is.na(final_res[1, 1])) {
                        if (nrow(final_res) == nrow(curdata)) {
                            Confidence <- as.numeric(as.character(final_res[, 1]))
                            curdata <- final_res
                        } else {
                            warning(paste("Row mismatch for", cur_chemid, "- skipping."))
                            return(data.frame(Confidence = 0, chemical_ID = cur_chemid))
                        }
                    } else {
                        Confidence <- 0
                    }
                } else {
                    Confidence <- 0
                }

                curdata <- cbind(Confidence, curdata)
                curdata <- as.data.frame(curdata)
                curdata <- curdata[, c("Confidence", "chemical_ID")]
                curdata <- unique(curdata)
                return(curdata)
            }, error = function(e) {
                warning(paste("Failed to process chemical ID:", chemids[c], "-", e$message))
                return(data.frame(Confidence = 0, chemical_ID = chemids[c]))
            })
        }
    }

    chemscoremat_conf_levels <- as.data.frame(chemscoremat_conf_levels)
    return(chemscoremat_conf_levels)
}

# Helper: merge confidence with Stage 3 data
.merge_confidence_with_stage3 <- function(chemscoremat_conf_levels, chemscoremat_highconf) {
    library(data.table)

    setDT(chemscoremat_conf_levels)
    setDT(chemscoremat_highconf)

    curated_res <- merge(
        chemscoremat_conf_levels,
        chemscoremat_highconf,
        by = "chemical_ID",
        all = FALSE
    )

    curated_res <- as.data.frame(curated_res)
    colnames(curated_res) <- as.character(colnames(curated_res))

    return(curated_res)
}

# Helper: compute and add delta ppm
.add_delta_ppm <- function(curated_res) {
    curated_res <- as.data.frame(curated_res)

    curated_res$mz <- as.numeric(as.character(curated_res$mz))
    curated_res$theoretical.mz <- as.numeric(as.character(curated_res$theoretical.mz))

    curated_res_temp <- curated_res[, c("mz", "theoretical.mz")]
    curated_res_temp <- apply(curated_res_temp, 1, as.numeric)
    curated_res_temp <- t(curated_res_temp)
    curated_res_temp <- as.data.frame(curated_res_temp)

    delta_ppm <- apply(curated_res_temp, 1, function(x) {
        ppmerror <- 10 ^ 6 * abs(x[2] - x[1]) / (x[2])
        return(ppmerror)
    })
    delta_ppm <- round(delta_ppm, 2)

    curated_res <- cbind(curated_res[, 1:8], delta_ppm, curated_res[, 9:dim(curated_res)[2]])
    curated_res <- curated_res[order(curated_res$Confidence, decreasing = TRUE),]

    return(curated_res)
}

# Helper: optionally boost IDs if provided
.apply_boost_ids <- function(curated_res, boostIDs, max.mz.diff, max_diff_rt) {
    if (is.na(boostIDs)[1] == TRUE) {
        return(curated_res)
    }

    if (length(colnames(boostIDs)) > 1) {
        curated_res_mzrt <- curated_res[, c("mz", "time")]
        validated_mzrt <- boostIDs[, c("mz", "time")]

        ghilicpos <- getVenn(
            curated_res_mzrt,
            name_a = "exp",
            validated_mzrt,
            name_b = "boost",
            mz.thresh = max.mz.diff,
            time.thresh = max_diff_rt,
            alignment.tool = NA,
            xMSanalyzer.outloc = getwd(),
            use.unique.mz = FALSE,
            plotvenn = FALSE
        )

        save(ghilicpos, file = "ghilicpos.Rda")

        g1 <- ghilicpos$common
        rm(ghilicpos)
        gc()

        t1 <- table(curated_res$Confidence, curated_res$chemical_ID)
        cnames <- colnames(t1)
        cnames <- cnames[which(cnames %in% boostIDs$ID)]

        good_ind_1 <- NULL

        for (ind2 in 1:dim(g1)[1]) {
            temp_ind1 <- g1$index.A[ind2]
            temp_ind2 <- g1$index.B[ind2]

            if ((curated_res$chemical_ID[temp_ind1] %in% boostIDs$ID[temp_ind2])[1]) {
                good_ind_1 <- c(good_ind_1, g1$index.A[ind2])
            }
        }

        overlap_mz_time_id <- good_ind_1

        curated_res$Confidence[overlap_mz_time_id] <- 4
        curated_res$score[overlap_mz_time_id] <- curated_res$score[overlap_mz_time_id] * 100
        t1 <- table(curated_res$Confidence[overlap_mz_time_id], curated_res$chemical_ID[overlap_mz_time_id])

        cnames1 <- colnames(t1)
        cnames2 <- cnames1[which(t1 > 0)]
        good_ind <- NULL
        if (length(good_ind) > 0) {
            curated_res$Confidence[good_ind] <- 4
            curated_res$score[good_ind] <- curated_res$score[good_ind] * 100
        }
    } else {
        good_ind <- which(curated_res$chemical_ID %in% boostIDs)
        if (length(good_ind) > 0) {
            curated_res$Confidence[good_ind] <- 4
            curated_res$score[good_ind] <- curated_res$score[good_ind] * 100
        }
    }

    return(curated_res)
}

# Helper: annotate match category (Unique vs Multiple)
.annotate_match_category <- function(curated_res) {
    t2 <- table(curated_res$mz)
    same1 <- which(t2 == 1)
    uniquemz <- names(same1)

    curated_res$MatchCategory <- rep("Multiple", dim(curated_res)[1])
    curated_res$MatchCategory[which(curated_res$mz %in% uniquemz)] <- "Unique"

    return(curated_res)
}

# Helper: write outputs
.write_stage4_outputs <- function(curated_res, outloc) {
    outloc3 <- outloc
    suppressWarnings(dir.create(outloc3))
    setwd(outloc3)

    write.csv(curated_res, file = "Stage4.csv", row.names = FALSE)
}

# Helper: print summaries
.print_stage4_summaries <- function(curated_res) {
    curated_res <- curated_res[order(curated_res$Confidence, decreasing = TRUE),]

    print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
    print(table(curated_res$Confidence[-which(duplicated(curated_res$chemical_ID) == TRUE)]))

    print("Stage 4 confidence level distribution for unique chemical/metabolite formulas")
    print(table(curated_res$Confidence[-which(duplicated(curated_res$Formula) == TRUE)]))
}

# Main function: refactored to orchestrate helpers
multilevelannotationstep4 <- function(outloc,
                                      max.mz.diff = 5,
                                      max.rt.diff = 30,
                                      adduct_weights = NA,
                                      filter.by = NA,
                                      min_ions_perchem = 1,
                                      boostIDs = NA,
                                      max_isp = 5,
                                      dbAllinf = NA,
                                      num_nodes = 2) {
    max_diff_rt <- max.rt.diff

    # Load and standardize Stage 3 data
    message("[Stage 4] Loading Stage3.csv from: ", normalizePath(outloc, winslash = "/", mustWork = FALSE))
    chemscoremat_highconf <- .load_stage3_data(outloc)

    # Discover unique chemicals
    chemids <- .get_unique_chemids(chemscoremat_highconf)
    message("[Stage 4] Found ", length(chemids), " unique chemical_IDs")

    # Initialize adduct configurations
    adduct_cfg <- .prepare_adduct_settings(adduct_weights)
    adduct_table <- adduct_cfg$adduct_table
    adduct_weights <- adduct_cfg$adduct_weights

    # Compute confidence levels (parallel-capable)
    message("[Stage 4] Computing confidence levels with up to ", num_nodes, " workers ...")
    chemscoremat_conf_levels <- .compute_confidence_levels(
        chemscoremat_highconf = chemscoremat_highconf,
        chemids = chemids,
        max_diff_rt = max_diff_rt,
        adduct_weights = adduct_weights,
        filter.by = filter.by,
        max_isp = max_isp,
        min_ions_perchem = min_ions_perchem,
        adduct_table = adduct_table,
        num_nodes = num_nodes
    )

    chemscoremat_highconf <- unique(chemscoremat_highconf)
    chemscoremat_conf_levels <- as.data.frame(chemscoremat_conf_levels)

    # Merge confidence with Stage 3
    curated_res <- .merge_confidence_with_stage3(
        chemscoremat_conf_levels = chemscoremat_conf_levels,
        chemscoremat_highconf = chemscoremat_highconf
    )

    rm(chemscoremat_highconf)
    gc()

    # Add delta ppm and sort by confidence
    curated_res <- .add_delta_ppm(curated_res)

    # Boost IDs if provided
    curated_res <- .apply_boost_ids(
        curated_res = curated_res,
        boostIDs = boostIDs,
        max.mz.diff = max.mz.diff,
        max_diff_rt = max_diff_rt
    )

    # Annotate unique/multiple matches
    curated_res <- .annotate_match_category(curated_res)

    # Persist outputs
    .write_stage4_outputs(curated_res, outloc)

    # Print summaries
    .print_stage4_summaries(curated_res)

    return(as.data.frame(curated_res))
}
