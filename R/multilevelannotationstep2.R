# Helper: Prepare environment and adduct/isotope settings
.prepare_adduct_isotope_context <- function(outloc1) {
    # Ensure output directory exists
    if (!dir.exists(outloc1)) {
        suppressWarnings(dir.create(outloc1, recursive = TRUE))
    }
    setwd(outloc1)
    message("[Stage 2] Using output directory: ", normalizePath(outloc1, winslash = "/", mustWork = FALSE))

    # Validate required inputs and load
    if (!file.exists("step1_results.Rda")) {
        stop(paste("step1_results.Rda not found in:", normalizePath(outloc1, winslash = "/", mustWork = FALSE)))
    }
    if (!file.exists("global_cor.Rda")) {
        stop(paste("global_cor.Rda not found in:", normalizePath(outloc1, winslash = "/", mustWork = FALSE)))
    }
    tryCatch(load("step1_results.Rda"), error = function(e) stop(paste("Failed to load step1_results.Rda:", e$message)))
    tryCatch(load("global_cor.Rda"), error = function(e) stop(paste("Failed to load global_cor.Rda:", e$message)))
    unlink("allmatches_with_isotopes.txt")

    # Ensure adduct weights exist with same default behavior
    if (is.na(adduct_weights)[1] == TRUE) {
        data(adduct_weights)
        adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
        adduct_weights1[1, ] <- c("M+H", 1)
        adduct_weights1[2, ] <- c("M-H", 1)
        adduct_weights <<- as.data.frame(adduct_weights1)
        colnames(adduct_weights) <<- c("Adduct", "Weight")
    }

    # Standardize mchemdata types/columns used downstream
    cnames <- colnames(mchemdata)
    cnames[2] <- "time"
    colnames(mchemdata) <<- as.character(cnames)
    mchemdata$mz <<- as.numeric(as.character(mchemdata$mz))
    mchemdata$time <<- as.numeric(as.character(mchemdata$time))

    # Respect externally set max.rt.diff
    if (is.na(max.rt.diff) == FALSE) {
        max_diff_rt <<- max.rt.diff
    }
}

# Helper: Build isotope search space per feature for a single chemical ID
.build_isotope_search_space <- function(curmchemdata, isop_res_md, mass_defect_window) {
    # Normalize RT cluster labels
    curmchemdata$Module_RTclust <- gsub(curmchemdata$Module_RTclust,
                                        pattern = "_[0-9]*", replacement = "")
    isop_res_md$Module_RTclust <- gsub(isop_res_md$Module_RTclust,
                                       pattern = "_[0-9]*", replacement = "")

    isp_masses_mz_data <- lapply(1:length(curmchemdata$mz), function(m) {
        isotope_group <- as.character(curmchemdata$ISgroup[m])
        module_rt_group <- as.character(curmchemdata$Module_RTclust[m])
        module_rt_group <- gsub(module_rt_group, pattern = "_[0-9]*", replacement = "")
        query_md <- curmchemdata$mz[m] - round(curmchemdata$mz[m])
        put_isp_masses_curmz_data <- isop_res_md[which(abs((isop_res_md$MD) - (query_md)) <
                                                         mass_defect_window &
                                                         isop_res_md$Module_RTclust == module_rt_group), ]
        put_isp_masses_curmz_data <- as.data.frame(put_isp_masses_curmz_data)
        return(put_isp_masses_curmz_data)
    })

    isp_masses_mz_data <- ldply(isp_masses_mz_data, rbind)
    cnames <- colnames(isp_masses_mz_data)
    cnames[5] <- "AvgIntensity"
    colnames(isp_masses_mz_data) <- cnames

    # Ensure numeric types used downstream
    isp_masses_mz_data$mz <- as.numeric(as.character(isp_masses_mz_data$mz))
    isp_masses_mz_data$time <- as.numeric(as.character(isp_masses_mz_data$time))

    return(as.data.frame(isp_masses_mz_data))
}

# Helper: Annotate a single chemical ID and compute intermediate score
.annotate_single_chemid <- function(chemid,
                                    mchemdata,
                                    isop_res_md,
                                    mass_defect_window,
                                    mass_defect_mode,
                                    corthresh,
                                    global_cor,
                                    mzid,
                                    max_diff_rt,
                                    adduct_table,
                                    adduct_weights,
                                    filter.by,
                                    max_isp,
                                    MplusH.abundance.ratio.check,
                                    outloc) {
    curmchemdata <- mchemdata[which(mchemdata$chemical_ID == chemid), ]
    curmchemdata$mz <- as.numeric(as.character(curmchemdata$mz))
    curmchemdata$time <- as.numeric(as.character(curmchemdata$time))
    curmchemdata <- as.data.frame(curmchemdata)

    isp_masses_mz_data <- .build_isotope_search_space(
        curmchemdata = curmchemdata,
        isop_res_md = isop_res_md,
        mass_defect_window = mass_defect_window
    )

    if (is.na(mass_defect_mode) == TRUE) {
        mass_defect_mode <- "pos"
    }

    chem_score <- get_chemscorev1.6.71(
        chemicalid = chemid,
        mchemicaldata = curmchemdata,
        corthresh = corthresh,
        global_cor = global_cor,
        mzid,
        max_diff_rt = max_diff_rt,
        level_module_isop_annot = isp_masses_mz_data,
        adduct_table = adduct_table,
        adduct_weights = adduct_weights,
        filter.by = filter.by, # preserve empty vector vs NULL behavior
        max_isp = max_isp,
        MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
        mass_defect_window = mass_defect_window,
        mass_defect_mode = mass_defect_mode,
        outlocorig = outloc
    )

    if (length(chem_score) > 0 && !is.null(chem_score)) {
        if (chem_score$chemical_score >= (-100)) {
            chem_score$filtdata <- chem_score$filtdata[order(chem_score$filtdata$mz), ]
            if (nrow(chem_score$filtdata) > 0) {
                cur_chem_score <- chem_score$chemical_score
                chemscoremat <- cbind(cur_chem_score, chem_score$filtdata)
                chemscoremat <- na.omit(chemscoremat)
                chemscoremat <- as.data.frame(chemscoremat)
                return(chemscoremat)
            }
        }
    }

    return(NULL)
}

# Helper: Pairwise annotation across a list of chem IDs (optionally parallel)
.run_pairwise_annotation <- function(chemid_indices,
                                     chemids,
                                     parallel = FALSE,
                                     num_nodes = 1,
                                     outloc) {
    local_cluster <- NULL

    if (parallel) {
        local_cluster <- makePSOCKcluster(num_nodes)
        on.exit(try(stopCluster(local_cluster), silent = TRUE), add = TRUE)
        # Load required packages on workers
        clusterEvalQ(local_cluster, { library(Rdisop); library(plyr) })
        # Export necessary symbols to workers
        clusterExport(local_cluster, varlist = c(
            "chemids", "mchemdata", "isop_res_md", "mass_defect_window",
            "mass_defect_mode", "corthresh", "global_cor", "mzid",
            "max_diff_rt", "adduct_table", "adduct_weights", "filter.by",
            "max_isp", "MplusH.abundance.ratio.check", "outloc",
            ".build_isotope_search_space", ".annotate_single_chemid"
        ), envir = environment())

        # Execute in parallel
        result_list <- parLapply(local_cluster, chemid_indices, function(j) {
            tryCatch({
                chemid <- chemids[j]
                .annotate_single_chemid(
                    chemid = chemid,
                    mchemdata = mchemdata,
                    isop_res_md = isop_res_md,
                    mass_defect_window = mass_defect_window,
                    mass_defect_mode = mass_defect_mode,
                    corthresh = corthresh,
                    global_cor = global_cor,
                    mzid = mzid,
                    max_diff_rt = max_diff_rt,
                    adduct_table = adduct_table,
                    adduct_weights = adduct_weights,
                    filter.by = filter.by,
                    max_isp = max_isp,
                    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
                    outloc = outloc
                )
            }, error = function(e) {
                warning(paste("[Stage 2] Failed to annotate chemid index", j, "(", chemids[j], "):", conditionMessage(e)))
                NULL
            })
        })

        # Ensure the cluster is stopped when parallel mode was used
        try(stopCluster(local_cluster), silent = TRUE)
    } else {
        # Sequential fallback
        result_list <- lapply(chemid_indices, function(j) {
            tryCatch({
                chemid <- chemids[j]
                .annotate_single_chemid(
                    chemid = chemid,
                    mchemdata = mchemdata,
                    isop_res_md = isop_res_md,
                    mass_defect_window = mass_defect_window,
                    mass_defect_mode = mass_defect_mode,
                    corthresh = corthresh,
                    global_cor = global_cor,
                    mzid = mzid,
                    max_diff_rt = max_diff_rt,
                    adduct_table = adduct_table,
                    adduct_weights = adduct_weights,
                    filter.by = filter.by,
                    max_isp = max_isp,
                    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
                    outloc = outloc
                )
            }, error = function(e) {
                warning(paste("[Stage 2] Failed to annotate chemid index", j, "(", chemids[j], "):", conditionMessage(e)))
                NULL
            })
        })
    }

    return(result_list)
}

# Helper: Aggregate intermediate scoring results and persist
.aggregate_and_persist_scores <- function(chem_score_list, list_number) {
    chem_score2 <- chem_score_list[which(chem_score_list != "NULL")]
    curchemscoremat <- ldply(chem_score2, rbind)

    # Stop if no results were produced for this list
    if (is.null(curchemscoremat) || length(curchemscoremat) == 0 || nrow(curchemscoremat) == 0) {
        stop("Stage2_results is NULL")
    }

    cur_fname <- paste("chem_score", list_number, ".Rda", sep = "")
    save(curchemscoremat, file = cur_fname)

    return(curchemscoremat)
}

multilevelannotationstep2 <- function(outloc1, list_number,
                                      parallel = FALSE,
                                      num_nodes = 1) {
    # Prepare context and data (adduct/isotope checks)
    .prepare_adduct_isotope_context(outloc1)

    outloc <- outloc1

    # Validate list_number bounds with original behavior
    if (list_number > length(chemids_split)) {
        list_number <- length(chemids_split)
        return(0)
    }

    if (list_number > num_sets) {
        list_number <- length(chemids_split)
        return(0)
    }

    # Stage 2 output directory
    outloc_stage2 <- paste(outloc, "/stage2/", sep = "")
    suppressWarnings(dir.create(outloc_stage2))
    setwd(outloc_stage2)
    message("[Stage 2] Processing list ", list_number, " of ", length(chemids_split),
            "; parallel=", parallel, ", workers=", num_nodes)

    # Pairwise annotation
    chem_score <- .run_pairwise_annotation(
        chemid_indices = chemids_split[[list_number]],
        chemids = chemids,
        parallel = parallel,
        num_nodes = num_nodes,
        outloc = outloc
    )

    # Aggregate results and persist to file (intermediate scoring)
    curchemscoremat <- .aggregate_and_persist_scores(
        chem_score_list = chem_score,
        list_number = list_number
    )
    message("[Stage 2] Saved intermediate scores for list ", list_number)

    # Selective cleanup to preserve original memory management
    rm(chem_score)
    gc()

    # Mimic original cleanup of heavy globals
    rm("mchemdata", "chemids", "adduct_table", "global_cor", "mzid",
       "max_diff_rt", "isop_res_md", "corthresh", "level_module_isop_annot",
       "chemids_split", "corthresh", "max.mz.diff", "outloc", "num_sets",
       "db_name", "num_nodes", "adduct_weights", "filter.by")
    suppressWarnings(rm(hmdbCompMZ))
    suppressWarnings(rm(hmdbAllinf))
    gc()

    # Return updated chemscoremat_highconf (use curchemscoremat here)
    chemscoremat_highconf <- curchemscoremat
    return(invisible(chemscoremat_highconf))
}
