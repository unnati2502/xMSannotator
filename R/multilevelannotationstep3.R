# Helper: Load prerequisites and set working directory
.load_step3_context <- function(outloc1) {
    if (!dir.exists(outloc1)) {
        suppressWarnings(dir.create(outloc1, recursive = TRUE))
    }
    setwd(outloc1)
    if (!file.exists("step1_results.Rda")) {
        stop("step1_results.Rda not found in outloc1: ", normalizePath(outloc1, winslash = "/", mustWork = FALSE))
    }
    if (!file.exists("chemCompMZ.Rda")) {
        stop("chemCompMZ.Rda not found in outloc1: ", normalizePath(outloc1, winslash = "/", mustWork = FALSE))
    }
    load("step1_results.Rda")
    load("chemCompMZ.Rda")
    # Cleanup globals as in original
    rm(global_cor)
    suppressWarnings(try(rm(global_cor, env = .GlobalEnv), silent = TRUE))
    gc()
}

# Helper: Read all Stage2 chem_score*.Rda and construct base graph table
.build_graph_from_stage2 <- function(outloc1, adduct_weights, num_sets) {
    # adduct weights default handling (preserve behavior)
    if (is.na(adduct_weights)[1] == TRUE) {
        data(adduct_weights)
        adduct_weights <- as.data.frame(adduct_weights)
        adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
        adduct_weights1[1, ] <- c("M+H", 1)
        adduct_weights1[2, ] <- c("M-H", 1)
        adduct_weights <- as.data.frame(adduct_weights1)
        colnames(adduct_weights) <- c("Adduct", "Weight")
    }

    # Determine number of sets
    num_sets_total <- length(chemids_split)
    if (!is.na(num_sets)) {
        num_sets_total <- num_sets
    }
    if (num_sets_total >= length(chemids_split)) {
        num_sets_total <- length(chemids_split)
    }

    # Read all chem_score*.Rda from Stage2
    stage2_dir <- file.path(outloc1, "stage2")
    suppressWarnings(dir.create(stage2_dir, recursive = TRUE))
    setwd(stage2_dir)

    message("[Stage 3] Reading Stage 2 score files from ", normalizePath(stage2_dir, winslash = "/", mustWork = FALSE))
    chemscore_list <- lapply(1:num_sets_total, function(sind) {
        cur_fname <- paste("chem_score", sind, ".Rda", sep = "")
        ok <- try(load(cur_fname), silent = TRUE)
        if (inherits(ok, "try-error") || !exists("curchemscoremat")) {
            return(NULL)
        }
        cur <- as.data.frame(curchemscoremat)
        # Normalize formula/adduct fields for graph-building
        cur$Formula <- gsub(cur$Formula, pattern = "_.*", replacement = "")
        return(cur)
    })

    chemscore_list <- chemscore_list[!vapply(chemscore_list, is.null, logical(1))]
    if (length(chemscore_list) == 0) {
        stop("No Stage 2 chem_score*.Rda files could be read. Stage 3 cannot proceed.")
    }

    chemscoremat <- ldply(chemscore_list, rbind)
    chemscoremat <- as.data.frame(chemscoremat)

    # Preserve original adduct (with isotope suffix) in temp column
    tempadduct <- chemscoremat$Adduct
    chemscoremat$Adduct <- gsub(chemscoremat$Adduct, pattern = "_.*", replacement = "")
    chemscoremat <- cbind(chemscoremat, tempadduct)

    # Merge with chemCompMZ (subset of columns 2:4, 6)
    library(data.table)
    setDT(chemscoremat)
    setDT(chemCompMZ)
    chemCompMZ_subset <- chemCompMZ[, .SD, .SDcols = c(2:4, 6)]
    chemscoremat <- merge(
        chemscoremat,
        chemCompMZ_subset,
        by = c("Formula", "Adduct"),
        allow.cartesian = TRUE
    )
    chemscoremat <- as.data.frame(chemscoremat)

    # Select and rename columns per original layout
    chemscoremat <- chemscoremat[, c("cur_chem_score", "Module_RTclust",
                                     "mz", "time", "MatchCategory", "theoretical.mz",
                                     "chemical_ID.y", "Name.y", "Formula", "MonoisotopicMass",
                                     "tempadduct", "ISgroup", "mean_int_vec", "MD")]
    colnames(chemscoremat) <- c("cur_chem_score", "Module_RTclust",
                                "mz", "time", "MatchCategory", "theoretical.mz",
                                "chemical_ID", "Name", "Formula", "MonoisotopicMass",
                                "Adduct", "ISgroup", "mean_int_vec", "MD")

    # Remove known bad HMDB IDs (as in original)
    hmdbbad <- c("HMDB29244", "HMDB29245", "HMDB29246")
    if (length(which(chemscoremat$chemical_ID %in% hmdbbad)) > 0) {
        chemscoremat <- chemscoremat[-which(chemscoremat$chemical_ID %in% hmdbbad), ]
    }

    # Standardize column names
    cnames <- colnames(chemscoremat)
    cnames[1] <- "score"
    colnames(chemscoremat) <- cnames

    # Remove rows with missing values in critical fields (restore original behavior and user-facing message)
    critical_cols <- c("Formula", "Adduct", "mz", "time", "chemical_ID", "theoretical.mz")
    critical_cols <- critical_cols[critical_cols %in% colnames(chemscoremat)]
    if (length(critical_cols) > 0) {
        complete_idx <- stats::complete.cases(chemscoremat[, critical_cols, drop = FALSE])
        if (any(!complete_idx)) {
            cat("Removing rows with missing values\n")
            chemscoremat <- chemscoremat[complete_idx, , drop = FALSE]
        }
    }

    return(chemscoremat)
}

# Helper: Apply pathway-based network scoring adjustments (DB-specific)
.apply_pathway_network_scoring <- function(chemscoremat, db_name, adduct_weights, scorethresh, pathwaycheckmode, max_diff_rt) {
    pthresh <- 0.05

    if (db_name == "KEGG") {
        data(keggotherinf)
        m1 <- apply(keggotherinf, 1, function(x) {
            chemid <- x[1]
            g1 <- gregexpr(x[4], pattern = "map")
            regexp_check <- attr(g1[[1]], "match.length")
            if (regexp_check[1] < 0) {
                pathid = "-"
                return(cbind(chemid, pathid))
            } else {
                pathid <- strsplit(x = x[4], split = ";")
                pathid <- unlist(pathid)
                return(cbind(chemid, pathid))
            }
        })
        m2 <- ldply(m1, rbind)

        # Use data.table merge to append KEGG info
        setDT(chemscoremat)
        setDT(keggotherinf)
        chemscoremat <- merge(
          chemscoremat,
          keggotherinf,
          by.x = "chemical_ID",
          by.y = "KEGGID",
          all = FALSE,
          allow.cartesian = TRUE
        )
        chemscoremat <- as.data.frame(chemscoremat)

        chemids <- as.character(chemscoremat$chemical_ID)
        pathway_ids <- unique(as.character(m2[, 2]))
        module_num <- gsub(chemscoremat$Module_RTclust, pattern = "_[0-9]*", replacement = "")
        chemscoremat <- cbind(chemscoremat, module_num)
        chemscoremat_orig <- chemscoremat
        chemscoremat <- chemscoremat_orig
        total_chem_count <- length(unique(m2$chemid))

        if (is.na(pathwaycheckmode) == FALSE) {
            for (path_id in pathway_ids) {
                if (path_id != "-" & path_id != "map01100") {
                    pathway_chemicals <- m2[which(m2[, 2] %in% path_id), 1]
                    curmchemicaldata1 <- chemscoremat[which(chemscoremat$chemical_ID %in% pathway_chemicals &
                                                             chemscoremat$score >= scorethresh &
                                                             chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]
                    num_chems_inpath <- length(unique(curmchemicaldata1$chemical_ID))
                  all_cur_path_numchem <- length(unique(pathway_chemicals))
                    curmchemicaldata2 <- chemscoremat[which(chemscoremat$score >= scorethresh &
                                                             chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]
                    curmchemicaldata2 <- curmchemicaldata2[-which(curmchemicaldata2$chemical_ID %in% pathway_chemicals), ]
                  num_chems_notinpath <- length(unique(curmchemicaldata2$chemical_ID))
                    all_notcurpath_numchem <- length(m2[-which(m2[, 1] %in% curmchemicaldata2$chemical_ID), 1])
                    counts <- matrix(data = c(num_chems_inpath,
                                              all_cur_path_numchem - num_chems_inpath,
                                              num_chems_notinpath,
                    all_notcurpath_numchem), nrow = 2)
                    rm(curmchemicaldata2); gc()
                    p1 <- fisher.test(counts)$p.value
                    if (p1 <= pthresh) {
                    t1 <- table(curmchemicaldata1$module_num)
                    module_counts <- t1[which(t1 > 0)]
                    module_names <- names(module_counts)
                    pathway_chemicals_1 <- curmchemicaldata1$chemical_ID
                    for (c in pathway_chemicals_1) {
                            curmchemicaldata <- curmchemicaldata1[which(as.character(curmchemicaldata1$chemical_ID) == c), ]
                            chem_path_data <- m2[which(m2$chemid == c), ]
                      total_num_chem <- 0
                      t2 <- table(curmchemicaldata$module_num)
                            cur_module <- names(t2[which(t2 == max(t2)[1])])
                      total_num_chem <- length(pathway_chemicals)
                            mzid_cur <- paste(curmchemicaldata$mz, curmchemicaldata$time, sep = "_")
                      if (nrow(curmchemicaldata) > 0) {
                                # per-module counts and Fisher test
                            num_chems <- t1[as.character(cur_module)]
                          num_chems <- round(num_chems, 0)
                                a <- length(unique(chemscoremat[which(chemscoremat$module_num == cur_module &
                                                                       chemscoremat$chemical_ID %in% pathway_chemicals &
                                                                       chemscoremat$score >= scorethresh &
                                                                       chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]$chemical_ID))
                                b <- length(unique(chemscoremat[which(chemscoremat$module_num == cur_module &
                                                                       chemscoremat$score >= scorethresh &
                                                                       chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]$chemical_ID)) - a
                                c_count <- length(unique(chemscoremat[which(chemscoremat$module_num != cur_module &
                                                                             chemscoremat$chemical_ID %in% pathway_chemicals &
                                                                             chemscoremat$score >= scorethresh &
                                                                             chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]$chemical_ID))
                                d <- length(unique(chemscoremat[which(chemscoremat$module_num != cur_module &
                                                                      chemscoremat$score >= scorethresh &
                                                                      chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]$chemical_ID)) - c_count
                                counts_md <- matrix(data = c(a, c_count, b, d), nrow = 2)
                                p1 <- if (a > 1) fisher.test(counts_md)$p.value else 1
                                if (p1 <= 0.2) {
                            if (num_chems[1] < 3) {
                              num_chems <- 0
                            } else {
                                        if (is.na(curmchemicaldata$score[1]) == TRUE) {
                                            diff_rt <- max(curmchemicaldata$time) - min(curmchemicaldata$time)
                                if (diff_rt > max_diff_rt) {
                                                if (length(which(t2 > 1)) == 1) {
                                                    curmchemicaldata$score <- rep(0.1, length(curmchemicaldata$score))
                                  } else {
                                                    curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                                  }
                                } else {
                                                curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                                }
                              }
                                        if (curmchemicaldata$score[1] < scorethresh) {
                                            chemscoremat$score[which(as.character(chemscoremat$chemical_ID) == c &
                                                                     chemscoremat$Adduct %in% as.character(adduct_weights[, 1]))] <-
                                                as.numeric(chemscoremat$score[which(as.character(chemscoremat$chemical_ID) == c &
                                                                                    chemscoremat$Adduct %in% as.character(adduct_weights[, 1]))][1]) + num_chems
                              } else {
                                            chemscoremat$score[which(as.character(chemscoremat$chemical_ID) == c)] <-
                                                chemscoremat$score[which(as.character(chemscoremat$chemical_ID) == c)][1] + num_chems
                              }
                            }
                          }
                        }
                      }
                    }
                }
                rm(curmchemicaldata1); gc()
            }
            # Standardize column names after merge side-effects
        cnames <- colnames(chemscoremat)
        cnames <- gsub(cnames, pattern = ".x", replacement = "")
        colnames(chemscoremat) <- cnames
        }
    } else if (db_name == "HMDB") {
            data(hmdbAllinf)
        hmdbAllinfv3.5 <- hmdbAllinf[, -c(26:27)]
        suppressWarnings(try(rm(hmdbAllinf, envir = .GlobalEnv), silent = TRUE))
            gc()
            m1 <- apply(hmdbAllinfv3.5, 1, function(x) {
                chemid <- x[1]
                g1 <- gregexpr(x[14], pattern = "SM")
                regexp_check <- attr(g1[[1]], "match.length")
                if (regexp_check[1] < 0) {
                  pathid = "-"
                  return(cbind(chemid, pathid))
                } else {
                  pathid <- strsplit(x = x[14], split = ";")
                  pathid <- unlist(pathid)
                  return(cbind(chemid, pathid))
                }
            })
            m2 <- ldply(m1, rbind)
            setDT(chemscoremat)
            setDT(hmdbAllinfv3.5)
            chemscoremat <- merge(chemscoremat, hmdbAllinfv3.5,
                by.x = "chemical_ID", by.y = "HMDBID", allow.cartesian = TRUE)
            chemscoremat <- as.data.frame(chemscoremat)

            chemids <- as.character(chemscoremat$chemical_ID)
        pathway_ids <- unique(as.character(m2[, 2]))
        module_num <- gsub(chemscoremat$Module_RTclust, pattern = "_[0-9]*", replacement = "")
            chemscoremat <- cbind(chemscoremat, module_num)
            chemscoremat_orig <- chemscoremat
            chemscoremat <- chemscoremat_orig
            total_chem_count <- length(unique(m2$chemid))

            if (is.na(pathwaycheckmode) == FALSE) {
                for (path_id in pathway_ids) {
                  if (path_id != "-") {
                    pathway_chemicals <- m2[which(m2[, 2] %in% path_id), 1]
                    curmchemicaldata1 <- chemscoremat[which(chemscoremat$chemical_ID %in% pathway_chemicals &
                                                             chemscoremat$score >= scorethresh &
                                                             chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]
                    num_chems_inpath <- length(unique(curmchemicaldata1$chemical_ID))
                    all_cur_path_numchem <- length(unique(pathway_chemicals))
                    curmchemicaldata2 <- chemscoremat[which(chemscoremat$score >= scorethresh &
                                                             chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]
                    curmchemicaldata2 <- curmchemicaldata2[-which(curmchemicaldata2$chemical_ID %in% pathway_chemicals), ]
                    num_chems_notinpath <- length(unique(curmchemicaldata2$chemical_ID))
                    all_notcurpath_numchem <- length(m2[-which(m2[, 1] %in% curmchemicaldata2$chemical_ID), 1])
                    counts <- matrix(data = c(num_chems_inpath,
                                              all_cur_path_numchem - num_chems_inpath,
                                              num_chems_notinpath,
                      all_notcurpath_numchem), nrow = 2)
                    rm(curmchemicaldata2); gc()
                    p1 <- fisher.test(counts)$p.value
                    if (p1 <= pthresh) {
                      t1 <- table(curmchemicaldata1$module_num)
                      module_counts <- t1[which(t1 > 0)]
                      module_names <- names(module_counts)
                      pathway_chemicals_1 <- curmchemicaldata1$chemical_ID
                      for (chemname in pathway_chemicals) {
                            curmchemicaldata <- curmchemicaldata1[which(as.character(curmchemicaldata1$chemical_ID) == chemname), ]
                        if (nrow(curmchemicaldata) > 0) {
                          t2 <- table(curmchemicaldata$module_num)
                                cur_module <- names(t2[which(t2 == max(t2)[1])])
                                mzid_cur <- paste(curmchemicaldata$mz, curmchemicaldata$time, sep = "_")
                              num_chems <- t1[as.character(cur_module)]
                            total_num_chem <- length(pathway_chemicals)
                                num_chems <- round(num_chems, 0)
                                a <- length(unique(chemscoremat[which(chemscoremat$module_num == cur_module &
                                                                       chemscoremat$chemical_ID %in% pathway_chemicals &
                                                                       chemscoremat$score >= scorethresh &
                                                                       chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]$chemical_ID))
                                b <- length(unique(chemscoremat[which(chemscoremat$module_num == cur_module &
                                                                       chemscoremat$score >= scorethresh &
                                                                       chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]$chemical_ID)) - a
                                c_count <- length(unique(chemscoremat[which(chemscoremat$module_num != cur_module &
                                                                             chemscoremat$chemical_ID %in% pathway_chemicals &
                                                                             chemscoremat$score >= scorethresh &
                                                                             chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]$chemical_ID))
                                d <- length(unique(chemscoremat[which(chemscoremat$module_num != cur_module &
                                                                      chemscoremat$score >= scorethresh &
                                                                      chemscoremat$Adduct %in% as.character(adduct_weights[, 1])), ]$chemical_ID)) - c_count
                                counts_md <- matrix(data = c(a, c_count, b, d), nrow = 2)
                                p1 <- if (a > 1) fisher.test(counts_md)$p.value else 1
                                if (p1 <= 0.2) {
                              if (num_chems[1] < 3) {
                                num_chems <- 0
                                    } else {
                                        if (is.na(curmchemicaldata$score[1]) == TRUE) {
                                            diff_rt <- max(curmchemicaldata$time) - min(curmchemicaldata$time)
                                            if (diff_rt > max_diff_rt) {
                                                if (length(which(t2 > 1)) == 1) {
                                                    curmchemicaldata$score <- rep(0.1, length(curmchemicaldata$score))
                                                } else {
                                                    curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                                    }
                                  } else {
                                                curmchemicaldata$score <- rep(0, length(curmchemicaldata$score))
                                  }
                                }
                                        if (curmchemicaldata$score[1] < scorethresh) {
                                            chemscoremat$score[which(as.character(chemscoremat$chemical_ID) == chemname &
                                                                     chemscoremat$Adduct %in% as.character(adduct_weights[, 1]))] <-
                                                as.numeric(chemscoremat$score[which(as.character(chemscoremat$chemical_ID) == chemname &
                                                                                    chemscoremat$Adduct %in% as.character(adduct_weights[, 1]))][1]) + num_chems
                                } else {
                                            chemscoremat$score[which(as.character(chemscoremat$chemical_ID) == chemname)] <-
                                                chemscoremat$score[which(as.character(chemscoremat$chemical_ID) == chemname)][1] + num_chems
                                }
                              }
                            }
                          }
                        }
                    }
                  }
                rm(curmchemicaldata1); gc()
                }
            cnames <- colnames(chemscoremat)
            cnames <- gsub(cnames, pattern = ".x", replacement = "")
            colnames(chemscoremat) <- cnames
        }
    }

    return(chemscoremat)
}

# Helper: Extract high-confidence subset and write Stage3.csv
.extract_and_write_stage3 <- function(outloc1, chemscoremat, scorethresh) {
    # Standardize and subset columns
    chemscoremat <- chemscoremat
    good_ind <- which(chemscoremat$score >= scorethresh)
    chemscoremat_highconf <- chemscoremat

    chemscoremat_highconf <- chemscoremat_highconf[, c("chemical_ID",
        "score", "Module_RTclust", "mz", "time", "MatchCategory",
        "theoretical.mz", "Name", "Formula", "MonoisotopicMass",
        "Adduct", "ISgroup", "mean_int_vec", "MD")]

    stage3_path <- file.path(outloc1, "Stage3.csv")
    write.csv(chemscoremat_highconf, file = stage3_path, row.names = FALSE)
    message("[Stage 3] Wrote Stage3.csv to ", normalizePath(stage3_path, winslash = "/", mustWork = FALSE))

    return(chemscoremat_highconf)
}

multilevelannotationstep3 <- function(outloc1, adduct_weights = NA,
    num_sets = NA, boostIDs = NA, pathwaycheckmode = "p",
    dbAllinf = NA, scorethresh = 0.1) {

    .load_step3_context(outloc1)

    # Derive num_sets consistent with original behavior
    num_sets_1 <- if (!is.na(num_sets)) num_sets else NA
    num_sets_final <- length(chemids_split)
    if (!is.na(num_sets_1)) num_sets_final <- num_sets_1
    if (num_sets_final >= length(chemids_split)) num_sets_final <- length(chemids_split)

    # Graph construction
    chemscoremat <- .build_graph_from_stage2(outloc1, adduct_weights, num_sets_final)

    # Network scoring
    chemscoremat <- .apply_pathway_network_scoring(
        chemscoremat = chemscoremat,
        db_name = db_name,
        adduct_weights = adduct_weights,
        scorethresh = scorethresh,
        pathwaycheckmode = pathwaycheckmode,
        max_diff_rt = max_diff_rt
    )

    # Subset extraction + write Stage3.csv
    chemscoremat_highconf <- .extract_and_write_stage3(outloc1, chemscoremat, scorethresh)

    # Cleanup mirroring original intent
    suppressWarnings(try(rm(chemCompMZ), silent = TRUE))
    suppressWarnings(try(rm(mchemdata), silent = TRUE))
    suppressWarnings(try(rm(hmdbAllinf), silent = TRUE))
    suppressWarnings(try(rm(hmdbAllinfv3.6), silent = TRUE))
    suppressWarnings(try(rm(dbAllinf), silent = TRUE))
    gc()

    return(chemscoremat_highconf)
}
