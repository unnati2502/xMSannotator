# Helper: Sanity checks for inputs
.check_sanity <- function(dataA, simmat, adjacencyfromsimilarity) {
	if (is.na(dataA)[1] == TRUE) {
		stop("dataA must be provided")
	}
	if (!is.data.frame(dataA) && !is.matrix(dataA)) {
		stop("dataA must be a data.frame or matrix")
	}
	if (ncol(dataA) < 2) {
		stop("dataA must have at least two columns: mz and time")
	}
	if (adjacencyfromsimilarity == TRUE) {
		if (is.na(simmat)[1] == TRUE) {
			stop("simmat must be provided when adjacencyfromsimilarity=TRUE")
		}
	}
	invisible(TRUE)
}

# Helper: Prepare expression matrix and metadata (mz, time)
.prepare_feature_matrix <- function(dataA, column.rm.index, step1log2scale) {
	cnames <- colnames(dataA)
	cnames[1] <- "mz"
	cnames[2] <- "time"
	colnames(dataA) <- as.character(cnames)

	data_mzrt <- dataA[, c(1:2)]

	if (is.na(column.rm.index) == FALSE) {
		dataA <- dataA[, -c(column.rm.index)]
	}

	feat_inf <- paste(dataA[, 1], dataA[, 2], sep = "_")
	dataA_vals <- dataA[, -c(1:2)]
	data_m <- t(dataA_vals)
	colnames(data_m) <- feat_inf

	if (step1log2scale == TRUE) {
		data_m <- 2^(data_m)
	}

	list(data_m = data_m, data_mzrt = data_mzrt, feat_inf = feat_inf, data_vals = dataA_vals)
}

# Helper: Select soft-threshold power
.select_power <- function(data_m, simmat, adjacencyfromsimilarity, networktype) {
	powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
	if (adjacencyfromsimilarity == TRUE) {
		sft <- try(pickSoftThreshold.fromSimilarity(similarity = simmat, powerVector = powers, verbose = 0), silent = TRUE)
	} else {
		sft <- try(pickSoftThreshold(data = data_m, dataIsExpr = TRUE, powerVector = powers, verbose = 0), silent = TRUE)
	}
	if (is(sft, "try-error")) {
		power_val <- 6
	} else {
		power_val <- sft$powerEstimate
		if (is.na(power_val) == TRUE) power_val <- 6
	}
	power_val
}

# Helper: Compute adjacency matrix
.compute_adjacency <- function(data_m, simmat, adjacencyfromsimilarity, networktype, power_val, cormethod) {
	if (adjacencyfromsimilarity == TRUE) {
		ADJdataOne <- adjacency.fromSimilarity(similarity = simmat, power = power_val, type = networktype)
	} else {
		if (cormethod == "pearson") {
			ADJdataOne <- adjacency(datExpr = data_m, type = networktype, power = power_val, corOptions = "use = 'p'")
		} else {
			ADJdataOne <- adjacency(datExpr = data_m, type = networktype, power = power_val, corOptions = "use = 'p',method='spearman'")
		}
	}
	# Drop duplicated rownames if any
	rnames_simmat <- rownames(ADJdataOne)
	dup_names <- which(duplicated(rnames_simmat) == TRUE)
	if (length(dup_names) > 0) {
		ADJdataOne <- ADJdataOne[-c(dup_names), -c(dup_names)]
	}
	ADJdataOne
}

# Helper: Distance calculation from adjacency (TOM-based)
.compute_distance <- function(ADJdataOne) {
	dissTOMCormat <- TOMdist(ADJdataOne)
	simTOMCormat <- 1 - dissTOMCormat
	list(diss = dissTOMCormat, sim = simTOMCormat)
}

# Helper: Hierarchical clustering
.perform_hclust <- function(dissTOMCormat) {
	flashClust(as.dist(dissTOMCormat), method = "complete")
}

# Helper: Module extraction from dendrogram
.extract_modules <- function(hclust_obj, dissTOMCormat, deepsplit, minclustsize, cutheight, data_m, power_val, num_nodes) {
	colorhdataOne2 <- cutreeDynamic(hclust_obj, distM = dissTOMCormat, deepSplit = deepsplit, minClusterSize = minclustsize, pamRespectsDendro = FALSE, pamStage = TRUE)
	l2colors <- levels(as.factor(colorhdataOne2))
	if (length(l2colors) > 1) {
		m2 <- try(mergeCloseModules(data_m, colors = colorhdataOne2, cutHeight = cutheight), silent = TRUE)
		if (is(m2, "try-error")) {
			mod_list <- colorhdataOne2
		} else {
			mod_list <- as.numeric(m2$colors)
		}
	} else {
		mod_list <- colorhdataOne2
	}
	mod_list
}

# Helper: Build RT-based subgroups per module
.build_rt_groups <- function(data_mzrt, data_vals, mod_list, time_step, max.rt.diff) {
	dataA_full <- cbind(data_mzrt, data_vals)
	dataA_full <- as.data.frame(dataA_full)
	t1 <- table(mod_list)
	mod_names <- as.numeric(names(t1))
	diffmatB <- lapply(1:length(mod_names), function(i) {
		groupA_num <- mod_names[i]
		subdata <- dataA_full[which(mod_list == groupA_num), ]
		subdata <- subdata[order(subdata$time), ]
		groupB <- group_by_rt_histv1(subdata, time_step, max_diff_rt = max.rt.diff, groupnum = groupA_num)
		rownames(groupB) <- NULL
		groupB
	})
	ldply(diffmatB, rbind)
}

get_peak_blocks_modulesvhclust <- function(dataA = NA, simmat = NA,
	adjacencyfromsimilarity = FALSE, time_step = 3, max.rt.diff = 10,
	outloc, column.rm.index = NA, cor.thresh = NA, deepsplit = 2,
	minclustsize = 20, cutheight = 0.2, cormethod = "spearman",
	networktype = "unsigned", num_nodes = 2, step1log2scale = TRUE,
	mycl_metabs = NA) {

	# Basic checks
	.check_sanity(dataA, simmat, adjacencyfromsimilarity)

	# Prepare matrices
	prep <- .prepare_feature_matrix(dataA, column.rm.index, step1log2scale)
	data_m <- prep$data_m
	data_mzrt <- prep$data_mzrt
	data_vals <- prep$data_vals

	# Attempt blockwiseModules first (as in original)
	do_stepwise <- 0
	powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
	sft_try <- try(pickSoftThreshold(data = data_m, dataIsExpr = TRUE, powerVector = powers, verbose = 0), silent = TRUE)
	if (is(sft_try, "try-error")) {
		power_val <- 6
	} else {
		power_val <- sft_try$powerEstimate
		if (is.na(power_val) == TRUE) power_val <- 6
	}

	netclassA <- try(blockwiseModules(datExpr = (data_m),
		checkMissingData = FALSE, blocks = mycl_metabs,
		maxBlockSize = 5000, blockSizePenaltyPower = 100,
		randomSeed = 12345, loadTOM = FALSE, corType = "pearson",
		maxPOutliers = 1, quickCor = 0, pearsonFallback = "individual",
		cosineCorrelation = FALSE, power = power_val,
		networkType = "unsigned", TOMType = "signed",
		TOMDenom = "min", getTOMs = NULL, saveTOMs = FALSE,
		saveTOMFileBase = "blockwiseTOM", deepSplit = deepsplit,
		detectCutHeight = NULL, minModuleSize = minclustsize,
		maxCoreScatter = NULL, minGap = NULL, maxAbsCoreScatter = NULL,
		minAbsGap = NULL, minSplitHeight = NULL, minAbsSplitHeight = NULL,
		useBranchEigennodeDissim = FALSE, minBranchEigennodeDissim = cutheight,
		stabilityLabels = NULL, minStabilityDissim = NULL,
		pamStage = TRUE, pamRespectsDendro = FALSE, reassignThreshold = 1e-06,
		minCoreKME = 0.5, minCoreKMESize = minclustsize/3,
		minKMEtoStay = 0.3, mergeCutHeight = cutheight,
		impute = TRUE, trapErrors = FALSE, numericLabels = FALSE,
		nThreads = num_nodes, verbose = 0, indent = 0),
		silent = TRUE)

	if (is(netclassA, "try-error")) {
		do_stepwise <- 1
	} else {
		n1 <- unlist(netclassA$colors)
		n2 <- as.data.frame(n1, stringsAsFactors = TRUE)
		mod_list <- as.numeric(n2[, 1])
		Alldegrees1 <- softConnectivity(datExpr = data_m, power = power_val, minNSamples = 2)
		Alldegrees1 <- cbind(Alldegrees1, Alldegrees1, Alldegrees1, Alldegrees1)
	}

	if (do_stepwise == 1) {
		# Select power based on requested pathway
		power_val <- .select_power(data_m, simmat, adjacencyfromsimilarity, networktype)
		# Build adjacency
		ADJdataOne <- .compute_adjacency(data_m, simmat, adjacencyfromsimilarity, networktype, power_val, cormethod)
		# Distance calculation
		dist_list <- .compute_distance(ADJdataOne)
		dissTOMCormat <- dist_list$diss
		# Hierarchical clustering
		hierTOMCormat <- .perform_hclust(dissTOMCormat)
		# optional dendrogram plot (disabled to match original behavior)
		if (FALSE) {
			par(mfrow = c(1, 1))
			pdf("plot.pdf")
			plot(hierTOMCormat, labels = FALSE, main = "Dendrogram")
			dev.off()
		}
		# Module extraction
		mod_list <- .extract_modules(hierTOMCormat, dissTOMCormat, deepsplit, minclustsize, cutheight, data_m, power_val, num_nodes)
		# Intramodular connectivity
		Alldegrees1 <- intramodularConnectivity(ADJdataOne, mod_list)
	}

	# Build RT-based groups and combine outputs
	diffmatB <- .build_rt_groups(data_mzrt = data_mzrt, data_vals = data_vals, mod_list = mod_list, time_step = time_step, max.rt.diff = max.rt.diff)
	diffmatB <- as.data.frame(diffmatB)
	diffmatB <- cbind(Alldegrees1[, c(1:4)], diffmatB)

	return(diffmatB)
}
