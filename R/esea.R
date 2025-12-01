#' Compute Edge Set Enrichment Analysis (ESEA)
#'
#' Perform an edge-based enrichment analysis on a ranked list of edge correlation
#' scores. This function supports two permutation schemes (`edge.labels` and
#' `gene.labels`) and returns normalized enrichment scores (NES), nominal
#' p-values, FDR-adjusted q-values, and detailed running-enrichment statistics for
#' significant or top-ranked pathways.
#'
#' @param EdgeCorScore A **named numeric vector** of edge scores (e.g.,
#'   correlation scores). Names must be edge identifiers in the form
#'   `"node1|node2"` (the function automatically converts `";"` to `"|"`).
#' @param pathway_edge_map A **list of character vectors**, where each element is
#'   a set of edges (same `"node1|node2"` format) defining a pathway or module.
#' @param weighted.score.type Numeric. Weighting exponent \eqn{\alpha} applied to
#'   edge scores when computing running enrichment statistics.  
#'   *0 = unweighted*, *1 = |score|*, *2 = score²*, or any non-negative numeric.
#' @param gs.size.threshold.min Minimum allowed size of an edge set after
#'   intersecting with available edges (default: 1).
#' @param gs.size.threshold.max Maximum allowed size of an edge set after
#'   intersecting with available edges (default: 1000).
#' @param reshuffling.type Character. Permutation scheme:
#'   *`"edge.labels"`* (shuffle ranks) or *`"gene.labels"`* (simulate a new score
#'   vector and re-rank).  
#'   Default: `"edge.labels"`.
#' @param nperm Integer. Number of permutations (default: 100).
#' @param p.val.threshold Numeric. Optional threshold to force inclusion of sets
#'   into the detailed output (default: -1 for disabled).
#' @param FDR.threshold Numeric. Threshold on FDR q-values to select sets for
#'   detailed output (default: 0.05).
#' @param topgs Integer. Number of top positively and negatively enriched sets to
#'   include in the detailed output (default: 1).
#' @param decreasing Logical. Whether to rank edge scores in decreasing order
#'   (default: TRUE).
#'
#' @details
#' The function:
#' * Normalizes all edge identifiers to use the `"|"` separator.  
#' * Ranks edges and intersects each pathway edge set with available edges.  
#' * Computes observed enrichment scores (ES) and auxiliary metrics.  
#' * Performs streaming permutation tests to estimate nominal p-values.  
#' * Computes normalized enrichment scores (NES) and FDR q-values.  
#' * Returns summary tables for positively and negatively enriched pathways, and
#'   detailed running-enrichment statistics for selected sets.
#'
#' The detailed output includes RES curves and core-enrichment membership of each
#' edge contributing to the enrichment signal.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{SUMMARY.RESULTS}}{
#'     A list with:
#'     \itemize{
#'       \item \code{SUMMARY.RESULTS.High_portion_of_the_rank}: data frame of
#'         positively enriched pathways (sorted by NES).
#'       \item \code{SUMMARY.RESULTS.Low_portion_of_the_rank}: data frame of
#'         negatively enriched pathways (sorted by NES).
#'     }
#'   }
#'   \item{\code{"Pathway results"}}{
#'     A named list of data frames, one per selected pathway, containing detailed
#'     running enrichment statistics (RES, edge IDs, positions, core enrichment).
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' scores <- rnorm(100)
#' names(scores) <- paste0("n", 1:100, "|n", 2:101)
#'
#' pathways <- list(
#'   P1 = sample(names(scores), 10),
#'   P2 = sample(names(scores), 15)
#' )
#'
#' res <- compute_esea(
#'   EdgeCorScore = scores,
#'   pathway_edge_map = pathways,
#'   nperm = 50
#' )
#'
#' res$SUMMARY.RESULTS$SUMMARY.RESULTS.High_portion_of_the_rank
#' }
#'
#' @export
compute_esea <- function(
    EdgeCorScore,
    pathway_edge_map,
    weighted.score.type = 1,
    gs.size.threshold.min = 1,
    gs.size.threshold.max = 1000,
    reshuffling.type = "edge.labels", # or "gene.labels"
    nperm = 100,
    p.val.threshold = -1,
    FDR.threshold = 0.05,
    topgs = 1,
    decreasing = TRUE
) {
  message("Running ESEA Analysis (list-based, node1|node2)…")
  
  stopifnot(is.numeric(EdgeCorScore), !is.null(names(EdgeCorScore)))
  stopifnot(is.list(pathway_edge_map))
  if (length(EdgeCorScore) == 0) stop("EdgeCorScore is empty.")
  if (anyNA(EdgeCorScore)) stop("EdgeCorScore contains NA values.")
  if (!is.logical(decreasing) || length(decreasing) != 1) stop("decreasing must be TRUE/FALSE.")
  
  # ---- Enforce/normalize the '|' separator everywhere (idempotent) ----
  names(EdgeCorScore) <- gsub(";", "|", names(EdgeCorScore), fixed = TRUE)
  pathway_edge_map <- lapply(pathway_edge_map, function(v) unique(gsub(";", "|", v, fixed = TRUE)))
  
  # ---- Helper: fast ES using tag positions (ranks) ----
  es_from_tag_locs <- function(tag_locs, N, alpha, correl_sorted) {
    Nh <- length(tag_locs)
    if (Nh == 0L) return(list(ES = 0, arg.ES = NA_integer_, tag_before = 0L, tag_after = 0L))
    
    if (alpha == 0) {
      w <- rep(1, Nh)
    } else if (alpha == 1) {
      w <- abs(correl_sorted[tag_locs])
    } else if (alpha == 2) {
      v <- correl_sorted[tag_locs]
      w <- abs(v * v)
    } else {
      w <- abs(correl_sorted[tag_locs])^alpha
    }
    w <- w / sum(w)
    
    Nm <- N - Nh
    step_miss <- 1 / Nm
    d <- c(tag_locs[1] - 1L, diff(tag_locs) - 1L)
    miss_cumsum <- cumsum(d * step_miss)
    hit_cumsum  <- cumsum(w)
    RES_peaks   <- hit_cumsum - miss_cumsum
    RES_valleys <- RES_peaks - w
    
    maxES <- max(RES_peaks)
    minES <- min(RES_valleys)
    if (maxES > -minES) {
      ES <- signif(maxES, 5)
      arg.ES <- tag_locs[ which.max(RES_peaks) ]
      tag_before <- sum(tag_locs <= arg.ES)
      tag_after  <- Nh - tag_before
    } else {
      ES <- signif(minES, 5)
      arg.ES <- tag_locs[ which.min(RES_valleys) ]
      tag_before <- sum(tag_locs < arg.ES)
      tag_after  <- Nh - tag_before
    }
    list(ES = ES, arg.ES = arg.ES, tag_before = tag_before, tag_after = tag_after)
  }
  
  # ---- Prepare ranked list ----
  N <- length(EdgeCorScore)
  ord <- order(EdgeCorScore, decreasing = decreasing)
  correl_sorted <- EdgeCorScore[ord]
  edge_labels   <- names(EdgeCorScore)
  labels_sorted <- edge_labels[ord]
  rank_pos <- seq_len(N); names(rank_pos) <- labels_sorted
  
  # ---- Build intersected edge sets as rank vectors ----
  set_names <- names(pathway_edge_map)
  if (is.null(set_names)) set_names <- paste0("SET_", seq_along(pathway_edge_map))
  
  tag_locs_list <- vector("list", length(pathway_edge_map))
  kept_names    <- character(0)
  kept_sizes    <- integer(0)
  
  for (i in seq_along(pathway_edge_map)) {
    edges <- unique(pathway_edge_map[[i]])
    locs <- unname(rank_pos[edges])
    locs <- locs[!is.na(locs)]
    s <- length(locs)
    if (s >= gs.size.threshold.min && s <= gs.size.threshold.max) {
      tag_locs_list[[length(kept_names) + 1L]] <- sort(as.integer(locs))
      kept_names <- c(kept_names, set_names[i])
      kept_sizes <- c(kept_sizes, s)
    }
  }
  
  Ng <- length(kept_names)
  if (Ng == 0L) {
    warning("No edge sets passed the size thresholds after intersecting with available edges.")
    return(list(
      SUMMARY.RESULTS = list(
        SUMMARY.RESULTS.High_portion_of_the_rank = data.frame(),
        SUMMARY.RESULTS.Low_portion_of_the_rank  = data.frame()
      ),
      "Pathway results" = list()
    ))
  }
  
  # ---- Observed ES per set ----
  alpha <- as.numeric(weighted.score.type)
  Obs.ES       <- numeric(Ng)
  Obs.arg.ES   <- integer(Ng)
  tag_before   <- integer(Ng)
  tag_after    <- integer(Ng)
  for (i in seq_len(Ng)) {
    res <- es_from_tag_locs(tag_locs_list[[i]], N, alpha, correl_sorted)
    Obs.ES[i]     <- res$ES
    Obs.arg.ES[i] <- ifelse(is.na(res$arg.ES), 0L, res$arg.ES)
    tag_before[i] <- res$tag_before
    tag_after[i]  <- res$tag_after
  }
  
  # ---- Heuristics (signal strength, tag/edge fractions) ----
  edge_frac <- ifelse(Obs.ES >= 0,
                      pmax(Obs.arg.ES, 1L) / N,
                      pmax(N - Obs.arg.ES + 1L, 1L) / N)
  tag_frac  <- ifelse(Obs.ES >= 0,
                      tag_before / kept_sizes,
                      tag_after  / kept_sizes)
  signal_strength <- signif(tag_frac * (1 - edge_frac) * (N / (N - kept_sizes)), 3)
  
  # ---- Permutations (streaming; memory-light) ----
  ge_counts <- integer(Ng)      # for ES >= Obs.ES (pos side)
  le_counts <- integer(Ng)      # for ES <= Obs.ES (neg side)
  pos_sum <- numeric(Ng); pos_cnt <- integer(Ng)
  neg_sum_abs <- numeric(Ng); neg_cnt <- integer(Ng)
  
  tag_idx <- tag_locs_list
  
  compute_all_ES <- function(rank_map, correl_vec) {
    out <- numeric(Ng)
    for (i in seq_len(Ng)) {
      locs <- rank_map[ labels_sorted[tag_idx[[i]]] ]
      locs <- sort(as.integer(unname(locs)))
      out[i] <- es_from_tag_locs(locs, N, alpha, correl_vec)$ES
    }
    out
  }
  
  if (reshuffling.type == "edge.labels") {
    for (r in seq_len(nperm)) {
      perm <- sample.int(N, N, replace = FALSE)
      rank_map <- setNames(perm, labels_sorted)
      es_perm <- compute_all_ES(rank_map, correl_sorted)
      pos_mask <- es_perm >= 0
      pos_sum[pos_mask] <- pos_sum[pos_mask] + es_perm[pos_mask]
      pos_cnt[pos_mask] <- pos_cnt[pos_mask] + 1L
      neg_mask <- !pos_mask
      neg_sum_abs[neg_mask] <- neg_sum_abs[neg_mask] + abs(es_perm[neg_mask])
      neg_cnt[neg_mask] <- neg_cnt[neg_mask] + 1L
      ge_counts <- ge_counts + as.integer(es_perm >= Obs.ES)
      le_counts <- le_counts + as.integer(es_perm <= Obs.ES)
    }
  } else if (reshuffling.type == "gene.labels") {
    for (r in seq_len(nperm)) {
      sim <- rnorm(N, mean = mean(EdgeCorScore), sd = stats::sd(EdgeCorScore))
      ord_sim <- order(sim, decreasing = decreasing)
      correl_sim <- sim[ord_sim]
      labels_sim <- edge_labels[ord_sim]
      rank_map <- setNames(seq_len(N), labels_sim)
      es_perm <- compute_all_ES(rank_map, correl_sim)
      pos_mask <- es_perm >= 0
      pos_sum[pos_mask] <- pos_sum[pos_mask] + es_perm[pos_mask]
      pos_cnt[pos_mask] <- pos_cnt[pos_mask] + 1L
      neg_mask <- !pos_mask
      neg_sum_abs[neg_mask] <- neg_sum_abs[neg_mask] + abs(es_perm[neg_mask])
      neg_cnt[neg_mask] <- neg_cnt[neg_mask] + 1L
      ge_counts <- ge_counts + as.integer(es_perm >= Obs.ES)
      le_counts <- le_counts + as.integer(es_perm <= Obs.ES)
    }
  } else {
    stop('reshuffling.type must be "edge.labels" or "gene.labels".')
  }
  
  # ---- Stats: p-values, NES, FDR ----
  p_nom <- numeric(Ng)
  pos_side <- Obs.ES >= 0
  p_nom[pos_side]  <- ge_counts[pos_side] / nperm
  p_nom[!pos_side] <- le_counts[!pos_side] / nperm
  # p_nom[pos_side]  <- (ge_counts[pos_side] + 1) / (nperm + 1) # If you don't want the FDR to be equal to zero
  # p_nom[!pos_side] <- (le_counts[!pos_side] + 1) / (nperm + 1)
  # p_nom <- signif(p_nom, 5)
  
  pos_mean <- ifelse(pos_cnt > 0, pos_sum / pmax(pos_cnt, 1L), NA_real_)
  neg_mean <- ifelse(neg_cnt > 0, neg_sum_abs / pmax(neg_cnt, 1L), NA_real_)
  NES <- numeric(Ng)
  NES[pos_side]  <- Obs.ES[pos_side]  / pos_mean[pos_side]
  NES[!pos_side] <- Obs.ES[!pos_side] / neg_mean[!pos_side]
  NES <- signif(NES, 5)
  
  FDR <- p.adjust(p_nom, method = "fdr")
  
  # ---- Summaries ----
  report <- data.frame(
    GS    = kept_names,
    SIZE  = kept_sizes,
    ES    = signif(Obs.ES, 5),
    NES   = NES,
    NOM_pval = p_nom,
    FDR_qval = FDR,
    Tag_pct  = signif(tag_frac, 3),
    Edge_pct = signif(edge_frac, 3),
    Signal   = signal_strength,
    stringsAsFactors = FALSE
  )
  
  pos_idx <- which(NES >= 0 & !is.na(NES))
  neg_idx <- which(NES <  0 & !is.na(NES))
  report_pos <- if (length(pos_idx)) report[pos_idx[order(-report$NES[pos_idx])], ] else report[0, ]
  report_neg <- if (length(neg_idx)) report[neg_idx[order(report$NES[neg_idx])], ]  else report[0, ]
  
  # ---- Choose sets for detailed output ----
  top_pos <- if (nrow(report_pos)) head(rownames(report_pos), n = min(topgs, nrow(report_pos))) else character(0)
  top_neg <- if (nrow(report_neg)) head(rownames(report_neg), n = min(topgs, nrow(report_neg))) else character(0)
  keep_detail <- rep(FALSE, Ng)
  if (p.val.threshold >= 0) keep_detail <- keep_detail | (p_nom <= p.val.threshold)
  keep_detail <- keep_detail | (FDR <= FDR.threshold)
  if (length(top_pos)) keep_detail[as.integer(top_pos)] <- TRUE
  if (length(top_neg)) keep_detail[as.integer(top_neg)] <- TRUE
  
  # ---- Detailed per-set outputs (RES only for chosen sets) ----
  detailed <- list()
  chosen_idx <- which(keep_detail)
  if (length(chosen_idx)) {
    if (alpha == 0) {
      w_all <- rep(1, N)
    } else if (alpha == 1) {
      w_all <- abs(correl_sorted)
    } else if (alpha == 2) {
      w_all <- abs(correl_sorted * correl_sorted)
    } else {
      w_all <- abs(correl_sorted)^alpha
    }
    
    for (i in chosen_idx) {
      set_name <- kept_names[i]
      locs <- tag_locs_list[[i]]
      Nh <- length(locs)
      Nm <- N - Nh
      
      tag_indicator <- integer(N); tag_indicator[locs] <- 1L
      norm_tag <- 1 / sum(w_all[locs])
      norm_miss <- 1 / Nm
      RES <- cumsum(tag_indicator * w_all * norm_tag - (1 - tag_indicator) * norm_miss)
      
      ES <- Obs.ES[i]
      argES <- Obs.arg.ES[i]
      phen_tag <- if (ES >= 0) "High_part_of_the_rank" else "Low_part_of_the_rank"
      
      if (ES >= 0) {
        order_ranks <- seq_len(N)
        core_yes <- which(locs <= argES)
      } else {
        order_ranks <- N:1
        core_yes <- which(locs > argES)
      }
      idx_in_order <- which(tag_indicator[order_ranks] == 1L)
      ranks_kept <- order_ranks[idx_in_order]
      edge_ids <- labels_sorted[ranks_kept]
      df <- data.frame(
        `#` = seq_along(ranks_kept),
        EdgeID = edge_ids,
        `List Loc` = ranks_kept,
        EdgeCorScore = signif(correl_sorted[ranks_kept], 3),
        RES = signif(RES[ranks_kept], 3),
        CORE_ENRICHMENT = ifelse(ranks_kept %in% locs[core_yes], "YES", "NO"),
        stringsAsFactors = FALSE
      )
      detailed[[paste0(set_name, ".", phen_tag)]] <- df
    }
  }
  
  list(
    SUMMARY.RESULTS = list(
      SUMMARY.RESULTS.High_portion_of_the_rank = report_pos,
      SUMMARY.RESULTS.Low_portion_of_the_rank  = report_neg
    ),
    "Pathway results" = detailed
  )
}
