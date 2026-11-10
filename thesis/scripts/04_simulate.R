# 04_simulate.R — main generator of the pseudo-datasets (PROTOCOLLO §5.2, §6 T4; ADEMP §2.1–2.4, §8).
# Generation only: scripts/05_deconvolve.R (CPU methods) reads the files written here (cell2location on
# the GPU reads the same files).
#
# Design (fixed by ADEMP v1.2; implementation choices written 2026-11-05 before the first run):
#   datasets of one replicate, per tissue: anchor + arms E, C−, E′, O, P (G reserved: it enters only
#   by a DECISIONS row, ADEMP §2.1) × cells per spot 2 / 10 × depth 175 / 4,900 UMI per cell
#   = 24 datasets per tissue (72 over the three tissues; ADEMP §2.2); 500 spots each.
#   - E: μ_E (training) on the spot types; C−: μ_E with gene labels permuted within each type, one new
#     permutation per replicate (DECISIONS 2026-10-30); E′: μ of the E′ donor; O: μ of the held-out
#     donor (mu_oracle()); P: μ^P and φ^P drawn in each replicate (R/arm_p.R); every arm with the
#     common NB noise model, φ_g from training (P: φ by rank); anchor = real held-out cells.
#   - Common random numbers (ADEMP §2.4): in a replicate and a cells-per-spot level, one composition
#     table and one cell table (type, library size) for every arm, and the same compositions for the
#     anchor. Choice: the two depth levels thin the SAME raw spots (one generation, two targets), so
#     the depth comparison is also paired.
#   - One μ^P, one P reference and one C− permutation per replicate and tissue, shared by its four
#     conditions: the P reference is trained once per tissue and replicate (DECISIONS 2026-11-04,
#     cell2location budget: "one complete replicate = 72 datasets + the training of the P reference
#     in each tissue").
#   - The P reference is not stored as a matrix (brain: ~27,000 cells per replicate): world_p.rds
#     stores μ^P, φ^P, the markers and the RNG state of the reference, and reference_p_from_world()
#     (R/arm_p.R) rebuilds it exactly (ADEMP §8: "or their seeds and states").
# Seeds (ADEMP §8): study seed 42, L'Ecuyer-CMRG. Tissue t (colon 1, lymphoid 2, brain 3) uses the
# streams (t − 1) × 1,000 + rep, so tissues never share a stream. Inside a replicate every random step
# has its own sub-stream, so adding or cutting an arm (G, or a cut by the cut rule) leaves every
# other dataset unchanged: 1 = μ^P and φ^P ties; 2 = P reference; 3 = C− permutation;
# 100 j = compositions and cells of cells-per-spot level j; 100 j + 10 a = counts of dataset a
# (anchor 1: held-out cells drawn; E 2, C− 3, E′ 4, O 5, P 6, G 7 reserved); 100 j + 10 a + d =
# thinning of dataset a to depth level d. The replicate stream is saved with each output.
#
# Pre-registered checks (written 2026-11-05 before the first run; the job stops if one fails):
#   S1. truth (cell shares): rows sum to 1 (1e-12), zero on reference types absent from the spots;
#       UMI-share truth (arms: drawn L; anchor: totals of the real cells) rows sum to 1 (1e-9).
#   S2. every dataset has the same 500 spot names; real arms and the anchor on the reference genes,
#       P on its own genes P00001…; arms share one cell table (one truth per condition).
#   S3. after thinning no spot exceeds the target; a spot at or above the target ends exactly at it,
#       a spot below it is unchanged.
#   S4. C− μ has the column sums of E μ (1e-9).
#   S5. the held-out cells come from one donor, the held-out donor; no reference cell comes from the
#       held-out or the E′ donor; held-out, E′ and training donors are disjoint; μ_E′ and μ_O columns
#       in (0.9, 1].
#   S6. μ^P: genes P00001…, columns sum to 1 (1e-12), markers per type = calibration counts.
# Check mode (run_mode = "check", replicates 1–2, nothing written to the output folder):
#   R1. replicate 1 generated twice from its stream: identical datasets.
#   R2. replicate 1 without P and E′: the anchor, E, C− and O datasets identical to the full run
#       (sub-streams independent, so adding G later changes nothing else).
#   R3. P reference rebuilt twice from world_p: identical; same cells, types and donors as the real
#       reference; column totals / L mean within 1% of 1.
#   R4. replicates 1 and 2 differ in μ^P, in the C− permutation and in the compositions.
#   K1. counter-proof (condition 10 cells, 4,900 UMI per cell, replicate 1): mean over spots of the
#       Pearson correlation (log1p CP10k) between spot s of a dataset and spot s of the anchor:
#       expected O > E > C− and O > E′ > C−; E vs E′ reported without expectation (P has other genes).
# Reported: seconds per replicate, MB on disk per replicate (input to the cost check, ADEMP §6);
# spots at and below the target per dataset (DECISIONS 2026-11-02); reuse per held-out cell per type,
# observed vs expected 500 × cells × p_k / n_k (ADEMP §2.2).
#
# The main run (reference "main") reads the held-out donor and is blocked until the ADEMP is approved:
# set ADEMP_APPROVED to Luca's approval date (a visible, committed change).
#
# Usage (working directory = thesis/):
#   Rscript scripts/04_simulate.R <tissue> [pilot|main] [n_rep] [workers] [run|check]
# Output: <TESI_DATA>/sim/<pilot|main>/<tissue>/r<rep>/{world_p.rds, c<cells>_d<depth>.rds};
# summaries <tissue>_sim_depth.csv, <tissue>_sim_reuse.csv, <tissue>_sim_jobs.csv next to the folder.

ADEMP_APPROVED <- NA_character_

args <- commandArgs(trailingOnly = TRUE)
tissue <- if (length(args) >= 1) args[1] else "colon"
ref_mode <- if (length(args) >= 2) args[2] else "pilot"
n_rep <- if (length(args) >= 3) as.integer(args[3]) else 81L
workers <- if (length(args) >= 4) as.integer(args[4]) else 8L
run_mode <- if (length(args) >= 5) args[5] else "run"
stopifnot(tissue %in% c("colon", "lymphoid", "brain"), ref_mode %in% c("pilot", "main"),
          run_mode %in% c("run", "check"), n_rep >= 1)
if (ref_mode == "main" && is.na(ADEMP_APPROVED)) stop("Main run blocked: the ADEMP is not approved yet")

data_dir <- Sys.getenv("TESI_DATA")
stopifnot("TESI_DATA non impostata" = nzchar(data_dir), dir.exists(data_dir))
stopifnot("Run from the thesis/ folder" = file.exists("thesis.Rproj"))
suppressPackageStartupMessages(library(Matrix))   # S4 subsetting of the sparse counts
for (f in c("R/simulate_spots.R", "R/arm_p.R")) source(f)
ref_file <- file.path(data_dir, "derived", "reference", sprintf("%s_%s.rds", tissue, ref_mode))
calib_file <- file.path(data_dir, "derived", "p_calibration", sprintf("%s_p_calibration.rds", tissue))
out_dir <- file.path(data_dir, "sim", ref_mode, tissue)
log_file <- file.path(dirname(out_dir), paste0(tissue, "_sim_progress.log"))

study_seed <- 42
n_spots <- 500
levels_cells <- c(2, 10)
levels_depth <- c(175, 4900)   # UMI per cell; target = cells x UMI per cell (ADEMP §2.2)
slot <- c(anchor = 1, E = 2, Cminus = 3, Eprime = 4, O = 5, P = 6, G = 7)
arms_all <- c("anchor", "E", "Cminus", "Eprime", "O", "P")

R0 <- readRDS(ref_file)
calib <- readRDS(calib_file)
ts <- R0$types_spot; types <- R0$types_ref
n_train <- lengths(R0$L_by_type[ts])
# Dirichlet concentration c per tissue and cells per spot, alpha_k = c p_k (ADEMP §2.2; DECISIONS
# 2026-11-06, Luca): values in dirichlet_c(), R/simulate_spots.R. One alpha per cells level.
alpha_by_n <- lapply(levels_cells, function(nc) dirichlet_c(tissue, nc) * n_train / sum(n_train))

# S5 and static μ
hd <- unique(R0$held_meta$donor_id)
if (length(hd) != 1 || hd != R0$donors$heldout) stop("S5: held-out cells not from the held-out donor only")
if (any(R0$ref_meta$donor_id %in% c(R0$donors$heldout, R0$donors$eprime))) stop("S5: held-out/E′ cell in the reference")
if (length(intersect(R0$donors$train, c(R0$donors$heldout, R0$donors$eprime))) || R0$donors$heldout == R0$donors$eprime)
  stop("S5: donors of the splits not disjoint")
mu_fix <- list(E = R0$mu_E[, ts, drop = FALSE], Eprime = R0$mu_Ep[, ts, drop = FALSE], O = mu_oracle(R0$mu_held, ts))
cs <- colSums(mu_fix$Eprime)
if (anyNA(mu_fix$Eprime) || any(cs <= 0.9 | cs > 1 + 1e-9)) stop("S5: μ_E′ columns out of (0.9, 1]")
stopifnot(identical(rownames(R0$mu_Ep), rownames(R0$mu_E)), identical(rownames(R0$mu_held), rownames(R0$mu_E)),
          all(types %in% names(calib$markers)))
mean_real <- pooled_mean(R0$mu_E, R0$L_by_type)
held_ok <- R0$held_meta$cell_type %in% ts
held <- R0$held[, held_ok]; held_type <- R0$held_meta$cell_type[held_ok]
held_L <- Matrix::colSums(held)
genes <- rownames(R0$mu_E)

# Streams
RNGkind("L'Ecuyer-CMRG")
set.seed(study_seed)
s <- .Random.seed
for (k in seq_len((match(tissue, c("colon", "lymphoid", "brain")) - 1) * 1000)) s <- parallel::nextRNGStream(s)
rep_seed <- vector("list", n_rep)
for (i in seq_len(n_rep)) { s <- parallel::nextRNGStream(s); rep_seed[[i]] <- s }
substream <- function(seed, k) { for (j in seq_len(k)) seed <- parallel::nextRNGSubStream(seed); seed }
use_stream <- function(seed) assign(".Random.seed", seed, envir = .GlobalEnv)

# μ^P and φ^P of one replicate (sub-stream 1); the reference state is sub-stream 2.
make_world_p <- function(sr) {
  use_stream(substream(sr, 1))
  P <- simulate_mu_p(calib, length(genes), types)
  phi <- phi_by_rank(R0$phi, mean_real, P$base)
  if (!all(grepl("^P\\d{5}$", rownames(P$mu))) || any(abs(colSums(P$mu) - 1) > 1e-12) ||
      !all(lengths(P$markers)[types] == lengths(calib$markers)[types])) stop("S6: μ^P")
  list(mu = P$mu, phi = phi, base = P$base, markers = P$markers, ref_seed = substream(sr, 2))
}

# One replicate: returns the world of P and one entry per condition (counts of every dataset, truths,
# depth table, reuse table).
make_replicate <- function(i, arms = arms_all) {
  sr <- rep_seed[[i]]
  world <- if ("P" %in% arms) make_world_p(sr)
  mu_C <- NULL
  if ("Cminus" %in% arms) {
    use_stream(substream(sr, 3)); mu_C <- permute_mu(mu_fix$E)
    if (any(abs(colSums(mu_C) - colSums(mu_fix$E)) > 1e-9)) stop("S4")
  }
  res <- list()
  for (j in seq_along(levels_cells)) {
    nc <- levels_cells[j]
    use_stream(substream(sr, 100 * j))
    comp <- draw_compositions(n_spots, nc, ts, alpha_by_n[[j]])
    cells <- draw_cells(comp, R0$L_by_type[ts])
    sp <- rownames(comp)
    truth <- matrix(0, n_spots, length(types), dimnames = list(sp, types))
    truth[, ts] <- comp / nc
    if (any(abs(rowSums(truth) - 1) > 1e-12) || any(truth[, setdiff(types, ts)] != 0)) stop("S1")
    raw <- list(); truth_umi <- list(); reuse <- NULL
    for (a in arms) {
      use_stream(substream(sr, 100 * j + 10 * slot[[a]]))
      if (a == "anchor") {
        an <- anchor_spots(held, held_type, comp)
        raw$anchor <- an$counts
        truth_umi$anchor <- umi_truth(an$spot, held_type[an$cells_used], held_L[an$cells_used], sp, types)
        n_held <- table(factor(held_type, ts))
        drawn <- table(factor(held_type[an$cells_used], ts))
        reuse <- data.frame(type = ts, n_held = as.vector(n_held), drawn = as.vector(drawn),
                            reuse_obs = as.vector(drawn / n_held),
                            reuse_exp = n_spots * nc * (n_train / sum(n_train)) / as.vector(n_held))
      } else {
        mu <- switch(a, Cminus = mu_C, P = world$mu[, ts, drop = FALSE], mu_fix[[a]])
        phi <- if (a == "P") world$phi else R0$phi
        raw[[a]] <- simulate_arm(mu, phi, cells, sp)
      }
    }
    tu_arm <- umi_truth(cells$spot, cells$type, cells$L, sp, types)
    for (a in setdiff(arms, "anchor")) truth_umi[[a]] <- tu_arm
    if (any(vapply(truth_umi, function(t) any(abs(rowSums(t) - 1) > 1e-9), logical(1)))) stop("S1 (UMI truth)")
    for (a in arms) {
      x <- raw[[a]]
      g_ok <- if (a == "P") identical(rownames(x), rownames(world$mu)) else identical(rownames(x), genes)
      if (ncol(x) != n_spots || !identical(colnames(x), sp) || !g_ok) stop("S2: ", a)
    }
    for (d in seq_along(levels_depth)) {
      target <- nc * levels_depth[d]
      counts <- list(); dep <- list()
      for (a in arms) {
        use_stream(substream(sr, 100 * j + 10 * slot[[a]] + d))
        x <- thin_to_depth(raw[[a]], target)
        u0 <- Matrix::colSums(raw[[a]]); u1 <- Matrix::colSums(x)
        if (any(u1 > target) || any(u1[u0 >= target] != target) || any(u1[u0 < target] != u0[u0 < target]))
          stop("S3: ", a)
        counts[[a]] <- x
        dep[[a]] <- data.frame(dataset = a, target = target, n_at_target = sum(u0 >= target),
                               n_below = sum(u0 < target), median_raw_umi = median(u0))
      }
      key <- sprintf("c%02d_d%04d", nc, levels_depth[d])
      res[[key]] <- list(counts = counts, truth = truth, truth_umi = truth_umi[arms],
                         depth = cbind(tissue = tissue, rep = i, n_cells = nc, depth = levels_depth[d], do.call(rbind, dep)),
                         reuse = cbind(tissue = tissue, rep = i, n_cells = nc, reuse),
                         meta = list(tissue = tissue, ref_mode = ref_mode, rep = i, n_cells = nc,
                                     depth = levels_depth[d], target = target, arms = arms, rep_seed = sr,
                                     types_spot = ts, types_ref = types, study_seed = study_seed, alpha = alpha_by_n[[j]]))
    }
  }
  list(world = world, conditions = res)
}

save_atomic <- function(x, f) { tmp <- paste0(f, ".tmp"); saveRDS(x, tmp); file.rename(tmp, f) }
log_line <- function(msg) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n", file = log_file, append = TRUE)

run_job <- function(i) {
  rd <- file.path(out_dir, sprintf("r%03d", i))
  files <- c(file.path(rd, "world_p.rds"),
             file.path(rd, sprintf("c%02d_d%04d.rds", rep(levels_cells, each = 2), rep(levels_depth, 2))))
  if (all(file.exists(files))) return(list(rep = i, skipped = TRUE))   # restart: replicate already on disk
  dir.create(rd, recursive = TRUE, showWarnings = FALSE)
  log_line(sprintf("start rep %d", i))
  t0 <- proc.time()[["elapsed"]]
  r <- make_replicate(i)
  save_atomic(r$world, files[1])
  for (k in names(r$conditions)) save_atomic(r$conditions[[k]], file.path(rd, paste0(k, ".rds")))
  sec <- proc.time()[["elapsed"]] - t0
  mb <- sum(file.size(files)) / 2^20
  log_line(sprintf("done rep %d in %.1f min, %.0f MB", i, sec / 60, mb))
  list(rep = i, skipped = FALSE, seconds = sec, mb = mb,
       depth = do.call(rbind, lapply(r$conditions, `[[`, "depth")),
       reuse = do.call(rbind, lapply(r$conditions, `[[`, "reuse")))
}

# ---------------------------------------------------------------------------------------------------
if (run_mode == "check") {
  res <- character(0)
  ok <- function(cond, msg) { r <- if (isTRUE(cond)) "PASS" else "FAIL"; res <<- c(res, r); cat(r, msg, "\n") }
  t0 <- proc.time()[["elapsed"]]
  r1 <- make_replicate(1)
  sec1 <- proc.time()[["elapsed"]] - t0
  tmp <- tempfile(fileext = ".rds"); saveRDS(r1$world, tmp); mb <- file.size(tmp)
  for (k in names(r1$conditions)) { saveRDS(r1$conditions[[k]], tmp); mb <- mb + file.size(tmp) }
  unlink(tmp)
  cat(sprintf("Replicate 1: %.1f s, %.0f MB on disk (24 datasets of one tissue, no P reference)\n", sec1, mb / 2^20))
  same_counts <- function(a, b, ds) all(vapply(names(a$conditions), function(k)
    all(vapply(ds, function(d) identical(a$conditions[[k]]$counts[[d]], b$conditions[[k]]$counts[[d]]), logical(1))), logical(1)))
  r1b <- make_replicate(1)
  ok(same_counts(r1, r1b, arms_all) && identical(r1$world, r1b$world), "R1 replicate 1 twice from its stream: identical")
  r1s <- make_replicate(1, arms = c("anchor", "E", "Cminus", "O"))
  ok(same_counts(r1, r1s, c("anchor", "E", "Cminus", "O")), "R2 without P and E′: anchor, E, C−, O identical")
  L_ref <- Matrix::colSums(R0$ref)
  t2 <- proc.time()[["elapsed"]]
  RP <- reference_p_from_world(r1$world, R0$ref_meta, L_ref, colnames(R0$ref))
  t_ref <- proc.time()[["elapsed"]] - t2
  RP2 <- reference_p_from_world(r1$world, R0$ref_meta, L_ref, colnames(R0$ref))
  ratio <- Matrix::colSums(RP$ref) / L_ref
  ok(identical(RP$ref, RP2$ref) && ncol(RP$ref) == ncol(R0$ref) && identical(RP$ref_meta, R0$ref_meta) &&
     identical(rownames(RP$ref), rownames(r1$world$mu)) && abs(mean(ratio) - 1) < 0.01,
     sprintf("R3 P reference rebuilt twice identical, %d cells, totals/L mean %.4f (%.1f s)", ncol(RP$ref), mean(ratio), t_ref))
  r2 <- make_replicate(2)
  k10 <- "c10_d4900"
  ok(!identical(r1$world$mu, r2$world$mu) &&
     !identical(r1$conditions[[k10]]$counts$Cminus, r2$conditions[[k10]]$counts$Cminus) &&
     !identical(r1$conditions[[k10]]$truth, r2$conditions[[k10]]$truth), "R4 replicates 1 and 2 differ (μ^P, C−, compositions)")
  x <- r1$conditions[[k10]]$counts
  lcp <- function(m) { m <- as.matrix(m); log1p(sweep(m, 2, colSums(m), "/") * 1e4) }
  la <- lcp(x$anchor)
  k1 <- vapply(c("O", "E", "Eprime", "Cminus"), function(d) { lx <- lcp(x[[d]])
    mean(vapply(seq_len(n_spots), function(s) cor(lx[, s], la[, s]), 1)) }, 1)
  cat(sprintf("K1 mean per-spot Pearson with the anchor (10 cells, 4,900 UMI): O %.4f, E %.4f, E′ %.4f, C− %.4f -> expected O > E > C− and O > E′ > C−: %s\n",
              k1["O"], k1["E"], k1["Eprime"], k1["Cminus"],
              if (k1["O"] > k1["E"] && k1["E"] > k1["Cminus"] && k1["O"] > k1["Eprime"] && k1["Eprime"] > k1["Cminus"])
                "as expected" else "AGAINST the expectation"))
  rs <- r1$conditions[[k10]]$reuse
  cat("Reuse per held-out cell, replicate 1, 10 cells per spot (observed vs expected):\n")
  print(transform(rs[, c("type", "n_held", "drawn", "reuse_obs", "reuse_exp")], reuse_obs = round(reuse_obs, 2), reuse_exp = round(reuse_exp, 2)), row.names = FALSE)
  cat("Spots at / below the target, replicate 1:\n")
  print(do.call(rbind, lapply(r1$conditions, `[[`, "depth"))[, c("n_cells", "depth", "dataset", "n_at_target", "n_below")], row.names = FALSE)
  if (any(res == "FAIL")) stop(sum(res == "FAIL"), " check(s) failed")
  cat("All", length(res), "checks PASS\n")
  quit(save = "no")
}

# ---------------------------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[%s, %s] %d replicates x 24 datasets, %d workers\nOutput: %s\n", tissue, ref_mode, n_rep, workers, out_dir))
log_line(sprintf("run start, %s %s, %d replicates, %d workers", tissue, ref_mode, n_rep, workers))
t0 <- proc.time()[["elapsed"]]
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")   # one core per worker (ADEMP §6); workers inherit it
cl <- parallel::makeCluster(workers, setup_strategy = "sequential")   # parallel start-up drops workers on Windows
lib <- .libPaths()
invisible(parallel::clusterCall(cl, function(lib) { .libPaths(lib); RNGkind("L'Ecuyer-CMRG"); NULL }, lib))
parallel::clusterExport(cl, c("R0", "calib", "ts", "types", "n_train", "alpha_by_n", "mu_fix", "mean_real", "held", "held_type",
  "held_L", "genes", "rep_seed", "substream", "use_stream", "make_world_p", "make_replicate", "save_atomic", "log_line",
  "run_job", "out_dir", "log_file", "tissue", "ref_mode", "study_seed", "n_spots", "levels_cells", "levels_depth",
  "slot", "arms_all"))
invisible(parallel::clusterEvalQ(cl, { suppressPackageStartupMessages(library(Matrix)); for (f in c("R/simulate_spots.R", "R/arm_p.R")) source(f); NULL }))
done <- parallel::parLapplyLB(cl, seq_len(n_rep), function(i)
  tryCatch(run_job(i), error = function(e) list(rep = i, error = conditionMessage(e))))
parallel::stopCluster(cl)
errs <- Filter(function(z) !is.null(z$error), done)
if (length(errs)) { for (z in errs) cat("ERROR rep", z$rep, z$error, "\n"); stop(length(errs), " replicates failed a check") }
new <- Filter(function(z) !z$skipped, done)
cat(sprintf("Done in %.1f min; %d replicates written, %d already on disk\n",
            (proc.time()[["elapsed"]] - t0) / 60, length(new), length(done) - length(new)))
if (length(new)) {
  sum_file <- function(x) file.path(dirname(out_dir), sprintf("%s_sim_%s.csv", tissue, x))
  app <- function(df, f) write.table(df, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
  app(do.call(rbind, lapply(new, `[[`, "depth")), sum_file("depth"))
  app(do.call(rbind, lapply(new, `[[`, "reuse")), sum_file("reuse"))
  app(data.frame(tissue = tissue, rep = vapply(new, `[[`, 1, "rep"), seconds = vapply(new, `[[`, 1, "seconds"),
                 mb = vapply(new, `[[`, 1, "mb"), workers = workers), sum_file("jobs"))
  cat(sprintf("Per replicate: %.1f min, %.0f MB (mean over the new ones)\n",
              mean(vapply(new, `[[`, 1, "seconds")) / 60, mean(vapply(new, `[[`, 1, "mb"))))
  dp <- do.call(rbind, lapply(new, `[[`, "depth"))
  cat("Spots at the target (thinned) and below it, summed over the new replicates:\n")
  print(aggregate(cbind(n_at_target, n_below) ~ dataset + depth + n_cells, dp, sum))
}
