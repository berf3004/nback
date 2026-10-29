# =============================================================================
# Tâche n-back : cotation à partir du fichier des essais (nback_essais_*.csv)
# R de base seulement, aucun paquet requis.
#
# Produit une ligne par séance (participant x temps de mesure) avec, pour chaque
# niveau (0-back à 3-back), les indices utilisés dans la littérature :
#   - nombres de hits, omissions, fausses alarmes, rejets corrects, non-réponses
#   - proportions de hits (H) et de fausses alarmes (FA), exactitude
#   - d' et c avec correction log-linéaire (Hautus, 1995)
#   - d' avec la correction 1/(2N) de Stanislaw et Todorov (1999),
#     utilisée par Pelegrina et coll. (2015)
#   - Pr = H - FA (Snodgrass et Corwin, 1988)
#   - TR moyen, médian et écart type des hits ; TR moyen des rejets corrects
#   - score d'efficience inverse (IES = TR moyen des hits / H)
#   - fausses alarmes aux leurres et aux non-leurres (si des leurres sont prévus)
#   - d' par bloc et coefficient de fidélité test-retest entre blocs (bloc 1 vs 2)
#   - effets de charge : d'(2-back) - d'(1-back), d'(2-back) - d'(0-back)
#
# Règles appliquées (modifiables plus bas) :
#   - blocs de test seulement (phase == "test")
#   - essais de départ exclus (start_trial == 1)
#   - TR : réponses de moins de RT_MIN ms retirées des moyennes de TR seulement
#     (réponses anticipées) ; les comptes de hits et de FA ne changent pas.
# =============================================================================

FICHIER_ESSAIS <- file.choose()   # ou "chemin/vers/nback_essais_AAAA-MM-JJ_HHMM.csv"
FICHIER_SORTIE <- sub("\\.csv$", "_scores.csv", FICHIER_ESSAIS)
RT_MIN <- 150                     # ms ; mettre 0 pour ne rien retirer

# Lecture : détecte le format (virgule ou point-virgule avec virgule décimale)
premiere <- readLines(FICHIER_ESSAIS, n = 1, encoding = "UTF-8")
if (grepl(";", premiere)) {
  d <- read.csv2(FICHIER_ESSAIS, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
} else {
  d <- read.csv(FICHIER_ESSAIS, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
}

d <- d[d$phase == "test" & d$start_trial == 0, ]
d$cle <- paste(d$session_id)

ll_d <- function(h, nt, fa, nn) qnorm((h + .5) / (nt + 1)) - qnorm((fa + .5) / (nn + 1))
ll_c <- function(h, nt, fa, nn) -(qnorm((h + .5) / (nt + 1)) + qnorm((fa + .5) / (nn + 1))) / 2
st_d <- function(h, nt, fa, nn) {
  H <- h / nt; F <- fa / nn
  H <- ifelse(H <= 0, 1 / (2 * nt), ifelse(H >= 1, 1 - 1 / (2 * nt), H))
  F <- ifelse(F <= 0, 1 / (2 * nn), ifelse(F >= 1, 1 - 1 / (2 * nn), F))
  qnorm(H) - qnorm(F)
}
moy <- function(x) if (length(x)) mean(x) else NA_real_
med <- function(x) if (length(x)) median(x) else NA_real_
et  <- function(x) if (length(x) > 1) sd(x) else NA_real_

coter <- function(a) {
  T <- a[a$is_target == 1, ]; N <- a[a$is_target == 0, ]
  h <- sum(T$sdt_outcome == "hit"); fa <- sum(N$sdt_outcome == "false_alarm")
  nt <- nrow(T); nn <- nrow(N)
  rt_h <- a$rt_ms[a$sdt_outcome == "hit" & !is.na(a$rt_ms) & a$rt_ms >= RT_MIN]
  rt_c <- a$rt_ms[a$sdt_outcome == "correct_rejection" & !is.na(a$rt_ms) & a$rt_ms >= RT_MIN]
  lu <- N[N$lure_type != "none", ]; nl <- N[N$lure_type == "none", ]
  H <- if (nt) h / nt else NA
  out <- c(
    trials = nrow(a), targets = nt, nontargets = nn,
    hits = h, misses = nt - h, false_alarms = fa,
    correct_rejections = sum(N$sdt_outcome == "correct_rejection"),
    no_responses = sum(a$sdt_outcome == "no_response"),
    accuracy = mean(a$correct),
    hit_rate = H, fa_rate = if (nn) fa / nn else NA,
    dprime = if (nt && nn) ll_d(h, nt, fa, nn) else NA,
    criterion_c = if (nt && nn) ll_c(h, nt, fa, nn) else NA,
    dprime_st = if (nt && nn) st_d(h, nt, fa, nn) else NA,
    pr = if (nt && nn) H - fa / nn else NA,
    rt_hit_mean_ms = moy(rt_h), rt_hit_median_ms = med(rt_h), rt_hit_sd_ms = et(rt_h),
    rt_cr_mean_ms = moy(rt_c),
    ies_ms = if (!is.na(H) && H > 0) moy(rt_h) / H else NA,
    lures = nrow(lu),
    lure_fa_rate = if (nrow(lu)) mean(lu$sdt_outcome == "false_alarm") else NA,
    nonlure_fa_rate = if (nrow(nl)) mean(nl$sdt_outcome == "false_alarm") else NA
  )
  for (b in sort(unique(a$block_num))) {
    ab <- a[a$block_num == b, ]
    tb <- sum(ab$is_target == 1); nb <- sum(ab$is_target == 0)
    out[paste0("dprime_block", b)] <- ll_d(sum(ab$sdt_outcome == "hit"), tb,
                                           sum(ab$sdt_outcome == "false_alarm"), nb)
  }
  out
}

id_cols <- c("participant_id", "session_label", "group_label", "age_years", "session_id",
             "session_start", "session_status", "stimulus_set", "response_mode", "yes_side")
id_cols <- intersect(id_cols, names(d))
lignes <- lapply(split(d, d$cle), function(s) {
  r <- s[1, id_cols, drop = FALSE]
  for (L in sort(unique(s$n_level))) {
    sc <- coter(s[s$n_level == L, ])
    r[paste0("n", L, "_", names(sc))] <- as.list(sc)
  }
  r
})
noms <- unique(unlist(lapply(lignes, names)))
res <- do.call(rbind, lapply(lignes, function(r) { r[setdiff(noms, names(r))] <- NA; r[noms] }))

# Effets de charge
if (all(c("n2_dprime", "n1_dprime") %in% names(res))) res$load_dprime_2v1 <- res$n2_dprime - res$n1_dprime
if (all(c("n2_dprime", "n0_dprime") %in% names(res))) res$load_dprime_2v0 <- res$n2_dprime - res$n0_dprime
if (all(c("n2_rt_hit_mean_ms", "n1_rt_hit_mean_ms") %in% names(res))) res$load_rt_2v1 <- res$n2_rt_hit_mean_ms - res$n1_rt_hit_mean_ms

num <- sapply(res, is.numeric)
res[num] <- lapply(res[num], function(x) round(x, 3))
rownames(res) <- NULL
write.csv(res, FICHIER_SORTIE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat("Scores écrits dans :", FICHIER_SORTIE, "\n")

# Fidélité entre blocs (corrélation de Pearson et Spearman-Brown) pour chaque niveau
for (L in 0:3) {
  b1 <- paste0("n", L, "_dprime_block1"); b2 <- paste0("n", L, "_dprime_block2")
  if (all(c(b1, b2) %in% names(res)) && sum(complete.cases(res[, c(b1, b2)])) >= 5) {
    r <- cor(res[[b1]], res[[b2]], use = "complete.obs")
    cat(sprintf("%d-back : r(bloc 1, bloc 2) = %.2f ; Spearman-Brown = %.2f (n = %d)\n",
                L, r, 2 * r / (1 + r), sum(complete.cases(res[, c(b1, b2)]))))
  }
}
