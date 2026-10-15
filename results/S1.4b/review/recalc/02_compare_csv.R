.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

OUTD <- "results/S1.4b/review/recalc"
V <- read.csv("results/S1.4b/S1.4b_verdicts.csv", stringsAsFactors = FALSE)
AR <- read.csv(file.path(OUTD, "recalc_archetype_B.csv"), stringsAsFactors = FALSE)
B4 <- read.csv(file.path(OUTD, "recalc_B4_B8_B9.csv"), stringsAsFactors = FALSE)
vb <- V[V$section == "B" & V$id %in% c("B-1", "B-2", "B-3", "B-7"), ]
vb$variant <- ifelse(grepl("^V-c", vb$metric), "V-c", "V-a")
vb$folds <- sub(".*\\(fold ([PWF]{5})\\).*", "\\1", vb$metric)
mm <- merge(vb, AR, by.x = c("archetype", "variant", "id"), by.y = c("archetype", "variant", "id"), suffixes = c(".csv", ".rec"))
mm$dist_match <- ifelse(mm$id == "B-7", NA, mm$value == mm$dists)
mm$ok <- mm$folds.csv == mm$folds.rec & mm$outcome.csv == mm$outcome.rec
cat("righe B confrontate:", nrow(mm), " fold+esito identici:", sum(mm$ok), " distanze (2 cifre) identiche:", sum(mm$dist_match, na.rm=TRUE), "/", sum(!is.na(mm$dist_match)), "\n")
print(mm[!mm$ok | (!is.na(mm$dist_match) & !mm$dist_match), c("archetype","variant","id","value","dists","folds.csv","folds.rec","outcome.csv","outcome.rec")])
# B-7 values sim | real
b7 <- mm[mm$id == "B-7", ]
b7$sim_csv <- as.numeric(sub(" \\|.*", "", b7$value)); b7$real_csv <- as.numeric(sub(".*\\| ", "", b7$value))
cat("B-7 max |sim diff|:", max(abs(b7$sim_csv - round(b7$b7_sim_mean, 3))), " max |real diff|:", max(abs(b7$real_csv - round(b7$b7_real_mean, 3))), "\n")
# B-4
v4 <- V[V$id == "B-4", ]; v4$variant <- ifelse(grepl("^V-c", v4$metric), "V-c", "V-a")
v4$b4_csv <- as.numeric(sub(" .*", "", v4$value)); v4$loss_csv <- as.numeric(sub(".*\\((.*)\\)", "\\1", v4$value))
m4 <- merge(v4, B4, by = c("archetype", "variant"))
cat("B-4 max |diff| b4:", max(abs(m4$b4_csv - round(m4$b4, 4))), " loss:", max(abs(m4$loss_csv - round(m4$loss, 4))), " esiti identici:", all(m4$outcome.x == m4$outcome.y), "\n")
# folds.csv del pipeline vs mio
FF <- read.csv("results/S1.4b/S1.4b_folds.csv", stringsAsFactors = FALSE); cat("colonne folds.csv:", paste(names(FF), collapse=","), "\n"); print(head(FF, 3))
