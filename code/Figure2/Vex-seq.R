#' ---
#' title: "Benchmark on Vex-seq test data"
#' author: "Jun Cheng"
#' ---

source('code/config.R')

DATADIR <- "data/vexseq/"
FIGDIR <- "data/"

##' ## Read all predictions
MMSplice <- fread(paste0(DATADIR, "MMSplice_test_pred.csv"))
MMSplice <- MMSplice[, .(ID, reference, variant, side, HepG2_ref_psi, HepG2_delta_psi, MMSplice_Delta_PSI_Pred)]
setnames(MMSplice, 'MMSplice_Delta_PSI_Pred', 'MMSplice')
MMSplice[, HepG2_delta_psi := HepG2_delta_psi / 100]
MMSplice[, isSNP := (nchar(reference)==1) & (nchar(variant)==1)]

HAL <- fread(paste0(DATADIR, "HAL_test_pred.tsv"))
HAL <- HAL[, .(ID, DPSI_PRED)]
setnames(HAL, 'DPSI_PRED', 'HAL')
HAL$HAL[abs(HAL$HAL)<1e-10] = 0

SPANR <- fread(paste0(DATADIR, "Spidex_Test_Pred.csv"))
SPANR <- SPANR[!duplicated(SPANR$ID)]
SPANR <- SPANR[, .(ID, dpsi_max_tissue)]
setnames(SPANR, 'dpsi_max_tissue', 'SPANR')
sum(is.na(SPANR$SPANR))
SPANR[, SPANR := SPANR/100]

MaxEnt <- fread(paste0(DATADIR, "MaxEnt_Test_Pred.csv"))
MaxEnt <- MaxEnt[, .(ID, DONOR_DIFF, ACCEPTOR_DIFF)]
setnames(MaxEnt, c('DONOR_DIFF', 'ACCEPTOR_DIFF'), c('MaxEnt5', 'MaxEnt3'))

merged <- merge(MMSplice, HAL, by='ID')
merged <- merge(merged, SPANR, by='ID', all.x=TRUE)
merged <- merge(merged, MaxEnt, by='ID', all.x=TRUE)
merged[, Donor := MaxEnt5 != 0]
merged[, Acceptor := MaxEnt3 != 0]
merged[, SPANR_NA := is.na(SPANR)]
merged[is.na(merged)] <- 0

merged <- melt(merged, 
               measure.vars = c('MMSplice', 'HAL', 'SPANR', 'MaxEnt5', 'MaxEnt3'), 
               value.name = 'Prediction', 
               variable.name = 'Method')

##' ## Calculate correlation, with and without out un-scored variants
##' ### All variants
# format_cor_all <- function(clt){
#   bquote(italic(R) == .(paste(format.pval(clt), "All variants")))
# }
# 
# correlation <- merged[, metric(HepG2_delta_psi, Prediction), 
#                       by='Method']$V1 %>% round(2)
# correlation <- sapply(correlation, format_cor_all)
# 
# ## Bootstrap to calculate p-value
# correlation_bp <- merged[, replicate(999, .SD[sample(nrow(.SD), replace = T)][, cor(Prediction, HepG2_delta_psi)]), by='Method']
# (1 + sum(correlation_bp[Method=='HAL', V1] > correlation_bp[Method=='MMSplice', V1])) / 1000 # 0.001
# (1 + sum(correlation_bp[Method=='SPANR', V1] > correlation_bp[Method=='MMSplice', V1])) / 1000 # 0.001

##' ## only capable points
MMSplice_cor <- merged[Method=='MMSplice'][, metric(HepG2_delta_psi, Prediction)]
hal_cor <- merged[side==''][(isSNP)][Method=='HAL'][, metric(HepG2_delta_psi, Prediction)]
spanr_cor <- merged[!(SPANR_NA)][Method=='SPANR'][, metric(HepG2_delta_psi, Prediction)]
#maxent5_cor <- merged[(Donor)][Method=='MaxEnt5'][, cor(Prediction, HepG2_delta_psi)]
correlation_p <- c(MMSplice_cor, hal_cor, spanr_cor) %>% round(2)
format_cor_valid <- function(clt){
  bquote(italic(R) == .(paste(format.pval(clt), "Scored variants")))
}
correlation_p <- sapply(correlation_p, format_cor_valid)

##' ## Plot
paper_theme <- theme_classic()

p_All <- ggplot(merged[!Method %in% c("MaxEnt5", "MaxEnt3")], aes(Prediction, HepG2_delta_psi)) +
  geom_hex(bins = 70) + 
  geom_abline(slope = 1) +
  scale_fill_distiller(name = "count", trans = "log10", palette = "Spectral") +
  facet_wrap(~Method, scales = 'free', ncol = 2)+paper_theme+
  theme(legend.justification = c(0, 1), legend.position = c(0.5, 0.5)) +
  annotate('text', x = -0.5, y=0.55, label=correlation[1:3], parse=T, size=3) +
  annotate('text', x = -0.45, y=0.45, label=correlation_p, parse=T, size=3) +
  labs(x = expression(Delta*Psi*" predicted"),
       y = expression(Delta*Psi*" measured"))
p_All
# ggsave(paste0(FIGDIR, 'CAGI_All.pdf'), p_All, width=6, height=5)


##' ## Exonic Variants
exon_dt <- merged[side==''][!Method %in% c("MaxEnt5", "MaxEnt3")]
cor_exon <- exon_dt[, cor(Prediction, HepG2_delta_psi), 
                    by='Method']$V1 %>% round(2)
format_cor_exon <- function(clt){
  bquote(italic(R) == .(format.pval(clt)))
}
cor_exon <- sapply(cor_exon, format_cor_exon)

p_Exon <- ggplot(exon_dt, aes(Prediction, HepG2_delta_psi)) +
  geom_hex(bins = 70) + 
  paper_theme +
  theme(legend.justification = c(0, 1), legend.position = c(0.5, 0.5)) +
  scale_fill_distiller(name = "count", trans = "log10", palette = "Spectral") +
  facet_wrap(~Method, scales = 'free', ncol=2) +
  annotate('text', x = -0.52, y=0.40, label=cor_exon, parse=T, size=2) +
  labs(x = expression(Delta*Psi*" predicted"),
       y = expression(Delta*Psi*" measured"))
# ggsave(paste0(FIGDIR, 'CAGI_Exon.pdf'), p_Exon, width=7, height=6)
p_Exon

##' ## Donor Variants
donor_dt <- merged[(Donor)][Method!="MaxEnt3"]
cor_donor <- donor_dt[, cor(Prediction, HepG2_delta_psi), 
                      by='Method']$V1 %>% round(2)
cor_donor <- sapply(cor_donor, format_cor_exon)

p_Donor <- ggplot(donor_dt, aes(Prediction, HepG2_delta_psi)) +
  geom_point(size=1) +
  facet_wrap(~Method, scales = 'free_x') +
  paper_theme +
  annotate('text', x = -0.4, y=0.40, label=cor_donor, parse=T, size=2) +
  labs(x = expression(Delta*Psi*" predicted"),
       y = expression(Delta*Psi*" measured"))
p_Donor
# ggsave(paste0(FIGDIR, 'CAGI_Donor.pdf'), p_Donor, width=5, height=5)

##' ## Acceptor Variants
acceptor_dt <- merged[(Acceptor)][Method!="MaxEnt5"]
cor_acceptor <- acceptor_dt[, cor(Prediction, HepG2_delta_psi), 
                            by='Method']$V1 %>% round(2)
cor_acceptor <- sapply(cor_acceptor, format_cor_exon)

p_Acceptor <- ggplot(merged[(Acceptor)][Method!="MaxEnt5"], aes(Prediction, HepG2_delta_psi)) +
  geom_point(size=1) +
  facet_wrap(~Method, scales = 'free_x') +
  paper_theme +
  annotate('text', x = -0.4, y=0.40, label=cor_acceptor, parse=T, size=2) +
  labs(x = expression(Delta*Psi*" predicted"),
       y = expression(Delta*Psi*" measured"))
p_Acceptor
# ggsave(paste0(FIGDIR, 'CAGI_Acceptor.pdf'), p_Acceptor, width=5, height=5)


##' ## Plot MMSplice prediction with simpy sum of 5 modules
MMSplice <- fread(paste0(DATADIR, "MMSplice_test_pred.csv"))
MMSplice <- MMSplice[, .(ID, HepG2_delta_psi, MMSplice_Delta_PSI_Pred_Sum, MMSplice_Delta_PSI_Pred_Sum3)]
setnames(MMSplice, c("MMSplice_Delta_PSI_Pred_Sum", "MMSplice_Delta_PSI_Pred_Sum3"), c("Sum Prediction Exon 5' Module", "Sum Prediction Exon 3' Module"))
MMSplice <- melt(MMSplice, 
                 measure.vars = c("Sum Prediction Exon 5' Module", "Sum Prediction Exon 3' Module"),
                 variable.name = 'Model',
                 value.name = "Prediction")

cor_sum <- MMSplice[, cor(Prediction, HepG2_delta_psi), by='Model']$V1 %>% round(2)
cor_sum <- sapply(cor_sum, format_cor_exon)
MMSplice[, HepG2_delta_psi := HepG2_delta_psi/100]

p_All_Sum <- ggplot(MMSplice, aes(Prediction, HepG2_delta_psi)) +
  geom_hex(bins = 70) + 
  scale_fill_distiller(name = "count", trans = "log10", palette = "Spectral") +
  facet_wrap(~Model, scales = 'free', ncol = 2)+
  paper_theme +
  annotate('text', x = -0.6, y=0.40, label=cor_sum, parse=T, size=3) +
  labs(x = expression(Delta*Psi*" predicted"),
       y = expression(Delta*Psi*" measured"))
# ggsave(paste0(FIGDIR, 'CAGI_All_Sum.pdf'), p_All_Sum, width = 8, height = 4)
p_All_Sum

##' ## Plot Exon prediction with two exon modules on training
ExonPrime3 <- fread(paste0(DATADIR, 'Exon_Pred_Prime3.txt'))
setnames(ExonPrime3, c("Prediction", "Measurement"))
ExonPrime3[, Module := "Exon 5' Module"]
ExonPrime5 <- fread(paste0(DATADIR, 'Exon_Pred_Prime5.txt'))
setnames(ExonPrime5, c("Prediction", "Measurement"))
ExonPrime5[, Module := "Exon 3' Module"]

ExonPred <- rbind(ExonPrime5, ExonPrime3)

cor_exons <- ExonPred[, cor(Prediction, Measurement), by='Module']$V1 %>% round(2)
cor_exons <- sapply(cor_exons, format_cor_exon)
#ExonPred[, Measurement := Measurement/100]
#ExonPred[, Prediction := Prediction/100]

p_exon_comp <- ggplot(ExonPred, aes(Prediction, Measurement)) +
  geom_hex(bins = 70) + 
  scale_fill_distiller(name = "count", trans = "log10", palette = "Spectral") +
  paper_theme +
  facet_wrap(~Module, scale='free') +
  annotate('text', x = -0.07, y=0.40, label=cor_exons, parse=T, size=3) +
  labs(x = expression(Delta*Psi*" predicted"),
       y = expression(Delta*Psi*" measured"))
# ggsave(paste0(FIGDIR, 'CAGI_Exon_Compare.pdf'), width = 8, height = 4)
p_exon_comp

