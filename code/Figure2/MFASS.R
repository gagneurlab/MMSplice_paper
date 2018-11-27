#' ---
#' title: "Evaluation Models on MFASS"
#' author: "Jun Cheng"
#' date: "May 3rd, 2018"
#' ---

# ```{r global_options, include=FALSE}
# knitr::opts_knit$set(echo = TRUE, warning = FALSE)
# ```

library(pROC)
library(data.table)
library(ggplot2)
require(pracma)

DATADIR <- "data/mfass/"
FIGDIR <- "data/"

MMSplice = fread(paste0(DATADIR, "MMSplice_pred.txt"), header=T)
setnames(MMSplice, c("label", "score", "Dpsi", "region", "id"))

## PR from Kimberly D. Insigne
spanr <- fread(paste0(DATADIR, 'snv_SPANR_scores_capped.txt'),
               sep = '\t', header = T)
hal <- fread(paste0(DATADIR, 'SNV_HAL_scores.txt'),
             sep = '\t', header = T)

## Implement a Precision-Recall curve
precision_recall <- function(pred, label, threshold){
  # label should be TRUE and FALSE
  tab <- table(label, pred < threshold)
  tp <- tab[2,2]
  fp <- tab[1,2]
  fn <- tab[2,1]
  precision <- tp*100 / (tp + fp)
  recall <- tp*100 / (tp + fn)
  return(c(precision, recall))
}

precision_recall_curve <- function(pred, label, by=0.05){
  ranges <- seq(min(pred)+1e-7, max(pred)-1e-7, by = by)
  pr <- sapply(ranges, function(i) precision_recall(pred, label, i))
  pr <- data.table(t(pr))
  pr
}

##' ## Overvall performance
##' ### Calculate PR for MMSplice
PR_MMSplice <- precision_recall_curve(MMSplice$score, MMSplice$label, by=0.005) #data.table(t(PR_MMSplice))
setnames(PR_MMSplice, c("Precision", "Recall"))
# PR_MMSplice[, threshold := ranges]
PR_MMSplice[, method := "MMSplice"]
PR_MMSplice[, threshold := seq(min(MMSplice$score)+1e-7, max(MMSplice$score)-1e-7, by = 0.005)]

##' ### Calculate PR for HAL
hal[, DPSI_pred := as.numeric(DPSI_pred)]
hal <- hal[strong_lof != 'NA']
PR_HAL <- precision_recall_curve(hal$DPSI_pred, hal$strong_lof, by = 0.005)
setnames(PR_HAL, c("Precision", "Recall"))
# PR_MMSplice[, threshold := ranges]
PR_HAL[, method := "HAL"]
PR_HAL[, threshold := seq(min(hal$DPSI_pred)+1e-7, max(hal$DPSI_pred)-1e-7, by = 0.005)]

##' ### Calculate PR for SPANR
PR_SPANR <- precision_recall_curve(spanr$dpsi_spanr_capped, spanr$strong_lof, by = 0.005)
setnames(PR_SPANR, c("Precision", "Recall"))
PR_SPANR[, method := "SPANR"]
PR_SPANR[, threshold := seq(min(spanr$dpsi_spanr_capped)+1e-7, max(spanr$dpsi_spanr_capped)-1e-7, by = 0.005)]

PR_all <- rbindlist(list(PR_HAL, PR_MMSplice, PR_SPANR))
# PR_all[is.na(PR_all)] = 100

auPR_All <- PR_all[, trapz(Recall, Precision)/10000 , by='method']$V1 %>% round(2)
auPR_All <- paste0("auPR=", auPR_All)
auPR_All <- paste(c("HAL:", "MMSplice:", "SPANR:"), auPR_All)
model <- gsub("MMSplice", auPR_All[2], PR_all$method)
model <- gsub("HAL", auPR_All[1], model)
model <- gsub("SPANR", auPR_All[3], model)
PR_all$method <- model

##' ## Calculate precision recall with bootstrap
pr_boot <- function(dt, y_hat, y){
  dt <- dt[sample(nrow(dt), replace = TRUE)]
  pr <- precision_recall_curve(dt[, get(y_hat)], dt[, get(y)])
  setnames(pr, c("Precision", "Recall"))
  pr[, trapz(Recall, Precision)/10000]
}

N <- 999
mmsplice_prs <- replicate(N, pr_boot(MMSplice, "score", "label"))
hal_prs <- replicate(N, pr_boot(hal, "DPSI_pred", "strong_lof"))
spanr_prs <- replicate(N, pr_boot(spanr, "dpsi_spanr_capped", "strong_lof"))

saveRDS(mmsplice_prs, paste0(DATADIR, "mmsplice_prs.rds"))
saveRDS(hal_prs, paste0(DATADIR, "hal_prs.rds"))
saveRDS(spanr_prs, paste0(DATADIR, "spanr_prs.rds"))

(1 + sum(hal_prs > mmsplice_prs)) / 1000  #0.001
(1 + sum(spanr_prs > mmsplice_prs)) / 1000  #0.001

##' ## Compare with SPANR, all regions
library(cowplot)
cbPalette <- c("#0072B2", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")
All <- ggplot(PR_all, aes(Recall, Precision, color=method)) +
  geom_line() + scale_colour_manual(values=cbPalette) +
  #scale_y_log10() + 
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("All Variants") +
  theme(legend.justification = c(0, 1), legend.position = c(0.4, 0.9), axis.title = element_text(size=15)) 
All

##' ## Exon performance
##' ### Calculate PR for MMSplice
MMSplice_Exon <- MMSplice[region == 'exon']
PR_MMSplice_Exon <- precision_recall_curve(MMSplice_Exon$score, MMSplice_Exon$label, by=0.005) #data.table(t(PR_MMSplice))
setnames(PR_MMSplice_Exon, c("Precision", "Recall"))
# PR_MMSplice[, threshold := ranges]
PR_MMSplice_Exon[, method := "MMSplice"]
PR_MMSplice_Exon[, threshold := seq(min(MMSplice_Exon$score)+1e-7, max(MMSplice_Exon$score)-1e-7, by = 0.005)]

##' ### Calculate PR for SPANR
SPANR_Exon <- spanr[location=='exonic']
PR_SPANR_Exon <- precision_recall_curve(SPANR_Exon$dpsi_spanr_capped, SPANR_Exon$strong_lof, by = 0.005)
setnames(PR_SPANR_Exon, c("Precision", "Recall"))
PR_SPANR_Exon[, method := "SPANR"]
PR_SPANR_Exon[, threshold := seq(min(SPANR_Exon$dpsi_spanr_capped)+1e-7, max(SPANR_Exon$dpsi_spanr_capped)-1e-7, by = 0.005)]
PR_exon <- rbindlist(list(PR_HAL, PR_MMSplice_Exon, PR_SPANR_Exon))

auPR_Exon <- PR_exon[, trapz(Recall, Precision)/10000, by='method']$V1 %>% round(2) 
auPR_Exon <- paste0("auPR=", auPR_Exon)
auPR_Exon <- paste(c("HAL:", "MMSplice:", "SPANR:"), auPR_Exon)
model <- gsub("MMSplice", auPR_Exon[2], PR_exon$method)
model <- gsub("HAL", auPR_Exon[1], model)
model <- gsub("SPANR", auPR_Exon[3], model)
PR_exon$method <- model

N <- 100
mmsplice_exon_prs <- replicate(N, pr_boot(MMSplice_Exon, "score", "label"))
#hal_exon_prs <- replicate(N, pr_boot(hal, "DPSI_pred", "strong_lof"))
spanr_exon_prs <- replicate(N, pr_boot(SPANR_Exon, "dpsi_spanr_capped", "strong_lof"))

saveRDS(mmsplice_exon_prs, paste0(DATADIR, "mmsplice_exon_prs.rds"))
#saveRDS(hal_prs, paste0(DATADIR, "hal_prs.rds"))
saveRDS(spanr_exon_prs, paste0(DATADIR, "spanr_exon_prs.rds"))

t.test(mmsplice_exon_prs[1:100], spanr_exon_prs[1:100]) #  < 2.2e-16
t.test(mmsplice_exon_prs[1:100], hal_prs[1:100]) # 2.946e-06

##' ## Compare with SPANR, HAL, Exon
Exon <- ggplot(PR_exon, aes(Recall, Precision, color=method)) +
  geom_line() + scale_colour_manual(values=cbPalette) +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Exonic Variants") +
  theme(legend.justification = c(0, 1), legend.position = c(0.5, 0.9))
ggsave(Exon, filename = 'MFASS_Exon.pdf', width = 5, height = 4)
Exon

##' ## Intron performance
##' ### Calculate PR for MMSplice
MMSplice_Intron <- MMSplice[region %in% c('upstr_intron', 'downstr_intron')]
PR_MMSplice_Intron <- precision_recall_curve(MMSplice_Intron$score, MMSplice_Intron$label, by=0.005) #data.table(t(PR_MMSplice))
setnames(PR_MMSplice_Intron, c("Precision", "Recall"))
# PR_MMSplice[, threshold := ranges]
PR_MMSplice_Intron[, method := "MMSplice"]
PR_MMSplice_Intron[, threshold := seq(min(MMSplice_Intron$score)+1e-7, max(MMSplice_Intron$score)-1e-7, by = 0.005)]

##' ## Calculate PR for SPANR
SPANR_Intron <- spanr[location=='intronic']
PR_SPANR_Intron <- precision_recall_curve(SPANR_Intron$dpsi_spanr_capped, SPANR_Intron$strong_lof, by = 0.005)
setnames(PR_SPANR_Intron, c("Precision", "Recall"))
PR_SPANR_Intron[, method := "SPANR"]
PR_SPANR_Intron[, threshold := seq(min(SPANR_Intron$dpsi_spanr_capped)+1e-7, max(SPANR_Intron$dpsi_spanr_capped)-1e-7, by = 0.005)]

PR_Intron <- rbindlist(list(PR_MMSplice_Intron, PR_SPANR_Intron))
auPR_Intron <- PR_Intron[, trapz(Recall, Precision)/10000, by='method']$V1 %>% round(2) 
auPR_Intron <- paste0("auPR=", auPR_Intron)
auPR_Intron <- paste(c("MMSplice:", "SPANR:"), auPR_Intron)
model <- gsub("MMSplice", auPR_Intron[1], PR_Intron$method)
model <- gsub("SPANR", auPR_Intron[2], model)
PR_Intron$method <- model

N <- 100
mmsplice_intron_prs <- replicate(N, pr_boot(MMSplice_Intron, "score", "label"))
spanr_intron_prs <- replicate(N, pr_boot(SPANR_Intron, "dpsi_spanr_capped", "strong_lof"))

saveRDS(mmsplice_intron_prs, paste0(DATADIR, "mmsplice_intron_prs.rds"))
saveRDS(spanr_intron_prs, paste0(DATADIR, "spanr_intron_prs.rds"))

t.test(mmsplice_intron_prs, spanr_intron_prs) #  < 2.2e-16

##' ## Compare with SPANR, HAL, Intron
Intron <- ggplot(PR_Intron, aes(Recall, Precision, color=method)) +
  geom_line() + scale_colour_manual(values=cbPalette) +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Intronic Variants") +
  theme(legend.justification = c(0, 1), legend.position = c(0.5, 0.9))
ggsave(Intron, filename = 'MFASS_Intron.pdf', width = 5, height = 4)
Intron

library(gridExtra)
ggsave(plot_grid(Exon, Intron, labels = c("A", "B")), 
       filename = 'MFASS_EI.pdf', 
       width = 10, height = 4)
