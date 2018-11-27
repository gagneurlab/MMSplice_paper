library(ggplot2)
library(data.table)
library(cowplot)
library(magrittr)

projectDIR <- "data/gtex/"

paper_theme <- theme(plot.title = element_text(hjust = 0.5)) + theme_cowplot() + theme(legend.position="none")
dot_size <- 0.8
cbPalette <- c("#0072B2", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")

##' ## SD variants
MMSplice_homo <- fread(paste0(projectDIR, "SD_HOMO.csv"))
MMSplice_homo[, dPSI_Measured := HOMO_MEAN - WT_MEAN]
MMSplice_hetero <- fread(paste0(projectDIR, "SD_HETERO.csv"))
MMSplice_hetero[, dPSI_Measured := HETERO_MEAN - WT_MEAN]
MMSplice <- rbindlist(list(homo=MMSplice_homo[, .(HOMO_dPSI, dPSI_Measured)], 
                           hetero=MMSplice_hetero[,.(HETERO_dPSI, dPSI_Measured)]), idcol='Genotype')
setnames(MMSplice, "HOMO_dPSI", "dPSI_Pred")

COSSMO_homo <- fread(paste0(projectDIR, "COSSMO_SD_HOMO.csv"))
COSSMO_homo[, dPSI_Measured := HOMO_MEAN - WT_MEAN]
COSSMO_hetero <- fread(paste0(projectDIR, "COSSMO_SD_HETERO.csv"))
COSSMO_hetero[, dPSI_Measured := HETERO_MEAN - WT_MEAN]
COSSMO <- rbindlist(list(homo=COSSMO_homo[, .(DPSI_PRED_HOMO, dPSI_Measured)],
                         hetero=COSSMO_hetero[, .(DPSI_PRED_HETERO, dPSI_Measured)]), idcol='Genotype')
setnames(COSSMO, "DPSI_PRED_HOMO", "dPSI_Pred")

HAL_homo <- fread(paste0(projectDIR, "HAL_SD_HOMO.csv"))
HAL_homo[, dPSI_Measured := HOMO_MEAN - WT_MEAN]
HAL_hetero <- fread(paste0(projectDIR, "HAL_SD_HETERO.csv"))
HAL_hetero[, dPSI_Measured := HETERO_MEAN - WT_MEAN]
HAL <- rbindlist(list(homo=HAL_homo[, .(HOMO_dPSI, dPSI_Measured)],
                      hetero=HAL_hetero[, .(HETERO_dPSI, dPSI_Measured)]), idcol='Genotype')
MaxEnt <- rbindlist(list(homo=HAL_homo[, .(dMaxEnt_homo, dPSI_Measured)],
                         hetero=HAL_hetero[, .(dMaxEnt_hetero, dPSI_Measured)]), idcol='Genotype') 
setnames(HAL, "HOMO_dPSI", "dPSI_Pred")
setnames(MaxEnt, "dMaxEnt_homo", "dPSI_Pred")

merged <- rbindlist(list(MMSplice=MMSplice, COSSMO=COSSMO, HAL=HAL, MaxEntScan=MaxEnt), idcol = "Method")
set.seed(222)
correlation_bp <- merged[, replicate(100, .SD[sample(nrow(.SD), replace = T)][, cor(dPSI_Pred, dPSI_Measured)]), by='Method']
t.test(correlation_bp[Method=='MMSplice', V1], correlation_bp[Method=='COSSMO', V1])
t.test(correlation_bp[Method=='MMSplice', V1], correlation_bp[Method=='HAL', V1])
t.test(correlation_bp[Method=='MMSplice', V1], correlation_bp[Method=='MaxEntScan', V1])
correlation_bp[, mean(V1), by='Method']

ggplot(merged, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point() +
  facet_wrap(~Method, scale='free_x') +
  theme_classic() +
  labs(x = expression(Delta*"PSI Predicted"),
       y = expression(Delta*"PSI Measured"))

format_cor <- function(clt){
  clt <- round(clt, 2)
  bquote(italic(R) == .(format.pval(clt)))
}

##' ### Combine SD Plots
mmsplice_cor <- MMSplice[, cor(dPSI_Pred, dPSI_Measured)] %>% format_cor
mmsplice_d <- ggplot(MMSplice, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point(size=dot_size) +
  paper_theme +
  labs(x = expression(Delta*Psi[3]*" Predicted"),
       y = expression(Delta*Psi[3]*" Measured")) +
  annotate('text', x = -0.6, y=0.8, label=c(mmsplice_cor), parse=T, size=4) +
  #theme(legend.justification = c(0, 1), legend.position = c(-50, 50)) +
  labs(title = 'MMSplice - A5SS') +
  scale_x_continuous(limits = c(-1, 1)) + 
  scale_y_continuous(limits = c(-1, 1)) +
  geom_abline(slope = 1, linetype="dotted") +
  scale_colour_manual(values=cbPalette)

cossmo_cor <- COSSMO[, cor(dPSI_Pred, dPSI_Measured)] %>% format_cor
cossmo_d <- ggplot(COSSMO, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point(size=dot_size) +
  paper_theme +
  labs(x = expression(Delta*Psi[3]*" Predicted"),
       y = expression(Delta*Psi[3]*" Measured")) +
  annotate('text', x = -0.6, y=0.8, label=c(cossmo_cor), parse=T, size=4) +
  labs(title = 'COSSMO - A5SS') +
  scale_x_continuous(limits = c(-1, 1)) + 
  scale_y_continuous(limits = c(-1, 1)) +
  geom_abline(slope = 1, linetype="dotted") +
  scale_colour_manual(values=cbPalette)

hal_cor <- HAL[, cor(dPSI_Pred, dPSI_Measured)] %>% format_cor
hal_d <- ggplot(HAL, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point(size=dot_size) +
  paper_theme +
  labs(x = expression(Delta*Psi[3]*" Predicted"),
       y = expression(Delta*Psi[3]*" Measured")) +
  annotate('text', x = -0.6, y=0.8, label=c(hal_cor), parse=T, size=4) +
  labs(title = 'HAL - A5SS') +
  scale_x_continuous(limits = c(-1, 1)) + 
  scale_y_continuous(limits = c(-1, 1)) +
  geom_abline(slope = 1, linetype="dotted") +
  scale_colour_manual(values=cbPalette)

maxent_cor <- MaxEnt[, cor(dPSI_Pred, dPSI_Measured)] %>% format_cor
maxent_d <- ggplot(MaxEnt, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point(size=dot_size) +
  paper_theme +
  labs(x = expression(Delta*"MaxEnt Score"),
       y = expression(Delta*Psi[3]*" Measured")) +
  annotate('text', x = -9, y=0.8, label=c(maxent_cor), parse=T, size=4) +
  labs(title = 'MaxEnt - A5SS') +
  # scale_x_continuous(limits = c(-1, 1)) + 
  # scale_y_continuous(limits = c(-1, 1)) +
  #geom_abline(slope = 1, linetype="dotted") +
  scale_colour_manual(values=cbPalette)

donor <- plot_grid(mmsplice_d, cossmo_d, hal_d, maxent_d, 
                   labels = c('A','B','C','D'), ncol=2)
# ggsave(paste0(projectDIR, 'Output/fig/donor_GTEx.pdf'), donor, width = 10, height = 7)

##' ## SA variants
MMSplice_homo <- fread(paste0(projectDIR, "SA_HOMO.csv"))
MMSplice_homo[, dPSI_Measured := HOMO_MEAN - WT_MEAN]
MMSplice_hetero <- fread(paste0(projectDIR, "SA_HETERO.csv"))
MMSplice_hetero[, dPSI_Measured := HETERO_MEAN - WT_MEAN]
MMSplice <- rbindlist(list(homo=MMSplice_homo[, .(HOMO_dPSI, dPSI_Measured)], 
                           hetero=MMSplice_hetero[,.(HETERO_dPSI, dPSI_Measured)]), idcol='Genotype')
setnames(MMSplice, "HOMO_dPSI", "dPSI_Pred")

COSSMO_homo <- fread(paste0(projectDIR, "COSSMO_SA_HOMO.csv"))
COSSMO_homo[, dPSI_Measured := HOMO_MEAN - WT_MEAN]
COSSMO_hetero <- fread(paste0(projectDIR, "COSSMO_SA_HETERO.csv"))
COSSMO_hetero[, dPSI_Measured := HETERO_MEAN - WT_MEAN]
COSSMO <- rbindlist(list(homo=COSSMO_homo[, .(DPSI_PRED_HOMO, dPSI_Measured)],
                         hetero=COSSMO_hetero[, .(DPSI_PRED_HETERO, dPSI_Measured)]), idcol='Genotype')
setnames(COSSMO, "DPSI_PRED_HOMO", "dPSI_Pred")

MaxEnt <- rbindlist(list(homo=COSSMO_homo[, .(dMaxEnt_homo, dPSI_Measured)],
                         hetero=COSSMO_hetero[, .(dMaxEnt_hetero, dPSI_Measured)]), idcol='Genotype') 
setnames(MaxEnt, "dMaxEnt_homo", "dPSI_Pred")

merged <- rbindlist(list(MMSplice=MMSplice, COSSMO=COSSMO, MaxEntScan=MaxEnt), idcol = "Method")
correlation_bp <- merged[, replicate(100, .SD[sample(nrow(.SD), replace = T)][, cor(dPSI_Pred, dPSI_Measured)]), by='Method']
t.test(correlation_bp[Method=='MMSplice', V1], correlation_bp[Method=='COSSMO', V1])
t.test(correlation_bp[Method=='MMSplice', V1], correlation_bp[Method=='MaxEntScan', V1])
correlation_bp[, mean(V1), by='Method']

ggplot(merged, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point(size=dot_size) +
  facet_wrap(~Method, scale='free_x') +
  theme_classic() +
  labs(x = expression(Delta*"PSI Predicted"),
       y = expression(Delta*"PSI Measured"))

##' ### Combine SD Plots
mmsplice_cor <- MMSplice[, cor(dPSI_Pred, dPSI_Measured)] %>% format_cor
mmsplice_a <- ggplot(MMSplice, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point(size=dot_size) +
  paper_theme +
  labs(x = expression(Delta*Psi[5]*" Predicted"),
       y = expression(Delta*Psi[5]*" Measured")) +
  annotate('text', x = -0.6, y=0.8, label=c(mmsplice_cor), parse=T, size=4) +
  #theme(legend.justification = c(0, 1), legend.position = c(-50, 50)) +
  labs(title = 'MMSplice - A3SS') +
  scale_x_continuous(limits = c(-1, 1)) + 
  scale_y_continuous(limits = c(-1, 1)) +
  geom_abline(slope = 1, linetype="dotted") +
  scale_colour_manual(values=cbPalette)

cossmo_cor <- COSSMO[, cor(dPSI_Pred, dPSI_Measured)] %>% format_cor
cossmo_a <- ggplot(COSSMO, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point(size=dot_size) +
  paper_theme +
  labs(x = expression(Delta*Psi[5]*" Predicted"),
       y = expression(Delta*Psi[5]*" Measured")) +
  annotate('text', x = -0.6, y=0.8, label=c(cossmo_cor), parse=T, size=4) +
  labs(title = 'COSSMO - A3SS') +
  scale_x_continuous(limits = c(-1, 1)) + 
  scale_y_continuous(limits = c(-1, 1)) +
  geom_abline(slope = 1, linetype="dotted") +
  scale_colour_manual(values=cbPalette)

maxent_cor <- MaxEnt[, cor(dPSI_Pred, dPSI_Measured)] %>% format_cor
maxent_a <- ggplot(MaxEnt, aes(dPSI_Pred, dPSI_Measured, color=Genotype)) +
  geom_point(size=dot_size) +
  paper_theme +
  labs(x = expression(Delta*"MaxEnt Score"),
       y = expression(Delta*Psi[5]*" Measured")) +
  annotate('text', x = -9, y=0.8, label=c(maxent_cor), parse=T, size=4) +
  labs(title = 'MaxEnt - A3SS') +
  # scale_x_continuous(limits = c(-1, 1)) + 
  # scale_y_continuous(limits = c(-1, 1)) +
  #geom_abline(slope = 1, linetype="dotted") +
  theme(legend.position=c(1.1, 0.8)) +
  scale_colour_manual(values=cbPalette)

acceptor <- plot_grid(mmsplice_a, cossmo_a, maxent_a, 
                      labels = c('A','B','C'), ncol=2)

ss <- plot_grid(mmsplice_d, cossmo_d, hal_d, maxent_d, mmsplice_a, cossmo_a, maxent_a,
                labels = c('A','B','C','D','E','F','G'), nrow = 2)
# ggsave(paste0(projectDIR, 'SS_GTEx.pdf'), ss, width = 15, height = 7)

