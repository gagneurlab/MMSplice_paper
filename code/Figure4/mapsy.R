library(data.table)
library(ggplot2)
library(cowplot)

projectDIR <- ""

## huber regression
format_cor <- function(clt){
  bquote(italic(R) == .(format.pval(clt)))
}

vivo <- fread(paste0(projectDIR, "Vivo_Pred_huber.txt"))
setnames(vivo, c("Measured", "Prediction"))
vivo_cor <- vivo[, cor(Prediction, Measured)] %>% round(2) %>% format_cor

vivo_plot <- ggplot(vivo, aes(Prediction, Measured)) +
  geom_hex(bins = 70) + 
  scale_fill_distiller(name = "count", trans = "log10", palette = "Spectral") +
  geom_abline(slope = 1, intercept = 0, linetype='dotted') +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  annotate('text', x = -4, y=5, label=c(vivo_cor), parse=T) +
  ggtitle("MaPSy in vivo data") + 
  labs(x='Predicted log2 M/W', y='Measured log2 M/W')
# ggsave(vivo_plot, filename = "MaPSy_vivo.pdf", width = 5, height = 4)
vivo_plot

vitro <- fread(paste0(projectDIR, "Vitro_Pred_huber.txt"))
setnames(vitro, c("Measured", "Prediction"))
vitro_cor <- vitro[, cor(Prediction, Measured)] %>% round(2) %>% format_cor

vitro_plot <- ggplot(vitro, aes(Prediction, Measured)) +
  geom_hex(bins = 70) + 
  scale_fill_distiller(name = "count", trans = "log10", palette = "Spectral") +
  geom_abline(slope = 1, intercept = 0, linetype='dotted') +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  annotate('text', x = -1.5, y=4, label=c(vitro_cor), parse=T) +
  ggtitle("MaPSy in vitro data") +
  labs(x='Predicted log2 M/W', y='Measured log2 M/W')
# ggsave(vitro_plot, filename = "MaPSy_vitro.pdf", width = 5, height = 4)
vitro_plot

## Merge
merged <- rbindlist(list(vivo=vivo, vitro=vitro), idcol="assay")
merged$assay <- list(merged$assay)
merged[, label := ifelse(assay=='vivo', "in vivo", 'in vitro')]

mapsy <- ggplot(merged, aes(Prediction, Measured)) +
  geom_hex(bins = 70) + facet_wrap(~label, scale='free') +
  scale_fill_distiller(name = "count", trans = "log10", palette = "Spectral") +
  geom_abline(slope = 1, intercept = 0, linetype='dotted') +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
  annotate('text', x = -2, y=4, label=c(vitro_cor, vivo_cor), parse=T) +
  labs(x='Predicted allelic ratio', y='Measured allelic ratio')

#mapsy <- plot_grid(vivo_plot, vitro_plot, labels = c("A", "B"))
# ggsave(mapsy, filename = paste0(projectDIR,"plots/MaPSy.pdf"), width = 8, height = 3.5)


