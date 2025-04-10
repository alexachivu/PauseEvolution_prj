load("/local/workdir/bab364/PauseEvolution/NucleosomePositioningWorkspace.Rdata.RData")

# load proseq from PauseEvoPlots
HUM_PE
DROS_PE
SCER_PE
SPOM_PE
# load procap
SPOM_cap <- import_bigWig ("/fs/cbsubscb17/storage/data/S.pombe/procap/GSM1974988_3870_7157_12385_C53ARACXX_pombe972h-ProCap_ATCACG_R1Normed_plus.bw", "/fs/cbsubscb17/storage/data/S.pombe/procap/GSM1974988_3870_7157_12385_C53ARACXX_pombe972h-ProCap_ATCACG_R1Normed_minus.bw")
seqlevelsStyle(SPOM_cap) <- seqlevelsStyle(SPOM_PE)
seqlevels(SPOM_cap) <- seqlevels(SPOM_PE)
SCER_cap <- import_bigWig("/fs/cbsubscb17/storage/data/sacser/procap/GSM1974987_3870_7157_12387_C53ARACXX_cerevisiaeW303-aProCap_TTAGGC_R1Normed_plus.bw", "/fs/cbsubscb17/storage/data/sacser/procap/GSM1974987_3870_7157_12387_C53ARACXX_cerevisiaeW303-aProCap_TTAGGC_R1Normed_minus.bw")
seqlevelsStyle(SCER_cap) <- seqlevelsStyle(SCER_5p)
DROS_cap <- import_bigWig("/fs/cbsubscb17/storage/data/dm3/s2/procap/PROcap.pl.bw", "/fs/cbsubscb17/storage/data/dm3/s2/procap/PROcap.mn.bw")
HUM_cap <- import_bigWig("/local/workdir/bab364/tmp/procap/K562_grocap_plus.bigWig","/local/workdir/bab364/tmp/procap/K562_grocap_minus.bigWig", genome = "hg19")

#load consensus annotations from "load5prime_FINAL.R"
sp.genes.og
sc.genes.og
dmel.genes.og
Hs.genes.og

#reannotate with procap data
sp.genes.filter <- unique(sp.genes.og %>% filter(width(sp.genes.og) > 400)) # remove short genes, set for pause index threshold to preserve extra genes
sp.genes.procap <- unique(Reannotate_maxPROcap(SPOM_cap, sp.genes.filter, binsize = 10, search.window = 500, field = "score")) # remove artificial duplicates
length(sp.genes.filter)-length(unique(promoters(sp.genes.procap, upstream = 0, downstream = 1))) # check for TSS duplicates
saveRDS(sp.genes.procap, file = "/local/workdir/bab364/Annotations/PROcap/Sp_procap_reann.rds")
sc.genes.filter <- unique(sc.genes.og %>% filter(width(sc.genes.og) > 400)) # remove short genes, set for pause index threshold to preserve extra genes
sc.genes.procap <- unique(Reannotate_maxPROcap(SCER_cap, sc.genes.filter, binsize = 10, search.window = 500, field = "score"))
saveRDS(sc.genes.procap, file = "/local/workdir/bab364/Annotations/PROcap/Sc_procap_reann.rds")
dm.genes.filter <- unique(dmel.genes.og %>% filter(width(dmel.genes.og) > 400))
dm.genes.procap <- unique(Reannotate_maxPROcap(DROS_cap, dm.genes.filter, binsize = 10, search.window = 500, field = "score"))
saveRDS(dm.genes.procap, file = "/local/workdir/bab364/Annotations/PROcap/dm3_procap_reann.rds")
hs.genes.filter <- unique(Hs.genes.og %>% filter(width(Hs.genes.og) > 400))
hs.genes.procap <- unique(Reannotate_maxPROcap(HUM_cap, hs.genes.filter, binsize = 10, search.window = 500, field = "score"))
saveRDS(hs.genes.procap, file = "/local/workdir/bab364/Annotations/PROcap/hg19_procap_reann.rds")
# or load finished annots:
sp.genes.procap <- readRDS(file = "/local/workdir/bab364/Annotations/PROcap/Sp_procap_reann.rds")
sc.genes.procap <- readRDS(file = "/local/workdir/bab364/Annotations/PROcap/Sc_procap_reann.rds")
dm.genes.procap <- readRDS(file = "/local/workdir/bab364/Annotations/PROcap/dm3_procap_reann.rds")
hs.genes.procap <- readRDS(file = "/local/workdir/bab364/Annotations/PROcap/hg19_procap_reann.rds")


#load mnase data
HUM.mnase <- import.bw("/fs/cbsubscb17/storage/data/hg19/k562/sydh_mnase/wgEncodeSydhNsomeK562Sig.bigWig")
HUM.mnase <- HUM.mnase %>% filter(HUM.mnase$score!=0) #corrects null value issue
#HUM.mnase.alt <- import_bedGraph("/local/workdir/bab364/ATACseq/GSE78984_k562.Pooled_WholeChromatin_bin500_tab.bedGraph")
#HUM.mnase.alt <- HUM.mnase.alt %>% filter(HUM.mnase.alt$score!=0) #corrects null value issue

DROS.mnase <- import.bw("/fs/cbsubscb17/storage/data/dm3/s2/mnase/Adelman_MNAse_Untreated.bigWig")
Scer.mnase <- import.bw("/local/workdir/bab364/PauseEvolution/GSM1299403_ISW2K215R_NR_MNase.bw")
seqlevelsStyle(Scer.mnase) <- seqlevelsStyle(SCER_PE)
Spom.mnase <- import.wig("/local/workdir/bab364/PauseEvolution/GSM1374060_WT.wig.gz")

MetaProfilePlotNonBRG(HUM.mnase, promoters(hs.genes.procap, upstream = 500, downstream = 500), TITLE="Human Mnase")
MetaProfilePlotNonBRG(DROS.mnase, promoters(dm.genes.procap, upstream = 500, downstream = 500), TITLE="Drosophlia Mnase")
MetaProfilePlotNonBRG(Scer.mnase, promoters(sc.genes.procap, upstream = 500, downstream = 500), TITLE="Saccer Mnase")
MetaProfilePlotNonBRG(Spom.mnase, promoters(spom.genes.procap, upstream = 500, downstream = 500), TITLE="Pombe Mnase")

#find bin containing max mnase signal for each promoter, starts from whole genes
mnaseMax <- function(mnase, annots, window = 350){
  counts <- getCountsByPositions(mnase, promoters(annots, upstream = 0, downstream = window), binsize = 5, expand_ranges = T)
  maxs <- max.col(counts, ties.method = "first")
  return(maxs[maxs != 1]*5)
}
firstMax <- function(x, window = 350){
  histogram <- hist(x, breaks = seq(0,window, length.out = (window/10)+1))$counts
  firstmax <- which.max(histogram[1:((window/10)-1)])
  if(firstmax<7){
    string2 <- c(rep(NA, firstmax+7), histogram[(firstmax+8):((window/10)-1)])
  } else if(firstmax>=(window/10)-8){
    string2 <- c(histogram[1:(firstmax-7)], rep(NA,(window/10)-firstmax-1))
  } else {
    string2 <- c(histogram[1:(firstmax-7)], rep(NA, 14), histogram[(firstmax+7):((window/10)-1)])
  }
  secondmax <- which.max(string2)
  c(firstmax,secondmax)*10
}

# filter for only expressed genes
sp.genes.procap.exp <- filterEXPfullREANNOTATE(SPOM_PE, sp.genes.procap, exp.level.minimum = .25, hard.minimum = 5, keep.exp = TRUE, filter.readthrough = F)
sc.genes.procap.exp <- filterEXPfullREANNOTATE(SCER_PE, sc.genes.procap, exp.level.minimum = .25, hard.minimum = 5, keep.exp = TRUE, filter.readthrough = F)
dm.genes.procap.exp <- filterEXPfullREANNOTATE(DROS_PE, dm.genes.procap, exp.level.minimum = .25, hard.minimum = 5, keep.exp = TRUE, filter.readthrough = F)
hs.genes.procap.exp <- filterEXPfullREANNOTATE(HUM_PE, hs.genes.procap, exp.level.minimum = .25, hard.minimum = 5, keep.exp = TRUE, filter.readthrough = F)
# check
MetaProfilePlotNonBRG(HUM.mnase, promoters(hs.genes.procap.exp, upstream = 500, downstream = 500), TITLE="Human Mnase")
MetaProfilePlotNonBRG(DROS.mnase, promoters(dm.genes.procap.exp, upstream = 500, downstream = 500), TITLE="Drosophlia Mnase")

# find nucleosome positions
pom.pos <- mnaseMax(Spom.mnase, sp.genes.procap.exp)
scer.pos <- mnaseMax(Scer.mnase, sc.genes.procap.exp)
dros.pos <- mnaseMax(DROS.mnase, dm.genes.procap.exp)
hum.pos <- mnaseMax(HUM.mnase, hs.genes.procap.exp)
# build df
mnase.df <- data.frame(distance = c(pom.pos, scer.pos, dros.pos, hum.pos),
                       species = c(rep("S. pombe", length(pom.pos)), rep("S. cerevisiae", length(scer.pos)), rep("D. melanogaster", length(dros.pos)), rep("H. sapiens", length(hum.pos))))
mnase.df$species <- factor(mnase.df$species, levels = c("H. sapiens","D. melanogaster","S. pombe","S. cerevisiae"))
mnase.fig <- ggplot(mnase.df, aes(x = distance, y = species)) +
  geom_violin(aes(fill = species), alpha = .5) + 
  stat_summary(fun = firstMax, color = c("black", "white","black","white", "black", "white","white", "black" ) )+
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),
        legend.position = "none", 
        legend.text = element_text( color="black", size=10),
        axis.text.y = element_text( color="black", size=10),
        axis.title.x = element_blank()) +  
  labs(title = "Nucleosome positioning", y = "Species", x = "Distance from TSS")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  geom_vline(xintercept = min(firstMax(pom.pos)), linetype = "dashed")
mnase.fig

# Plot proseq with MNase:
comboPlot <- function(nonBRG, BRG, ANNOTATIONS, TITLE, last.plot = TRUE, ncores = 8, scale.to = "max", X.lab="Distance from TSS"){
  ANNOTATIONS <- promoters(ANNOTATIONS, upstream =  150, downstream = 500)
  meansNONbrg <- metaSubsample(nonBRG, ANNOTATIONS, binsize = 5,lower = 0.25, upper = 0.75, ncores = ncores, expand_ranges = T)
  meansBRG <- metaSubsample(BRG, ANNOTATIONS, binsize = 5,lower = 0.25, upper = 0.75, ncores = ncores, expand_ranges = F)
  scaleBRG <- ifelse(scale.to == "max", max(meansNONbrg$upper)/max(meansBRG$upper),max(meansNONbrg$upper)/max(meansBRG$mean)) #if not max, then mean
  meansBRG$upper <- meansBRG$upper*scaleBRG
  meansBRG$lower <- meansBRG$lower*scaleBRG
  meansBRG$mean <- meansBRG$mean*scaleBRG
  means <- rbind(meansNONbrg, meansBRG)
  scalefactor <- round(max(means$upper), 3) # this establishes an upper bound for plotting. setting plus.minus.ratio sets the lower bound (-scalefactor*plus.minus.ratio)
  limits <- ifelse(scale.to == "max", scalefactor*(1.1), max(means$mean)*1.1)
  PLOT <- ggplot(means, aes(x, mean, color = sample.name)) +
    geom_line() + theme_classic() + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    geom_ribbon(aes(x, ymin = lower, ymax = upper,
                    color = NULL, fill = sample.name),
                alpha = 0.2) + ylim(0,scalefactor*(1.1))+
    labs(title = TITLE,
         x = X.lab ,
         y = "Mean Signal + 50% CI") + theme(legend.position = "none")+ 
    scale_x_continuous(labels = function(x) (x)-150, breaks = c(0, 150, 650))+
    coord_cartesian(ylim = c(0,limits))
  return(PLOT)
} #currently doesn't include labeled axes, but scales proseq signal to match mnase maximum

comboPlot(DROS.mnase, DROS_PE, dm.genes.procap.exp, TITLE="Drosophlia Mnase & PRO-seq")
#hs.mnase.genes <- hs.genes %>% filter(getCountsByRegions(HUM.mnase, resize(hs.genes, width = 500, fix = "end"), expand_ranges = T)>0)
comboPlot(HUM.mnase, HUM_PE, hs.genes.procap.exp, TITLE="Human Mnase & PRO-seq")
comboPlot(Scer.mnase, SCER_PE, sc.genes.procap.exp, TITLE="Saccer Mnase & PRO-seq")
comboPlot(Spom.mnase, SPOM_PE, sp.genes.procap.exp, TITLE="Pombe Mnase & PRO-seq")

### make combined clean plot for figure 2 lower pannels
# mnase fig
mnase.df$species <- factor(mnase.df$species, levels = c("S. cerevisiae","S. pombe","D. melanogaster","H. sapiens"))
mnase.fig.pub <- ggplot(mnase.df, aes(x = distance, y = species)) +
  geom_violin(alpha = .5) + 
  stat_summary(fun = firstMax, color = c("gray","black", "black","gray", "black", "gray", "black" ,"gray") )+
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),
        legend.position = "none", 
        legend.text = element_text( color="black", size=10),
        axis.text.y = element_blank(),
        axis.title.x = element_blank()) +  
  labs(title = NULL, y = NULL, x = "Distance from TSS")
mnase.fig.pub
# pi fig
pi20p <- calcPI(hs.genes.procap.exp,HUM_PE, output = "individual")
pi16p <- calcPI(dm.genes.procap.exp,DROS_PE, output = "individual")
pi7p <- calcPI(sp.genes.procap.exp,SPOM_PE, output = "individual")
pi8p <- calcPI(sc.genes.procap.exp,SCER_PE, output = "individual")
pi.df.cap <- data.frame(PI = c(pi20p,pi16p,pi7p,pi8p),
                        Species = c(rep("H. sapiens", length(pi20p)),rep("D. melanogaster", length(pi16p)),rep("S. pombe", length(pi7p)),rep("S. cerevisiae", length(pi8p))))
pi.df.cap$Species <- factor(pi.df.cap$Species, levels = c("S. cerevisiae","S. pombe","D. melanogaster","H. sapiens"))
pi.fig.pub <- ggplot(pi.df.cap, aes(x = log10(PI), y = Species)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),
        legend.position = "none", 
        legend.text = element_text( color="black", size=10),
        axis.text.y = element_text( color="black", size=10),
        axis.title.x = element_blank()) +  
  labs(title = NULL, y = "Species", x = "log10(Pausing Index)")
pi.fig.pub
# stacked pro-seq/mnase-seq
dm.comboplot <- comboPlot(DROS.mnase, DROS_PE, dm.genes.procap.exp, TITLE=NULL, X.lab = NULL)
hs.comboplot <- comboPlot(HUM.mnase, HUM_PE, hs.genes.procap.exp, TITLE=NULL, X.lab = NULL)
sc.comboplot <- comboPlot(Scer.mnase, SCER_PE, sc.genes.procap.exp, TITLE=NULL, X.lab = NULL)
sp.comboplot <- comboPlot(Spom.mnase, SPOM_PE, sp.genes.procap.exp, TITLE=NULL, X.lab = NULL)
laym <- rbind(c(1),
              c(2),
              c(3),
              c(4))
comb.stack <- grid.arrange(hs.comboplot, dm.comboplot, sp.comboplot, sc.comboplot, layout_matrix=laym)
# combine all
comb.all <- plot_grid(pi.fig.pub,mnase.fig.pub,comb.stack, align = "v", nrow = 1)
comb.all
ggsave2("/local/workdir/bab364/PauseEvolution/fig3new.pdf", plot = comb.all, device = "pdf", width = 8.2, height = 3.5, units = "in") #export to match size of OG figures, allow direct substitution


