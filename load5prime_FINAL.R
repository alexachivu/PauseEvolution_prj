####################################################################################################################################
# this script loads 5' proseq data and consensus annotations
# performs filtering & TSS re-annotation
# saves new annotations as .RDS in 5' reann folder
# dumps memory after each species 
# required functions at bottom of script
####################################################################################################################################
library(zoo)
library(parallel)
library(Biostrings)
library(BRGenomics)
library(rtracklayer)
library(BSgenome)
library(tidyverse)
library(plyranges)

# E.coli
#load TSS
ec.genes <- read.table( "/fs/cbsubscb17/storage/data/short_read_index/ec_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.46.gff3", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(ec.genes) = c( "chr", "Transcript.type", "start", "end", "strand" )
ec.genes = ec.genes[ ec.genes$Transcript.type == "gene", ]
ec.genes = as(ec.genes, "GRanges")
seqlevels(ec.genes) <- c("U00096.2")
#load5p
#ECOLI_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/E.coli/My_output-05_23_2023/Ec-C123_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/E.coli/My_output-05_23_2023/Ec-C123_minus.bw")
ECOLI_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ecoli/Ec-merge.bam",revcomp = F, trim.to = "5p", paired_end = T)
#refineTSS
ec.genes.filter <- ec.genes %>% filter(width(ec.genes) > 750) # remove short genes
ec.genes.filter <- filterEXPfullREANNOTATE(ECOLI_5p, ec.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
ec.genes.filter <- unique(ec.genes.filter)
paste0(length(ec.genes)-length(ec.genes.filter),"(",round((length(ec.genes)-length(ec.genes.filter))/length(ec.genes)*100, 2),"%)"," genes removed by filters. ", length(ec.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
ECOLI_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(ECOLI_5p, ec.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(ec.genes.filter)-length(unique(promoters(ECOLI_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
median(start(resize(ec.genes.filter,width = 1, fix = "start"))-start(resize(ECOLI_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(ECOLI_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Ecoli.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# H. mediterranei
#load TSS
hm.genes <- read.table( "/fs/cbsubscb17/storage/data/short_read_index/Haloferax_archaea/Haloferax.gff", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(hm.genes) = c( "chr", "Transcript.type", "start", "end", "strand" )
hm.genes = hm.genes[ hm.genes$Transcript.type == "gene", ]
hm.genes = as(hm.genes, "GRanges")
#load5p
#HALO_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/H.mediterranei/halo_5p_plus_R.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/H.mediterranei/halo_5p_minus_R.bw")
HALO_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Hmediterranei/Hmedi-merge.sort.bam",revcomp = T, trim.to = "5p", paired_end = T)
#refineTSS
hm.genes.filter <- hm.genes %>% filter(width(hm.genes) > 750) # remove short genes
hm.genes.filter <- filterEXPfullREANNOTATE(HALO_5p, hm.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
hm.genes.filter <- unique(hm.genes)
paste0(length(hm.genes)-length(hm.genes.filter),"(",round((length(hm.genes)-length(hm.genes.filter))/length(hm.genes)*100, 2),"%)"," genes removed by filters. ", length(hm.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
HALO_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(HALO_5p, hm.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(hm.genes.filter)-length(unique(promoters(HALO_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
mean(start(resize(hm.genes.filter,width = 1, fix = "start"))-start(resize(HALO_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
metaProfilePauseEvolution(hm.genes.filter, HALO_5p)
metaProfilePauseEvolution(HALO_5p_step_reannot, HALO_5p)
#save new TSS
saveRDS(HALO_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Hmediterranei.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# A. thaliana ##### no data maps to minus strand using RunOnBamToBW. Imported direct from BAM instead #####
#load TSS
athal.TSS = read.table( "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/A.thaliana/Arabidopsis_thaliana.TAIR10.56.gff3", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
#athal.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/A_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(athal.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
athal.TSS = athal.TSS[ athal.TSS$Transcript.type == "gene", ]
athal.TSS = as(athal.TSS, "GRanges")
#load5p
ARA_5plus <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/A.thaliana/My_output-05_23_2023/Arabidopsis_Gro_12_plus.bw")
#ARA_5minus <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/A.thaliana/My_output-05_23_2023/Arabidopsis_Gro_12_minus.bw") #Error in seqinfo(ranges) : UCSC library operation failed
#ARA_5p <- mergeGRangesData (ARA_5plus, ARA_5minus)
#rm(ARA_5plus, ARA_5minus)
ARA_5p <- import_bam("/fs/cbsubscb17/storage/projects/NASA_2020/Arabidopsis_Gro_1/Arabidopsis_Gro_1_QC.sort.bam",revcomp = F, trim.to = "5p")
ARA_3p <- import_bam("/fs/cbsubscb17/storage/projects/NASA_2020/Arabidopsis_Gro_1/Arabidopsis_Gro_1_QC.sort.bam",revcomp = F, trim.to = "3p")
ARA_5p_2 <- import_bam("/fs/cbsubscb17/storage/projects/NASA_2020/Arabidopsis_Gro_2/Arabidopsis_Gro_2_QC.sort.bam",revcomp = F, trim.to = "5p")
ARA_5p <- mergeGRangesData(ARA_5p,ARA_5p_2)
MetaProfilePlotMulti_sample(ARA_5p, ARA_5p_2, ARA_3p, ANNOTATIONS = promoters(athal.TSS, upstream = 500, downstream = 500), SPECIES = "A.thaliana, import tests")
#refineTSS
athal.genes.filter <- athal.TSS %>% filter(width(athal.TSS) > 750) # remove short genes
athal.genes.filter <- filterEXPfullREANNOTATE(ARA_5p, athal.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
athal.genes.filter <- unique(athal.genes.filter)
paste0(length(athal.TSS)-length(athal.genes.filter),"(",round((length(athal.TSS)-length(athal.genes.filter))/length(athal.TSS)*100, 2),"%)"," genes removed by filters. ", length(athal.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
ARA_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(ARA_5p, athal.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(athal.genes.filter)-length(unique(promoters(ARA_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
median(start(resize(athal.genes.filter,width = 1, fix = "start"))-start(resize(ARA_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
MetaProfilePlotMulti_sample(ARA_5p, ARA_3p, ANNOTATIONS = promoters(ARA_5p_step_reannot, upstream = 500, downstream = 500), SPECIES = "A.thaliana, import tests")

#save new TSS
saveRDS(ARA_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Athaliana.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# Z. mays ##### hard minimum kicks in on expression filter using 5map_IGV, not when using Map5p.end from james #####
#load TSS
ZMAYS.nc.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/Zea.mays/genome_assemblies_genome_gff/ncbi-genomes-2022-09-22/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(ZMAYS.nc.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
ZMAYS.nc.TSS = ZMAYS.nc.TSS[ ZMAYS.nc.TSS$Transcript.type == "gene", ]
ZMAYS.nc.TSS = as(ZMAYS.nc.TSS, "GRanges")
#MetaProfilePlotMany(ZMAYS_PE2, Z.mays_reann.TSS, ZMAYS.nc.TSS, SPECIES = "Z.mays, Reann vs Gencode", plus.minus.ratio = .5)
#load5p
#ZMAYS_5p<- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/Z.mays/ZMAYS_5map_IGV_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/Z.mays/ZMAYS_5map_IGV_minus.bw") #SE
ZMAYS_5p2 <- import_bigWig("/fs/cbsubscb17/storage/data/short_read_index/Zea.mays/Z.mays_Map5p.end_plus.bw", "/fs/cbsubscb17/storage/data/short_read_index/Zea.mays/Z.mays_Map5p.end_minus.bw") #SE
MetaProfilePlotMulti_sample(ZMAYS_PE2, ZMAYS_5p2, ANNOTATIONS = promoters(ZMAYS.nc.TSS, upstream = 500, downstream = 500), SPECIES = "Z.mays, final tests")

#refineTSS
ZMAYS.genes.filter <- ZMAYS.nc.TSS %>% filter(width(ZMAYS.nc.TSS) > 750) # remove short genes
ZMAYS.genes.filter <- filterEXPfullREANNOTATE(ZMAYS_5p2, ZMAYS.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
ZMAYS.genes.filter <- unique(ZMAYS.genes.filter)
paste0(length(ZMAYS.nc.TSS)-length(ZMAYS.genes.filter),"(",round((length(ZMAYS.nc.TSS)-length(ZMAYS.genes.filter))/length(ZMAYS.nc.TSS)*100, 2),"%)"," genes removed by filters. ", length(ZMAYS.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
ZMAYS_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(ZMAYS_5p2, ZMAYS.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(ZMAYS.genes.filter)-length(unique(promoters(ZMAYS_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(ZMAYS.genes.filter)-length(unique(ZMAYS_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(ZMAYS.genes.filter,width = 1, fix = "start"))-start(resize(ZMAYS_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(ZMAYS_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Zmays.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# O. sativa
#load TSS
osat.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/O.sativa/O.sativa_GCF_001433935.1_IRGSP-1.0_genomic.gff", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(osat.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
osat.TSS = osat.TSS[ osat.TSS$Transcript.type == "gene", ]
osat.TSS = as(osat.TSS, "GRanges")
#load5p
OSAT_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/O.sativa/Osat_5p_plus_R.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/O.sativa/Osat_5p_minus_R.bw") #SE
#refineTSS
osat.genes.filter <- osat.TSS %>% filter(width(osat.TSS) > 750) # remove short genes
osat.genes.filter <- filterEXPfullREANNOTATE(OSAT_5p, osat.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
osat.genes.filter <- unique(osat.genes.filter)
paste0(length(osat.TSS)-length(osat.genes.filter),"(",round((length(osat.TSS)-length(osat.genes.filter))/length(osat.TSS)*100, 2),"%)"," genes removed by filters. ", length(osat.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
OSAT_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(OSAT_5p, osat.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(osat.genes.filter)-length(unique(promoters(OSAT_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
median(start(resize(osat.genes.filter,width = 1, fix = "start"))-start(resize(OSAT_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(OSAT_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Osativa.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# D. discoideum ##### some TSS are consolidated by re-annotation. these have different PAS, so they are retained in current workflow #####
#load TSS
dicty.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/dictyDisco/bwa/usr/local/dicty/data/gff3/dicty.combined.gtf", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(dicty.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
dicty.TSS = dicty.TSS[ dicty.TSS$Transcript.type == "transcript", ]
dicty.TSS = as(dicty.TSS, "GRanges")
#load5p
#DIC_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/D.discoideum/5prime_output-12_20_2023/dicty_r1_QC_end.sort_plus.bw","/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/D.discoideum/5prime_output-12_20_2023/dicty_r1_QC_end.sort_minus.bw")
DIC_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ddiscoideum/Ddisco-merge.sort.bam",revcomp = F, trim.to = "5p", paired_end = T)
#refineTSS
dicty.genes.filter <- dicty.TSS %>% filter(width(dicty.TSS) > 750) # remove short genes
dicty.genes.filter <- filterEXPfullREANNOTATE(DIC_5p, dicty.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
dicty.genes.filter <- unique(dicty.genes.filter)
paste0(length(dicty.TSS)-length(dicty.genes.filter),"(",round((length(dicty.TSS)-length(dicty.genes.filter))/length(dicty.TSS)*100, 2),"%)"," genes removed by filters. ", length(dicty.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
DIC_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(DIC_5p, dicty.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(dicty.genes.filter)-length(unique(promoters(DIC_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(dicty.genes.filter)-length(unique(DIC_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(dicty.genes.filter,width = 1, fix = "start"))-start(resize(DIC_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(DIC_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Ddiscoideum.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# S. pombe ##### issues with 5' mapping. Seems to be the protocol used for this sample- because reverse complement is needed, runonbamtobw does not have a built in 5' mapper. confirmed that not all SE samples are affected #####
#load TSS
sp.genes.og <- read.table( "/fs/cbsubscb17/storage/data/short_read_index/S_pombe/Schizosaccharomyces_pombe.ASM294v2.46.gff3", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(sp.genes.og) = c( "chr", "Transcript.type", "start", "end", "strand" )
sp.genes.og = sp.genes.og[ sp.genes.og$Transcript.type == "mRNA", ]
sp.genes.og = as(sp.genes.og, "GRanges")
#sp.TSS <- promoters(sp.genes.og, upstream = 500, downstream = 500)
#load5p
SPOM_5p <- import_bam("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/S.pombe/Pombe_WT_12.bam",revcomp = T, trim.to = "5p", paired_end = F)
#refineTSS
sp.genes.filter <- sp.genes.og %>% filter(width(sp.genes.og) > 750) # remove short genes
sp.genes.filter <- filterEXPfullREANNOTATE(SPOM_5p, sp.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
sp.genes.filter <- unique(sp.genes.filter)
paste0(length(sp.genes.og)-length(sp.genes.filter),"(",round((length(sp.genes.og)-length(sp.genes.filter))/length(sp.genes.og)*100, 2),"%)"," genes removed by filters. ", length(sp.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
SPOM_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(SPOM_5p, sp.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(sp.genes.filter)-length(unique(promoters(SPOM_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
metaProfilePauseEvolution(SPOM_5p_step_reannot, SPOM_5p)
metaProfilePauseEvolution(SPOM_5p_step_reannot, SPOM_PE)
#verify
sp.TSS <- promoters(sp.genes.og, upstream = 500, downstream = 500)
sp.reann.TSS <- promoters(SPOM_5p_step_reannot, upstream = 500, downstream = 500)
MetaProfilePlotMany(SPOM_5p, sp.TSS, sp.reann.TSS, SPECIES = "S. pombe", plus.minus.ratio = .2, aspect.ratio = 1/2)
MetaProfilePlotMany(SPOM_PE, sp.TSS, sp.reann.TSS, SPECIES = "S. pombe", plus.minus.ratio = .2, aspect.ratio = 1/2)
MetaProfilePlotMany(SPOM_cap, sp.TSS, sp.reann.TSS, SPECIES = "S. pombe", plus.minus.ratio = .2, aspect.ratio = 1/2)
#save new annotations
saveRDS(SPOM_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Spombe.step.reann.rds")
#sp.load <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Spombe.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# S. cerevisiae
#load TSS
sc.genes.og <- read.table( "/fs/cbsubscb17/storage/data/short_read_index/S_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.99.gff3", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(sc.genes.og) = c( "chr", "Transcript.type", "start", "end", "strand" )
sc.genes.og = sc.genes.og[ sc.genes.og$Transcript.type == "gene", ]
sc.genes.og = as(sc.genes.og, "GRanges")
#load5p
SCER_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/S.cerevisiae/My_output-05_22_2023/S.cerevisiae_Gro_12_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/S.cerevisiae/My_output-05_22_2023/S.cerevisiae_Gro_12_minus.bw") #SE
#refineTSS
sc.genes.filter <- sc.genes.og %>% filter(width(sc.genes.og) > 750) # remove short genes
sc.genes.filter <- filterEXPfullREANNOTATE(SCER_5p, sc.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
sc.genes.filter <- unique(sc.genes.filter)
paste0(length(sc.genes.og)-length(sc.genes.filter),"(",round((length(sc.genes.og)-length(sc.genes.filter))/length(sc.genes.og)*100, 2),"%)"," genes removed by filters. ", length(sc.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
SCER_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(SCER_5p, sc.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(sc.genes.filter)-length(unique(promoters(SCER_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(sc.genes.filter)-length(unique(SCER_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(sc.genes.filter,width = 1, fix = "start"))-start(resize(SCER_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(SCER_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Scerevisiae.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# S. arctica ##### seqlevel mismatch corrected in .converted. files within short_read_index directory #####
#load TSS
sarct.TSS <- import.bed("/fs/cbsubscb17/storage/data/short_read_index/S_arctica/Sphaeroforma_arctica_jp610_gca_001186125.Spha_arctica_JP610_V1.46.converted.gff3_mRNA_peppro.bed")
#load5p
#SARCT_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/S.arctica/My_output-05_23_2023/sphaeroforma_r12_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/S.arctica/My_output-05_23_2023/sphaeroforma_r12_minus.bw")
SARCT_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Sarctica/Sarctica-merge.sort.bam",revcomp = F, trim.to = "5p", paired_end = T)
#refineTSS
sa.genes.filter <- sarct.TSS %>% filter(width(sarct.TSS) > 1000) # remove short genes
sa.genes.filter <- filterEXPfullREANNOTATE(SARCT_5p, sa.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
sa.genes.filter <- unique(sa.genes.filter)
paste0(length(sarct.TSS)-length(sa.genes.filter),"(",round((length(sarct.TSS)-length(sa.genes.filter))/length(sarct.TSS)*100, 2),"%)"," genes removed by filters. ", length(sa.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
SARCT_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(SARCT_5p, sa.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(sa.genes.filter)-length(unique(promoters(SARCT_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(sa.genes.filter)-length(unique(SARCT_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(sa.genes.filter,width = 1, fix = "start"))-start(resize(SARCT_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(SARCT_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Sarctica.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# C. fragrantissima ##### alternate annotations available from Michelle (not checked yet) #####
#load TSS
cfrag.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/C_fragrantissima/Creolimax_fragrantissima.gtf", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(cfrag.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
cfrag.TSS = cfrag.TSS[ cfrag.TSS$Transcript.type == "gene", ]
cfrag.TSS = as(cfrag.TSS, "GRanges")
#load5p
#CFRAG_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/C.fragrantissima/My_output-05_23_2023/creolimax_r12_plus.bw","/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/C.fragrantissima/My_output-05_23_2023/creolimax_r12_minus.bw")
CFRAG_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Cfragrantissima/Cfragrantissima-merge.sort.bam",revcomp = F, trim.to = "5p", paired_end = T)
#refineTSS
cf.genes.filter <- cfrag.TSS %>% filter(width(cfrag.TSS) > 1000) # remove short genes
cf.genes.filter <- filterEXPfullREANNOTATE(CFRAG_5p, cf.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
cf.genes.filter <- unique(cf.genes.filter)
paste0(length(cfrag.TSS)-length(cf.genes.filter),"(",round((length(cfrag.TSS)-length(cf.genes.filter))/length(cfrag.TSS)*100, 2),"%)"," genes removed by filters. ", length(cf.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
CFRAG_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(CFRAG_5p, cf.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(cf.genes.filter)-length(unique(promoters(CFRAG_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(cf.genes.filter)-length(unique(CFRAG_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(cf.genes.filter,width = 1, fix = "start"))-start(resize(CFRAG_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(CFRAG_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Cfragrantissima.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# C. owczarzaki
#load TSS
cap.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/C_owczarzaki/Capsaspora_owczarzaki_atcc_30864_gca_000151315.C_owczarzaki_V2.46.converted.gff3", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(cap.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
cap.TSS = cap.TSS[ cap.TSS$Transcript.type == "gene", ]
cap.TSS = as(cap.TSS, "GRanges")
#load5p
#CAP_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/C.owcz/My_output-05_23_2023/capsaspora_r12_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/C.owcz/My_output-05_23_2023/capsaspora_r12_minus.bw")
CAP_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Cowczarzaki/Cowczarzaki-merge.sort.bam",revcomp = F, trim.to = "5p", paired_end = T)

#refineTSS
co.genes.filter <- cap.TSS %>% filter(width(cap.TSS) > 1000) # remove short genes
co.genes.filter <- filterEXPfullREANNOTATE(CAP_5p, co.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
co.genes.filter <- unique(co.genes.filter)
paste0(length(cap.TSS)-length(co.genes.filter),"(",round((length(cap.TSS)-length(co.genes.filter))/length(cap.TSS)*100, 2),"%)"," genes removed by filters. ", length(co.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
CAP_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(CAP_5p, co.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(co.genes.filter)-length(unique(promoters(CAP_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(co.genes.filter)-length(unique(CAP_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(co.genes.filter,width = 1, fix = "start"))-start(resize(CAP_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(CAP_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Cowczarzaki.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# N. vectensis ##### alternate TSS set available #####
#load TSS
nemvect.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/nemVec1/NemVec.AM20922v1.46.gff3", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(nemvect.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
nemvect.TSS = nemvect.TSS[ nemvect.TSS$Transcript.type == "CDS", ]
nemvect.TSS = as(nemvect.TSS, "GRanges")
#load5p
#NEM_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/N.vectensis/NEM_5p_plus_R.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/N.vectensis/NEM_5p_minus_R.bw")
NEM_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Nvectensis/Nvectensis-merge.sort.bam",revcomp = F, trim.to = "5p", paired_end = T)

#refineTSS
nv.genes.filter <- nemvect.TSS %>% filter(width(nemvect.TSS) > 750) # remove short genes
nv.genes.filter <- filterEXPfullREANNOTATE(NEM_5p, nv.genes.filter, exp.level.minimum = 0.1, hard.minimum = 3, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
nv.genes.filter <- unique(nv.genes.filter)
paste0(length(nemvect.TSS)-length(nv.genes.filter),"(",round((length(nemvect.TSS)-length(nv.genes.filter))/length(nemvect.TSS)*100, 2),"%)"," genes removed by filters. ", length(nv.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
NEM_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(NEM_5p, nv.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(nv.genes.filter)-length(unique(promoters(NEM_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(nv.genes.filter)-length(unique(NEM_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(nv.genes.filter,width = 1, fix = "start"))-start(resize(NEM_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#save new TSS
saveRDS(NEM_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Nvectensis.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# C. elegans ##### many whole genes are collapsed by reannotation. These have been de-duplicated. Many TSS pileups remain #####
#load TSS
cel.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/C_elegans/ce6.ensGene.gtf", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(cel.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
cel.TSS = cel.TSS[ cel.TSS$Transcript.type == "transcript", ]
cel.TSS = as(cel.TSS, "GRanges")
seqlevelsStyle(cel.TSS) <- "Ensembl"
#load5p
CEL_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/C.elegans/Celegans_CHROMSmatchGFF_5p_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/C.elegans/Celegans_CHROMSmatchGFF_5p_minus.bw") #SE
#refineTSS
ce.genes.filter <- cel.TSS %>% filter(width(cel.TSS) > 750) # remove short genes
ce.genes.filter <- filterEXPfullREANNOTATE(CEL_5p, ce.genes.filter, exp.level.minimum = 0.1, hard.minimum = 3, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
ce.genes.filter <- unique(ce.genes.filter)
paste0(length(cel.TSS)-length(ce.genes.filter),"(",round((length(cel.TSS)-length(ce.genes.filter))/length(cel.TSS)*100, 2),"%)"," genes removed by filters. ", length(ce.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
CEL_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(CEL_5p, ce.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(ce.genes.filter)-length(unique(promoters(CEL_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(ce.genes.filter)-length(unique(CEL_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(ce.genes.filter,width = 1, fix = "start"))-start(resize(CEL_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
CEL_5p_step_reannot <- unique(CEL_5p_step_reannot) #consolidate collapsed annotations
#save new TSS
saveRDS(CEL_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Celegans.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# D. pulex ##### something went wrong with 5p mapping-no reads produced. Import 5' ends directly from bam, checked correct 5' ends reported #####
#load TSS
daphnia.TSS = read.table( "/fs/cbsubscb17/storage/data/short_read_index/dpulex_jgi060905/dpulex-all-genes-jgi060905.gff", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(daphnia.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
daphnia.TSS = daphnia.TSS[ daphnia.TSS$Transcript.type == "gene", ]
daphnia.TSS = as(daphnia.TSS, "GRanges")
#load5p
#DAPH_5p <- import_bam("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/D.pulex/Daphnia_Part12.bam",revcomp = F, trim.to = "5p", paired_end = T)
DAPH_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Dpulex/Dpulex-merge.sort.bam",revcomp = F, trim.to = "5p", paired_end = T)

#refineTSS
dp.genes.filter <- daphnia.TSS %>% filter(width(daphnia.TSS) > 750) # remove short genes
dp.genes.filter <- filterEXPfullREANNOTATE(DAPH_5p, dp.genes.filter, exp.level.minimum = 0.1, hard.minimum = 3, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
dp.genes.filter <- unique(dp.genes.filter)
paste0(length(daphnia.TSS)-length(dp.genes.filter),"(",round((length(daphnia.TSS)-length(dp.genes.filter))/length(daphnia.TSS)*100, 2),"%)"," genes removed by filters. ", length(dp.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
DAPH_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(DAPH_5p, dp.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(dp.genes.filter)-length(unique(promoters(DAPH_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(dp.genes.filter)-length(unique(DAPH_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(dp.genes.filter,width = 1, fix = "start"))-start(resize(DAPH_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#DAPH_5p_step_reannot <- unique(DAPH_5p_step_reannot) #consolidate collapsed annotations
#save new TSS
saveRDS(DAPH_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Dpulex.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# D. iulia ##### something went wrong with 5p mapping-no reads produced. Import 5' ends directly from bam, checked correct 5' ends reported #####
di.genes <- read.table( "/local/workdir/James/PauseEvolution/data/D.iulia/Dryas_julia.v0.7.CUSTOM.CDS.sorted.gtf", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(di.genes) = c( "chr", "Transcript.type", "start", "end", "strand" )
di.genes = di.genes[ di.genes$Transcript.type == "transcript", ]
di.genes = as(di.genes, "GRanges")
#load TSS
DIU_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Diulia/Diulia-merge.sort.bam",revcomp = F, trim.to = "5p", paired_end = T)

#refineTSS
di.genes.filter <- di.genes %>% filter(width(di.genes) > 750) # remove short genes
di.genes.filter <- filterEXPfullREANNOTATE(DIU_5p, di.genes.filter, exp.level.minimum = 0.1, hard.minimum = 3, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
di.genes.filter <- unique(di.genes.filter)
paste0(length(di.genes)-length(di.genes.filter),"(",round((length(di.genes)-length(di.genes.filter))/length(di.genes)*100, 2),"%)"," genes removed by filters. ", length(di.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
DIU_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(DIU_5p, di.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(di.genes.filter)-length(unique(promoters(DIU_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(di.genes.filter)-length(unique(DIU_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(di.genes.filter,width = 1, fix = "start"))-start(resize(DIU_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
MetaProfilePlotMulti_sample(DIU_5p, DIU_3p, ANNOTATIONS = promoters(DIU_5p_step_reannot, upstream = 500, downstream = 500), SPECIES = "D.iulia, reannotation tests")
DIU_5p_step_reannot <- unique(DIU_5p_step_reannot) #consolidate collapsed annotations
MetaProfilePlotMulti_sample(DIU_5p, DIU_3p, ANNOTATIONS = promoters(DIU_5p_step_reannot, upstream = 500, downstream = 500), SPECIES = "D.iulia, reannotation tests")
#save new TSS
saveRDS(DIU_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Diulia.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# D. melanogaster
#load TSS
dmel.genes.og <- read.table("/fs/cbsubscb17/storage/data/short_read_index/dm3/dm3.ensGene.gtf", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(dmel.genes.og) = c( "chr", "Transcript.type", "start", "end", "strand" )
dmel.genes.og = dmel.genes.og[ dmel.genes.og$Transcript.type == "transcript", ]
dmel.genes.og = as(dmel.genes.og, "GRanges")
#load5p
DROS_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/D.melanogaster/My_output-05_25_2023/Drosophila_Gro_QC.sort_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/D.melanogaster/My_output-05_25_2023/Drosophila_Gro_QC.sort_minus.bw") #SE
# from bam files of se proseq
dm.5p1 <- import_bam("/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_1/Drosophilia_01A_1_QC.sort.bam",revcomp = T, trim.to = "5p", paired_end = F)
dm.5p2 <- import_bam("/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_2/Drosophilia_01A_2_QC.sort.bam",revcomp = T, trim.to = "5p", paired_end = F)
dm.5p3 <- import_bam("/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_3/Drosophilia_01A_3_QC.sort.bam",revcomp = T, trim.to = "5p", paired_end = F)
dm.5p4 <- import_bam("/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_4/Drosophilia_01A_4_QC.sort.bam",revcomp = T, trim.to = "5p", paired_end = F)
DROS_5p <- mergeGRangesData(dm.5p1,dm.5p2,dm.5p3,dm.5p4)
rm(dm.5p1,dm.5p2,dm.5p3,dm.5p4)
#refineTSS
dmel.genes.filter <- dmel.genes.og %>% filter(width(dmel.genes.og) > 750) # remove short genes
dmel.genes.filter <- filterEXPfullREANNOTATE(DROS_5p, dmel.genes.filter, exp.level.minimum = 0.1, hard.minimum = 5, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
dmel.genes.filter <- unique(dmel.genes.filter)
paste0(length(dmel.genes.og)-length(dmel.genes.filter),"(",round((length(dmel.genes.og)-length(dmel.genes.filter))/length(dmel.genes.og)*100, 2),"%)"," genes removed by filters. ", length(dmel.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
DROS_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(DROS_5p, dmel.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(dmel.genes.filter)-length(unique(promoters(DROS_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(dmel.genes.filter)-length(unique(DROS_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(dmel.genes.filter,width = 1, fix = "start"))-start(resize(DROS_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
DROS_5p_step_reannot <- unique(DROS_5p_step_reannot) #consolidate collapsed annotations
#save new TSS
saveRDS(DROS_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Dmelanogaster.se.step.reann.rds")

saveRDS(DROS_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Dmelanogaster.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# S. purpuratus
#load TSS
spur.genes <- read.table( "/fs/cbsubscb17/storage/data/short_read_index/sp5.0/GCF_000002235.5/genomic.gtf", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(spur.genes) = c( "chr", "Transcript.type", "start", "end", "strand" )
spur.genes = spur.genes[ spur.genes$Transcript.type == "gene", ]
spur.genes = as(spur.genes, "GRanges")
spur.TSS <- promoters(spur.genes, upstream = 500, downstream = 500)
#load5p
URCH_5p_rc <- import_bam("/fs/cbsubscb17/storage/data/short_read_index/sp5.0/mapped_PROseq/merged/Sea_Urchin_Embryos_12_20.merge.sort.bam",revcomp = T, trim.to = "5p", paired_end = T)
#refineTSS
spur.genes.filter <- spur.genes %>% filter(width(spur.genes) > 750) # remove short genes
spur.genes.filter <- filterEXPfullREANNOTATE(URCH_5p_rc, spur.genes.filter, exp.level.minimum = 0.1, hard.minimum = 3, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
spur.genes.filter <- unique(spur.genes.filter)
paste0(length(spur.genes)-length(spur.genes.filter),"(",round((length(spur.genes)-length(spur.genes.filter))/length(spur.genes)*100, 2),"%)"," genes removed by filters. ", length(spur.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
URCH_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(URCH_5p_rc, spur.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(spur.genes.filter)-length(unique(promoters(URCH_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(spur.genes.filter)-length(unique(URCH_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(spur.genes.filter,width = 1, fix = "start"))-start(resize(URCH_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
URCH_5p_step_reannot <- unique(URCH_5p_step_reannot) #consolidate collapsed annotations
#save new TSS
saveRDS(URCH_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Spurpuratus5.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# P. marinus
#load TSS
Pm2.genes <- read.table( "/fs/cbsubscb17/storage/data/short_read_index/petMar2/Petromyzon_marinus.Pmarinus_7.0.99.gtf", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(Pm2.genes) = c( "chr", "Transcript.type", "start", "end", "strand" )
Pm2.genes = Pm2.genes[ Pm2.genes$Transcript.type == "gene", ]
Pm2.genes = as(Pm2.genes, "GRanges")
pm2.tss <- promoters(Pm2.genes, upstream = 500, downstream = 500)
#load5p
LAM_5p <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Pmarinus/Pmarinus-merge.sort.bam",revcomp = F, trim.to = "5p", paired_end = T)
#refineTSS
Pm2.genes.filter <- Pm2.genes %>% filter(width(Pm2.genes) > 750) # remove short genes
Pm2.genes.filter <- filterEXPfullREANNOTATE(LAM_5p, Pm2.genes.filter, exp.level.minimum = 0.1, hard.minimum = 3, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
Pm2.genes.filter <- unique(Pm2.genes.filter)
paste0(length(Pm2.genes)-length(Pm2.genes.filter),"(",round((length(Pm2.genes)-length(Pm2.genes.filter))/length(Pm2.genes)*100, 2),"%)"," genes removed by filters. ", length(Pm2.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
LAM_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(LAM_5p, Pm2.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(Pm2.genes.filter)-length(unique(promoters(LAM_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(Pm2.genes.filter)-length(unique(LAM_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(Pm2.genes.filter,width = 1, fix = "start"))-start(resize(LAM_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#LAM_5p_step_reannot <- unique(LAM_5p_step_reannot) #consolidate collapsed annotations
#save new TSS
saveRDS(LAM_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Pmarinus.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# M. musculus
#load TSS
mouse.TSS = read.table( "/fs/cbsubscb17/storage/data/mm10/gencode/gencode.vM20.annotation.gtf", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(mouse.TSS) = c( "chr", "Transcript.type", "start", "end", "strand" )
mouse.TSS = mouse.TSS[ mouse.TSS$Transcript.type == "gene", ]
mouse.TSS = as(mouse.TSS, "GRanges")
#load5p
MOUSE_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/M.musculus/My_output-05_23_2023/Abood123_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/M.musculus/My_output-05_23_2023/Abood123_minus.bw")
#refineTSS
mouse.genes.filter <- mouse.TSS %>% filter(width(mouse.TSS) > 750) # remove short genes
mouse.genes.filter <- filterEXPfullREANNOTATE(MOUSE_5p, mouse.genes.filter, exp.level.minimum = 0.1, hard.minimum = 3, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
mouse.genes.filter <- unique(mouse.genes.filter)
paste0(length(mouse.TSS)-length(mouse.genes.filter),"(",round((length(mouse.TSS)-length(mouse.genes.filter))/length(mouse.TSS)*100, 2),"%)"," genes removed by filters. ", length(mouse.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
mouse_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(MOUSE_5p, mouse.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(mouse.genes.filter)-length(unique(promoters(mouse_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(mouse.genes.filter)-length(unique(mouse_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(mouse.genes.filter,width = 1, fix = "start"))-start(resize(mouse_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
#mouse_5p_step_reannot <- unique(mouse_5p_step_reannot) #consolidate collapsed annotations
#save new TSS
saveRDS(mouse_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Mmusculus.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))

# H. sapiens (using ChROseq data, not Core/Martins data)
#load TSS
Hs.genes.og <- read.table( "/fs/cbsubscb17/storage/data/short_read_index/hg19/gencode.v19.annotation.gff3", header = FALSE, sep = '\t'  )[, c(1, 3:5, 7)]
colnames(Hs.genes.og) = c( "chr", "Transcript.type", "start", "end", "strand" )
Hs.genes.og = Hs.genes.og[ Hs.genes.og$Transcript.type == "gene", ]
Hs.genes.og = as(Hs.genes.og, "GRanges")
Hs.TSS <- promoters(Hs.genes.og,upstream = 500, downstream = 500)
#load5p
HUM_5p <- import_bam("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/H.sapiens/ChROseq_merged_0h.bam",revcomp = F, trim.to = "5p", paired_end = T)
#refineTSS
hum.genes.filter <- Hs.genes.og %>% filter(width(Hs.genes.og) > 750) # remove short genes
hum.genes.filter <- filterEXPfullREANNOTATE(HUM_5p, hum.genes.filter, exp.level.minimum = 0.1, hard.minimum = 3, keep.exp = T, filter.readthrough = F) #remove non-expressed genes (measured at TSS- really we're asking "is there enough data to refine this TSS" not "is this gene expressed")
hum.genes.filter <- unique(hum.genes.filter)
paste0(length(Hs.genes.og)-length(hum.genes.filter),"(",round((length(Hs.genes.og)-length(hum.genes.filter))/length(Hs.genes.og)*100, 2),"%)"," genes removed by filters. ", length(hum.genes.filter), " remain.")
searchWindow <- 500 #in Reannotate... search.window is the total region to search (so a window of 1kb means +/- 500bp from expected TSS)
hum_5p_step_reannot <- Reannotate_TSS_stepFunc_GR(HUM_5p, hum.genes.filter, binsize = 10, search.window = searchWindow, step = 20, name="Test", field = "score")
paste0(length(hum.genes.filter)-length(unique(promoters(hum_5p_step_reannot,upstream = 0,downstream = 1)))," alternative TSSs consolidated by re-annotation")
paste0(length(hum.genes.filter)-length(unique(hum_5p_step_reannot))," whole genes consolidated by re-annotation")
mean(start(resize(hum.genes.filter,width = 1, fix = "start"))-start(resize(hum_5p_step_reannot,width = 1, fix = "start"))) #report average move distance
hum_5p_step_reannot <- unique(hum_5p_step_reannot) #consolidate collapsed annotations
#save new TSS
saveRDS(hum_5p_step_reannot, file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Hsapiens.step.reann.rds")
#clear envt
rm(list = setdiff(ls(), lsf.str()))


####################################################################################################################################
# create step function & helper functions
####################################################################################################################################

# version of filterEXP that considers only the promoter region (-500 to 500) and outputs full genes
filterEXPfullREANNOTATE <- function(DATA, genes, exp.level.minimum = .5, hard.minimum = 1, keep.exp = TRUE, filter.readthrough = TRUE){
  promoter.region <- promoters(genes, upstream = 500, downstream = 500)
  tss.resize <- resize(promoter.region, width = 750, fix = "end")
  upstream <- resize(promoter.region, width = 500, fix = "start")
  DATA <- DATA %>% filter(DATA$score!=0)
  counts <- getCountsByRegions(DATA,tss.resize)
  cutoff <- quantile(counts,exp.level.minimum)
  expressed <- ifelse(counts >= cutoff, TRUE, FALSE) #is gene activity greater than percentile cutoff (exp.level.min)
  expressed <- ifelse(counts >= hard.minimum, expressed, FALSE) #is gene activity greater than absolute cutoff (hard.min)
  rt <- getCountsByRegions(DATA, upstream)
  upstream.pass <- rt < counts #is upstream activity less than gene activity
  if (keep.exp == TRUE & filter.readthrough == TRUE) {
    genes <- genes %>% filter(expressed == TRUE, upstream.pass == TRUE)
  } else if (keep.exp == TRUE & filter.readthrough == FALSE) {
    genes <- genes %>% filter(expressed == TRUE)
  } else {
    genes <- genes %>% filter(expressed == FALSE)
  }
  return(genes)
}


####################################################################################################################################
# DO NOT EDIT BELOW THIS LINE
####################################################################################################################################

#no_cores = detectCores() / 8
findLowestPValue <- function(matrix, windowSize, no_cores = 12) {
  num_rows <- nrow(matrix)
  num_cols <- ncol(matrix)
  # Calculate rolling medians for each row
  roll_medians <- t(apply(matrix, 1, function(row) rollapply(row, windowSize, median, align = "left", fill = NA)))
  # Preallocate the p_values matrix with high values
  p_values <- matrix(rep(1000000000000, num_rows * (num_cols - windowSize + 1)), nrow = num_rows)
  # We can now compare medians in a vectorized way by operating on slices of the roll_medians matrix
  left_medians <- roll_medians[, 1:(num_cols - 2 * windowSize + 1)]
  right_medians <- roll_medians[, (windowSize + 1):(num_cols - windowSize + 1)]
  # Indices where left median is less than the right median
  valid_indices <- which(left_medians < right_medians, arr.ind = TRUE)
  # Calculate the single index for the valid_indices in the p_values matrix
  single_indices <- (valid_indices[, 2] - 1) * num_rows + valid_indices[, 1]
  # Perform Wilcoxon test on valid indices only
  if (length(valid_indices) > 0) {
    # Extract values for Wilcoxon test
    p_values_vector <- mclapply(seq_len(nrow(valid_indices)), function(k) {
      i <- valid_indices[k, 1]
      j <- valid_indices[k, 2]
      left_dist <- matrix[i, (j):(j + windowSize - 1)]
      right_dist <- matrix[i, (j + windowSize):(j + 2 * windowSize - 1)]
      wilcox.test(left_dist, right_dist, alternative = "less")$p.value
    }, mc.cores = no_cores)
    # Assign the computed p-values back into the p_values matrix using the single index
    p_values[single_indices] <- unlist(p_values_vector)
  }
  # Find the minimum p-value and its index for each row
  lowest_p_values <- apply(p_values, 1, min, na.rm = TRUE)
  index.lowest_p_values <- max.col(-p_values, ties.method = "first")  # Get the index of the min value
  adj.index <- index.lowest_p_values + (windowSize-1) # account for windowing
  return(list(p_values = lowest_p_values, indices = adj.index))
}
#alternative vectorization of resizeGranges, now working, make standard.
resizeGRanges <- function(gr, resizeVector, windowSize = 10, stepsize = 20) {
  if (length(resizeVector) != length(gr))
    stop("The length of resizeVector must be equal to the number of entries in GRanges.")
  if (any(is.na(resizeVector))) {
    stop("resizeVector cannot contain NA values.")
  }
  if (any(is.na(start(gr))) || any(is.na(end(gr)))) {
    stop("GRanges object cannot have NA values in 'start' or 'end'.")
  }
  # Calculate resized centers
  resizedStarts <- start(resize(gr, width = 1, fix = "start"))
  # Adjust start and end positions based on strand information
  startPositions <- ifelse(as.character(strand(gr)) == "+", resizedStarts + (resizeVector*windowSize), resizedStarts - (resizeVector*windowSize))
  endPositions <- startPositions + 1 # set endPositions <- startPositions if working with odd-number-sized ranges
  # Construct and return the resized GRanges
  output <- GRanges(seqnames = seqnames(gr), ranges = IRanges(start = startPositions, end = endPositions), strand = strand(gr))
  return(resize(output, width = 1, fix = "start"))
}
Reannotate_TSS_stepFunc_GR <- function(GRANGES, old_regions, binsize = 10, search.window = 1000, step = 20, name="Test", field = "score"){
  options(warn=-1)
  #Resize the old_regions to a search.window-sized bin centered on the currently annotated TSS (allows full gene annotations as input)
  old_regions_res <- promoters(old_regions, upstream = search.window/2, downstream = search.window/2)
  #replaced to take 1kb window centered on TSS as input:
  #old_regions_res <- resize(old_regions, width = search.window, fix = "center")
  #Count signal
  GRANGES$score <- abs(GRANGES$score);
  #Count signal
  proseq <- getCountsByPositions(GRANGES, old_regions_res, binsize = binsize, field="score", FUN = "sum") #confirm if this considers strand info- is bin 1 the same regardless of strand?
  old_regions$proseq.cts <- rowSums(proseq);
  old_regions_res$proseq.cts <- rowSums(proseq);
  TSN_index <- findLowestPValue(matrix=proseq, windowSize = step)$indices #originally called proseq2
  #Resize the old annotations to the new TSN index
  new_annotations <- resizeGRanges(gr = old_regions_res, resizeVector = TSN_index, windowSize = binsize, stepsize = step)
  new_annotations_10bp <- promoters( new_annotations, upstream = 2, downstream = binsize) # originally upstream = 2, downstream = 10
  new_annotations_10bp_cts <- getCountsByPositions(GRANGES, new_annotations_10bp, binsize = 1, field="score", FUN = "sum");
  index.max_new_annotations_10bp_cts <- vector(mode = "numeric");
  for (i in 1:nrow(new_annotations_10bp_cts)) {
    index.max_new_annotations_10bp_cts[i] <- max(which(new_annotations_10bp_cts[i,] == max(new_annotations_10bp_cts[i,])))
  }
  new_annotations_within.10bp <- resizeGRanges(gr = new_annotations_10bp, resizeVector = index.max_new_annotations_10bp_cts, windowSize = 1) #output should be single base TSS
  #dif_TSS_in.bp = (start(new_annotations_within.10bp) - start(old_regions));       
  #removed line which outputs diff bw old and new annotations to a file
  #new_annotations_within.10bp <- promoters(new_annotations_within.10bp, upstream = 150, downstream = 500) #allows output of just promoter region- can set depending on need. alex plts want -150-500
  start(new_annotations_within.10bp) <- ifelse(strand(old_regions)=="-", start(old_regions), start(new_annotations_within.10bp))
  end(new_annotations_within.10bp) <- ifelse(strand(old_regions)=="-", end(new_annotations_within.10bp), end(old_regions)) #this pair of lines adds unchanged gene ends in a stranded manner to altered start site, allowing output of full gene with new tss
  return(new_annotations_within.10bp)
}

