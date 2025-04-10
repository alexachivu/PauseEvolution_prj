####################################################################################################################################
#save.image("/local/workdir/bab364/r_envts/PauseEvoPlots15may2024.RData")
# this script loads new annotations & 3' proseq data
# generates new metaprofile plots, pause index plots, some supplemental figures
# required functions at bottom of script
load("/local/workdir/bab364/r_envts/PauseEvoPlots15may2024.RData")
####################################################################################################################################
library(Biostrings)
library(BRGenomics)
library(rtracklayer)
library(BSgenome)
library(tidyverse)
library(plyranges) #all above packages added to paperpile
# plots only
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggsignif) #all above packages added to paperpile
# reannotate only
library(zoo)
library(parallel) #part of R, already cited
# motif only
library(seqLogo) #all above packages added to paperpile

# build dataframe to store species info
figure.data <- data.frame(
  species = c("E. coli","H. mediterranei","A. thaliana","Z. mays","O. sativa","D. discoideum", "S. pombe", "S. cerevisiae","S. arctica","C. fragrantissima", "C. owczarzaki","N. vectensis","C. elegans","D. pulex","D. iulia","D. melanogaster", "S. purpuratus", "P. marinus", "M. musculus", "H. sapiens"),
  NELF = c("None","None","None","None","None","NELF-B and C/D","None","None","NELF-B and C/D","NELF-B and C/D","All, or lacking NELF-E","All, or lacking NELF-E","None","All, or lacking NELF-E","All, or lacking NELF-E","All, or lacking NELF-E","All, or lacking NELF-E","All, or lacking NELF-E","All, or lacking NELF-E","All, or lacking NELF-E"),
  no_NELF = c(0,0,0,0,0,2,0,0,2,2,3,4,0,4,4,4,4,4,4,4),
  any_NELF = c("Absent","Absent","Absent","Absent","Absent","Present","Absent","Absent","Present","Present","Present","Present","Absent","Present","Present","Present","Present","Present","Present","Present"),
  HEXIM = c("Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent","Present","Absent","Present","Present","Present","Present","Present","Present","Present"),
  PI = rep(0, 20),
  pauseButtonScore = rep(0, 20),
  InitiatorScore = rep(0, 20)
  )

# E.coli ##### 
#load proseq
ECOLI_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ecoli/Ec-merge.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
ec.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Ecoli.step.reann.rds")
#make plot
p1 <- metaProfilePauseEvolution(ec.genes, ECOLI_3p_merge_bam)
pi1 <- calcPI(ec.genes,ECOLI_3p_merge_bam, output = "individual")
#calculate pause index, add DF
figure.data$PI[1] <- calcPI(ec.genes,ECOLI_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vectors
ec.pauses <- promoters(ec.genes, upstream = 0, downstream = 100)
ec.background <- random.flank(ec.pauses)
ec.pause.seqs <- extractSequenceFromFasta(ec.pauses, "/fs/cbsubscb17/storage/data/short_read_index/ec_mg1655/escherichia_coli_mg1655_01312020.fasta", rename.chrs = T)
ec.background.seqs <- extractSequenceFromFasta(ec.background, "/fs/cbsubscb17/storage/data/short_read_index/ec_mg1655/escherichia_coli_mg1655_01312020.fasta", rename.chrs = T)
figure.data$pauseButtonScore[1] <- motifEnrichment(ec.pause.seqs, ec.background.seqs, pause.button)
ec.tss.seqs <- extractSequenceFromFasta(promoters(ec.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/ec_mg1655/escherichia_coli_mg1655_01312020.fasta", rename.chrs = T)
figure.data$InitiatorScore[1] <- motifEnrichment(ec.tss.seqs, ec.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# H. mediterranei
#load proseq
HALO_3p_merge_manual <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Hmediterranei/Hmedi-merge.sort.bam",revcomp = T, trim.to = "3p", paired_end = T)
#load new TSS
hm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Hmediterranei.step.reann.rds")
#make plot
p2 <- metaProfilePauseEvolution(hm.genes, HALO_3p_merge_manual)
pi2 <- calcPI(hm.genes,HALO_3p_merge_manual, output = "individual")
#calculate pause index, add to protected vector
figure.data$PI[2] <- calcPI(hm.genes,HALO_3p_merge_manual,output = "mean")
#calculate motif enrichment, add to protected vector
hm.pauses <- promoters(hm.genes, upstream = 0, downstream = 100)
hm.background <- random.flank(hm.pauses)
hm.pause.seqs <- extractSequenceFromFasta(hm.pauses, "/fs/cbsubscb17/storage/data/short_read_index/Haloferax_archaea/Hm.sequence.fasta", rename.chrs = T)
hm.background.seqs <- extractSequenceFromFasta(hm.background, "/fs/cbsubscb17/storage/data/short_read_index/Haloferax_archaea/Hm.sequence.fasta", rename.chrs = T)
figure.data$pauseButtonScore[2] <- motifEnrichment(hm.pause.seqs, hm.background.seqs, pause.button)
hm.tss.seqs <- extractSequenceFromFasta(promoters(hm.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/Haloferax_archaea/Hm.sequence.fasta", rename.chrs = T)
figure.data$InitiatorScore[2] <- motifEnrichment(hm.tss.seqs, hm.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# A. thaliana ##### 
#load proseq
ARA_3p <- import_bam("/fs/cbsubscb17/storage/projects/NASA_2020/Arabidopsis_Gro_1/Arabidopsis_Gro_1_QC.sort.bam",revcomp = F, trim.to = "3p")
#load new TSS
athal.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Athaliana.step.reann.rds")
#make plot
p3 <- metaProfilePauseEvolution(athal.genes, ARA_3p)
pi3 <- calcPI(athal.genes,ARA_3p, output = "individual")
# pause index
figure.data$PI[3] <- calcPI(athal.genes,ARA_3p,output = "mean")
#calculate motif enrichment, add to protected vector
athal.pauses <- promoters(athal.genes, upstream = 0, downstream = 100)
athal.background <- random.flank(athal.pauses)
athal.pause.seqs <- extractSequenceFromFasta(athal.pauses, "/fs/cbsubscb17/storage/data/short_read_index/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa", rename.chrs = T)
athal.background.seqs <- extractSequenceFromFasta(athal.background, "/fs/cbsubscb17/storage/data/short_read_index/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa", rename.chrs = T)
figure.data$pauseButtonScore[3] <- motifEnrichment(athal.pause.seqs, athal.background.seqs, pause.button)
athal.tss.seqs <- extractSequenceFromFasta(promoters(athal.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa", rename.chrs = T)
figure.data$InitiatorScore[3] <- motifEnrichment(athal.tss.seqs, athal.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# Z. mays 
#load proseq
ZMAYS_PE2 <- import_bigWig("/fs/cbsubscb17/storage/data/short_read_index/Zea.mays/Z.mays_3p.end_plus.bw","/fs/cbsubscb17/storage/data/short_read_index/Zea.mays/Z.mays_3p.end_minus.bw")
#load new TSS
zm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Zmays.step.reann.rds")
#make plot
p4 <- metaProfilePauseEvolution(zm.genes, ZMAYS_PE2)
pi4 <- calcPI(zm.genes,ZMAYS_PE2, output = "individual")
# pause index
figure.data$PI[4] <- calcPI(zm.genes,ZMAYS_PE2,output = "mean")
#calculate motif enrichment, add to protected vector
zm.pauses <- promoters(zm.genes, upstream = 0, downstream = 100)
zm.background <- random.flank(zm.pauses)
zm.pause.seqs <- extractSequenceFromFasta(zm.pauses, "/local/storage/data/short_read_index/Zea.mays/ncbi-genomes-2022-09-20/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna", rename.chrs = T)
zm.background.seqs <- extractSequenceFromFasta(zm.background, "/local/storage/data/short_read_index/Zea.mays/ncbi-genomes-2022-09-20/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[4] <- motifEnrichment(zm.pause.seqs, zm.background.seqs, pause.button)
zm.tss.seqs <- extractSequenceFromFasta(promoters(zm.genes, upstream = 50, downstream = 50), "/local/storage/data/short_read_index/Zea.mays/ncbi-genomes-2022-09-20/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna", rename.chrs = T)
figure.data$InitiatorScore[4] <- motifEnrichment(zm.tss.seqs, zm.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# O. sativa
#load proseq
OSAT_PE <- import_bigWig("/fs/cbsubscb17/storage/data/short_read_index/O.sativa/O.sativa_3p.end_plus.bw","/fs/cbsubscb17/storage/data/short_read_index/O.sativa/O.sativa_3p.end_minus.bw")
#load new TSS
os.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Osativa.step.reann.rds")
#make plot
p5 <- metaProfilePauseEvolution(os.genes, OSAT_PE)
pi5 <- calcPI(os.genes,OSAT_PE, output = "individual")
figure.data$PI[5] <- calcPI(os.genes,OSAT_PE,output = "mean")
#calculate motif enrichment, add to protected vector
os.pauses <- promoters(os.genes, upstream = 0, downstream = 100)
os.background <- random.flank(os.pauses)
os.pause.seqs <- extractSequenceFromFasta(os.pauses, "/local/storage/data/short_read_index/O.sativa/O.sativa_GCF_001433935.1_IRGSP-1.0_genomic.fna", rename.chrs = T)
os.background.seqs <- extractSequenceFromFasta(os.background, "/local/storage/data/short_read_index/O.sativa/O.sativa_GCF_001433935.1_IRGSP-1.0_genomic.fna", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[5] <- motifEnrichment(os.pause.seqs, os.background.seqs, pause.button)
os.tss.seqs <- extractSequenceFromFasta(promoters(os.genes, upstream = 50, downstream = 50), "/local/storage/data/short_read_index/O.sativa/O.sativa_GCF_001433935.1_IRGSP-1.0_genomic.fna", rename.chrs = T)
figure.data$InitiatorScore[5] <- motifEnrichment(os.tss.seqs, os.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# D. discoideum ##### some TSS are consolidated by re-annotation. these have different PAS, so they are retained in current workflow #####
#load proseq
DIC_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ddiscoideum/Ddisco-merge.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
dd.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Ddiscoideum.step.reann.rds")
#make plot
p6 <- metaProfilePauseEvolution(dd.genes, DIC_3p_merge_bam)
pi6 <- calcPI(dd.genes,DIC_3p_merge_bam, output = "individual")
figure.data$PI[6] <- calcPI(dd.genes,DIC_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
dd.pauses <- promoters(dd.genes, upstream = 0, downstream = 100)
dd.background <- random.flank(dd.pauses)
dd.pause.seqs <- extractSequenceFromFasta(dd.pauses, "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Dictyostelium_discoideum/fasta/default/Dictyostelium_discoideum.fa", rename.chrs = T)
dd.background.seqs <- extractSequenceFromFasta(dd.background, "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Dictyostelium_discoideum/fasta/default/Dictyostelium_discoideum.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[6] <- motifEnrichment(dd.pause.seqs, dd.background.seqs, pause.button)
dd.tss.seqs <- extractSequenceFromFasta(promoters(dd.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Dictyostelium_discoideum/fasta/default/Dictyostelium_discoideum.fa", rename.chrs = T)
figure.data$InitiatorScore[6] <- motifEnrichment(dd.tss.seqs, dd.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# S. pombe ##### 
#load proseq
#SPOM_5p <- import_bam("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/S.pombe/Pombe_WT_12.bam",revcomp = T, trim.to = "5p", paired_end = F)
SPOM_PE <- import_bigWig("/local/workdir/James/PauseEvolution/data/S.pombe/Pombe_WT_merged_plus.rpm.bw","/local/workdir/James/PauseEvolution/data/S.pombe/Pombe_WT_merged_minus.rpm.bw")
SPOM_cap <- import_bigWig ("/fs/cbsubscb17/storage/data/S.pombe/procap/GSM1974988_3870_7157_12385_C53ARACXX_pombe972h-ProCap_ATCACG_R1Normed_plus.bw", "/fs/cbsubscb17/storage/data/S.pombe/procap/GSM1974988_3870_7157_12385_C53ARACXX_pombe972h-ProCap_ATCACG_R1Normed_minus.bw")
seqlevelsStyle(SPOM_cap) <- seqlevelsStyle(SPOM_PE)
#load new TSS
spom.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Spombe.step.reann.rds")
#POMBE.reannotated.metaplot <- MetaProfilePlotMulti_sample_supplement(SPOM_PE, SPOM_cap, ANNOTATIONS = promoters(spom.genes, upstream = 500, downstream = 500), SPECIES = "Pombe Final")
#POMBE.consensus.metaplot <- MetaProfilePlotMulti_sample_supplement(SPOM_PE, SPOM_cap, ANNOTATIONS = promoters(sp.genes.filter, upstream = 500, downstream = 500), SPECIES = "Pombe OG")
#make plot
p7 <- metaProfilePauseEvolution(spom.genes, SPOM_PE)
pi7 <- calcPI(spom.genes,SPOM_PE, output = "individual")
figure.data$PI[7] <- calcPI(spom.genes,SPOM_PE,output = "mean")
#calculate motif enrichment, add to protected vector
spom.pauses <- promoters(spom.genes, upstream = 0, downstream = 100)
spom.background <- random.flank(spom.pauses)
spom.pause.seqs <- extractSequenceFromFasta(spom.pauses, "/local/storage/data/short_read_index/S_pombe/Schizosaccharomyces_pombe.ASM294v2.simple.chrs.fa")
spom.background.seqs <- extractSequenceFromFasta(spom.background, "/local/storage/data/short_read_index/S_pombe/Schizosaccharomyces_pombe.ASM294v2.simple.chrs.fa", BACKGROUND = T)
figure.data$pauseButtonScore[7] <- motifEnrichment(spom.pause.seqs, spom.background.seqs, pause.button)
spom.tss.seqs <- extractSequenceFromFasta(promoters(spom.genes, upstream = 50, downstream = 50), "/local/storage/data/short_read_index/S_pombe/Schizosaccharomyces_pombe.ASM294v2.simple.chrs.fa")
figure.data$InitiatorScore[7] <- motifEnrichment(spom.tss.seqs, spom.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# S. cerevisiae
#load TSS
#load proseq
SCER_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/S.cerevisiae/My_output-05_22_2023/S.cerevisiae_Gro_12_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/S.cerevisiae/My_output-05_22_2023/S.cerevisiae_Gro_12_minus.bw") #SE
SCER_PE <- import_bigWig("/local/workdir/James/PauseEvolution/data/S.cerevisiae/GSM1974983_4474_6102_15483_H3HN2BGXX_SC_WT_w303a_Combined_PROseq_GATCAG_R1spikeNormed_plus.bw","/local/workdir/James/PauseEvolution/data/S.cerevisiae/GSM1974983_4474_6102_15483_H3HN2BGXX_SC_WT_w303a_Combined_PROseq_GATCAG_R1spikeNormed_minus.bw")
seqlevelsStyle(SCER_PE) <- seqlevelsStyle(SCER_5p)
SCER_cap <- import_bigWig("/fs/cbsubscb17/storage/data/sacser/procap/GSM1974987_3870_7157_12387_C53ARACXX_cerevisiaeW303-aProCap_TTAGGC_R1Normed_plus.bw", "/fs/cbsubscb17/storage/data/sacser/procap/GSM1974987_3870_7157_12387_C53ARACXX_cerevisiaeW303-aProCap_TTAGGC_R1Normed_minus.bw")
seqlevelsStyle(SCER_cap) <- seqlevelsStyle(SCER_5p)
#load new TSS
sc.genes <- readRDS(file ="/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Scerevisiae.step.reann.rds")
#CEREV.reannotated.metaplot <- MetaProfilePlotMulti_sample_supplement(SCER_PE, SCER_cap, ANNOTATIONS = promoters(sc.genes, upstream = 500, downstream = 500), SPECIES = "Cerevisiae Final")
#CEREV.consensus.metaplot <- MetaProfilePlotMulti_sample_supplement(SCER_PE, SCER_cap, ANNOTATIONS = promoters(sc.genes.filter, upstream = 500, downstream = 500), SPECIES = "Cerevisiae OG")
#make plot
p8 <- metaProfilePauseEvolution(sc.genes, SCER_PE)
pi8 <- calcPI(sc.genes,SCER_PE, output = "individual")
figure.data$PI[8] <- calcPI(sc.genes,SCER_PE,output = "mean")
#calculate motif enrichment, add to protected vector
sc.pauses <- promoters(sc.genes, upstream = 0, downstream = 100)
sc.background <- random.flank(sc.pauses)
sc.pause.seqs <- extractSequenceFromFasta(sc.pauses, "/fs/cbsubscb17/storage/data/short_read_index/S_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa", rename.chrs = T)
sc.background.seqs <- extractSequenceFromFasta(sc.background, "/fs/cbsubscb17/storage/data/short_read_index/S_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[8] <- motifEnrichment(sc.pause.seqs, sc.background.seqs, pause.button)
sc.tss.seqs <- extractSequenceFromFasta(promoters(sc.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/S_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa", rename.chrs = T)
figure.data$InitiatorScore[8] <- motifEnrichment(sc.tss.seqs, sc.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# S. arctica #####
#load proseq
SARCT_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Sarctica/Sarctica-merge.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
sa.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Sarctica.step.reann.rds")
#make plot
p9 <- metaProfilePauseEvolution(sa.genes, SARCT_3p_merge_bam)
pi9 <- calcPI(sa.genes,SARCT_3p_merge_bam, output = "individual")
figure.data$PI[9] <- calcPI(sa.genes,SARCT_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
sa.pauses <- promoters(sa.genes, upstream = 0, downstream = 100)
sa.background <- random.flank(sa.pauses)
sa.pause.seqs <- extractSequenceFromFasta(sa.pauses, "/fs/cbsubscb17/storage/data/short_read_index/S_arctica/Sphaeroforma_arctica_jp610.Spha_arctica_JP610_V1.dna.toplevel.fa", rename.chrs = T)
sa.background.seqs <- extractSequenceFromFasta(sa.background, "/fs/cbsubscb17/storage/data/short_read_index/S_arctica/Sphaeroforma_arctica_jp610.Spha_arctica_JP610_V1.dna.toplevel.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[9] <- motifEnrichment(sa.pause.seqs, sa.background.seqs, pause.button) #likely issue: sequences of entirely Ns
sa.tss.seqs <- extractSequenceFromFasta(promoters(sa.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/S_arctica/Sphaeroforma_arctica_jp610.Spha_arctica_JP610_V1.dna.toplevel.fa", rename.chrs = T)
figure.data$InitiatorScore[9] <- motifEnrichment(sa.tss.seqs, sa.background.seqs, initiator)
# this genome assembly contains many short contigs which cause issues with background region selection. 
# It also contains a high N content, so some selected background regions are partially or completely "N"
# originally 10690 sequences. filtered out contigs below 2kb: leaving 9531 sequences
# background regions containing fewer than 21 sequential non-N bases are ignored
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# C. fragrantissima ##### 
#load proseq
CRE_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Cfragrantissima/Cfragrantissima-merge.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
cf.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Cfragrantissima.step.reann.rds")
#make plot
p10 <- metaProfilePauseEvolution(cf.genes, CRE_3p_merge_bam)
pi10 <- calcPI(cf.genes,CRE_3p_merge_bam, output = "individual")
figure.data$PI[10] <- calcPI(cf.genes,CRE_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
cf.pauses <- promoters(cf.genes, upstream = 0, downstream = 100)
cf.background <- random.flank(cf.pauses)
cf.pause.seqs <- extractSequenceFromFasta(cf.pauses, "/fs/cbsubscb17/storage/data/short_read_index/C_fragrantissima/Creolimax_fragrantissima.genome.fasta", rename.chrs = T)
cf.background.seqs <- extractSequenceFromFasta(cf.background, "/fs/cbsubscb17/storage/data/short_read_index/C_fragrantissima/Creolimax_fragrantissima.genome.fasta", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[10] <- motifEnrichment(cf.pause.seqs, cf.background.seqs, pause.button)
cf.tss.seqs <- extractSequenceFromFasta(promoters(cf.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/C_fragrantissima/Creolimax_fragrantissima.genome.fasta", rename.chrs = T)
figure.data$InitiatorScore[10] <- motifEnrichment(cf.tss.seqs, cf.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# C. owczarzaki
#load proseq
CAP_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Cowczarzaki/Cowczarzaki-merge.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
co.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Cowczarzaki.step.reann.rds")
#make plot
p11 <- metaProfilePauseEvolution(co.genes, CAP_3p_merge_bam)
pi11 <- calcPI(co.genes,CAP_3p_merge_bam, output = "individual")
figure.data$PI[11] <- calcPI(co.genes,CAP_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
co.pauses <- promoters(co.genes, upstream = 0, downstream = 100)
co.background <- random.flank(co.pauses)
co.pause.seqs <- extractSequenceFromFasta(co.pauses, "/fs/cbsubscb17/storage/data/short_read_index/C_owczarzaki/Capsaspora_owczarzaki_atcc_30864.C_owczarzaki_V2.dna.toplevel.fa", rename.chrs = T)
co.background.seqs <- extractSequenceFromFasta(co.background, "/fs/cbsubscb17/storage/data/short_read_index/C_owczarzaki/Capsaspora_owczarzaki_atcc_30864.C_owczarzaki_V2.dna.toplevel.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[11] <- motifEnrichment(co.pause.seqs, co.background.seqs, pause.button)
co.tss.seqs <- extractSequenceFromFasta(promoters(co.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/C_owczarzaki/Capsaspora_owczarzaki_atcc_30864.C_owczarzaki_V2.dna.toplevel.fa", rename.chrs = T)
figure.data$InitiatorScore[11] <- motifEnrichment(co.tss.seqs, co.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# N. vectensis ##### 
#load proseq
NEM_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Nvectensis/Nvectensis-merge.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
nv.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Nvectensis.step.reann.rds")
#make plot
p12 <- metaProfilePauseEvolution(nv.genes, NEM_3p_merge_bam)
pi12 <- calcPI(nv.genes,NEM_3p_merge_bam, output = "individual")
figure.data$PI[12] <- calcPI(nv.genes,NEM_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
nv.pauses <- promoters(nv.genes, upstream = 0, downstream = 100)
nv.background <- random.flank(nv.pauses)
nv.pause.seqs <- extractSequenceFromFasta(nv.pauses, "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Nematostella_vectensis/fasta/default/Nematostella_vectensis.fa", rename.chrs = T)
nv.background.seqs <- extractSequenceFromFasta(nv.background, "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Nematostella_vectensis/fasta/default/Nematostella_vectensis.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[12] <- motifEnrichment(nv.pause.seqs, nv.background.seqs, pause.button)
nv.tss.seqs <- extractSequenceFromFasta(promoters(nv.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Nematostella_vectensis/fasta/default/Nematostella_vectensis.fa", rename.chrs = T)
figure.data$InitiatorScore[12] <- motifEnrichment(nv.tss.seqs, nv.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# C. elegans ##### many whole genes are collapsed by reannotation. These have been de-duplicated. Some TSS pileups remain #####
#load proseq # RERUN THIS BLOCK INCLUDING seqlevelStyle() BEFORE GENERATING FIGURE S1
CEL_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/C.elegans/Celegans_CHROMSmatchGFF_5p_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/C.elegans/Celegans_CHROMSmatchGFF_5p_minus.bw") #SE
CEL_PE_G <- import_bigWig("/local/workdir/James/PauseEvolution/data/C.elegans/Celegans_Gro_merged_rpm.plus.bw","/local/workdir/James/PauseEvolution/data/C.elegans/Celegans_Gro_merged_rpm.minus.bw",genome = "ce11")
CEL_PE_G <- CEL_PE_G %>% filter(CEL_PE_G$score!=0) #corrects null value issue
seqlevelsStyle(CEL_PE_G) <- seqlevelsStyle(CEL_5p)
#load new TSS # DO NOT RERUN REMAINING CODE IN BLOCK
ce.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Celegans.step.reann.rds")
#make plot
p13 <- metaProfilePauseEvolution(ce.genes, CEL_PE_G)
pi13 <- calcPI(ce.genes,CEL_PE_G, output = "individual") # may need to move line to avoid seqlevel mismatch
figure.data$PI[13] <- calcPI(ce.genes,CEL_PE_G,output = "mean")
#calculate motif enrichment, add to protected vector
ce.pauses <- promoters(ce.genes, upstream = 0, downstream = 100)
ce.background <- random.flank(ce.pauses)
CEL_PE_G <- import_bigWig("/local/workdir/James/PauseEvolution/data/C.elegans/Celegans_Gro_merged_rpm.plus.bw","/local/workdir/James/PauseEvolution/data/C.elegans/Celegans_Gro_merged_rpm.minus.bw",genome = "ce11")
CEL_PE_G <- CEL_PE_G %>% filter(CEL_PE_G$score!=0) #corrects null value issue
seqlevelsStyle(ce.pauses) <- seqlevelsStyle(CEL_PE_G) #corrects seqlevel missmatch
seqlevelsStyle(ce.background) <- seqlevelsStyle(CEL_PE_G) #corrects seqlevel missmatch
ce.pause.seqs <- extractSequenceFromFasta(ce.pauses, "/fs/cbsubscb17/storage/data/short_read_index/C_elegans/ce6.fa", rename.chrs = T)
ce.background.seqs <- extractSequenceFromFasta(ce.background, "/fs/cbsubscb17/storage/data/short_read_index/C_elegans/ce6.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[13] <- motifEnrichment(ce.pause.seqs, ce.background.seqs, pause.button)
ce.tss.seqs <- extractSequenceFromFasta(promoters(ce.pauses, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/C_elegans/ce6.fa", rename.chrs = T)
figure.data$InitiatorScore[13] <- motifEnrichment(ce.tss.seqs, ce.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# D. pulex ##### 
#load proseq
DAPH_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Dpulex/Dpulex-merge.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
dp.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Dpulex.step.reann.rds")
#make plot
p14 <- metaProfilePauseEvolution(dp.genes, DAPH_3p_merge_bam)
pi14 <- calcPI(dp.genes,DAPH_3p_merge_bam, output = "individual")
figure.data$PI[14] <- calcPI(dp.genes,DAPH_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
dp.pauses <- promoters(dp.genes, upstream = 0, downstream = 100)
dp.background <- random.flank(dp.pauses)
dp.pause.seqs <- extractSequenceFromFasta(dp.pauses, "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Daphnia/fasta/default/Daphnia.fa", rename.chrs = T)
dp.background.seqs <- extractSequenceFromFasta(dp.background, "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Daphnia/fasta/default/Daphnia.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[14] <- motifEnrichment(dp.pause.seqs, dp.background.seqs, pause.button)
dp.tss.seqs <- extractSequenceFromFasta(promoters(dp.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Daphnia/fasta/default/Daphnia.fa", rename.chrs = T)
figure.data$InitiatorScore[14] <- motifEnrichment(dp.tss.seqs, dp.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# D. iulia ##### 
#load proseq
DIU_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Diulia/Diulia-merge.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
di.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Diulia.step.reann.rds")
#make plot
p15 <- metaProfilePauseEvolution(di.genes, DIU_3p_merge_bam)
pi15 <- calcPI(di.genes,DIU_3p_merge_bam, output = "individual")
figure.data$PI[15] <- calcPI(di.genes,DIU_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
di.pauses <- promoters(di.genes, upstream = 0, downstream = 100)
di.background <- random.flank(di.pauses)
di.pause.seqs <- extractSequenceFromFasta(di.pauses, "/fs/cbsubscb17/storage/projects/D.iulia/D.iulia_v1.1/D.iulia_v1.1.fa", rename.chrs = T)
di.background.seqs <- extractSequenceFromFasta(di.background, "/fs/cbsubscb17/storage/projects/D.iulia/D.iulia_v1.1/D.iulia_v1.1.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[15] <- motifEnrichment(di.pause.seqs, di.background.seqs, pause.button)
di.tss.seqs <- extractSequenceFromFasta(promoters(di.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/projects/D.iulia/D.iulia_v1.1/D.iulia_v1.1.fa", rename.chrs = T)
figure.data$InitiatorScore[15] <- motifEnrichment(di.tss.seqs, di.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# D. melanogaster
#load proseq
#DROS_5p <- import_bigWig("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/D.melanogaster/My_output-05_25_2023/Drosophila_Gro_QC.sort_plus.bw", "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/SE_data/D.melanogaster/My_output-05_25_2023/Drosophila_Gro_QC.sort_minus.bw") #SE
DROS_A1 <- import_bigWig("/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_1/Drosophilia_01A_1_QC_plus.rpm.bw","/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_1/Drosophilia_01A_1_QC_minus.rpm.bw", genome = "dm3")
DROS_A2 <- import_bigWig("/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_2/Drosophilia_01A_2_QC_plus.rpm.bw","/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_2/Drosophilia_01A_2_QC_minus.rpm.bw", genome = "dm3")
DROS_A3 <- import_bigWig("/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_3/Drosophilia_01A_3_QC_plus.rpm.bw","/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_3/Drosophilia_01A_3_QC_minus.rpm.bw", genome = "dm3")
DROS_A4 <- import_bigWig("/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_4/Drosophilia_01A_4_QC_plus.rpm.bw","/fs/cbsubscb17/storage/projects/NASA_2020/Drosophilia_01A_4/Drosophilia_01A_4_QC_minus.rpm.bw", genome = "dm3")
DROS_PE <- mergeGRangesData(DROS_A1,DROS_A2,DROS_A3,DROS_A4)
rm(DROS_A1,DROS_A2,DROS_A3,DROS_A4)
DROS_cap <- import_bigWig("/fs/cbsubscb17/storage/data/dm3/s2/procap/PROcap.pl.bw", "/fs/cbsubscb17/storage/data/dm3/s2/procap/PROcap.mn.bw")
#load new TSS
dm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Dmelanogaster.step.reann.rds")
#DROS.reannotated.metaplot <- MetaProfilePlotMulti_sample_supplement(DROS_PE, DROS_cap, ANNOTATIONS = promoters(dm.genes, upstream = 500, downstream = 500), SPECIES = "Drosophila Final")
#DROS.consensus.metaplot <- MetaProfilePlotMulti_sample_supplement(DROS_PE, DROS_cap, ANNOTATIONS = promoters(dmel.genes.filter, upstream = 500, downstream = 500), SPECIES = "Drosophila OG")
#make plot
p16 <- metaProfilePauseEvolution(dm.genes, DROS_PE)
pi16 <- calcPI(dm.genes,DROS_PE, output = "individual")
figure.data$PI[16] <- calcPI(dm.genes,DROS_PE,output = "mean")
#calculate motif enrichment, add to protected vector
dm.pauses <- promoters(dm.genes, upstream = 0, downstream = 100)
dm.background <- random.flank(dm.pauses)
dm.pause.seqs <- extractSequenceFromFasta(dm.pauses, "/fs/cbsubscb17/storage/data/short_read_index/dm3/dm3.fa", rename.chrs = T)
dm.background.seqs <- extractSequenceFromFasta(dm.background, "/fs/cbsubscb17/storage/data/short_read_index/dm3/dm3.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[16] <- motifEnrichment(dm.pause.seqs, dm.background.seqs, pause.button)
dm.tss.seqs <- extractSequenceFromFasta(promoters(dm.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/dm3/dm3.fa", rename.chrs = T)
figure.data$InitiatorScore[16] <- motifEnrichment(dm.tss.seqs, dm.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# S. purpuratus # updated to use Spur 5.0 genome build # reminder to update citation to match
#load proseqros
URCH_3p_merge_bam <- import_bam("/fs/cbsubscb17/storage/data/short_read_index/sp5.0/mapped_PROseq/merged/Sea_Urchin_Embryos_12_20.merge.sort.bam",revcomp = T, trim.to = "3p", paired_end = T)
#load new TSS
sp.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Spurpuratus5.step.reann.rds")
#make plot
p17 <- metaProfilePauseEvolution(sp.genes, URCH_3p_merge_bam)
pi17 <- calcPI(sp.genes,URCH_3p_merge_bam, output = "individual")
figure.data$PI[17] <- calcPI(sp.genes,URCH_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
sp.pauses <- promoters(sp.genes, upstream = 0, downstream = 100)
sp.background <- random.flank(sp.pauses)
sp.pause.seqs <- extractSequenceFromFasta(sp.pauses, "/fs/cbsubscb17/storage/data/short_read_index/sp5.0/GCF_000002235.5/GCF_000002235.5_Spur_5.0_genomic.fa.gz", rename.chrs = T)
sp.background.seqs <- extractSequenceFromFasta(sp.background, "/fs/cbsubscb17/storage/data/short_read_index/sp5.0/GCF_000002235.5/GCF_000002235.5_Spur_5.0_genomic.fa.gz", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[17] <- motifEnrichment(sp.pause.seqs, sp.background.seqs, pause.button)
sp.tss.seqs <- extractSequenceFromFasta(promoters(sp.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/sp5.0/GCF_000002235.5/GCF_000002235.5_Spur_5.0_genomic.fa.gz", rename.chrs = T)
figure.data$InitiatorScore[17] <- motifEnrichment(sp.tss.seqs, sp.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# P. marinus
#load proseq
LAM_3p_merge_bam <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Pmarinus/Pmarinus-merge.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
#load new TSS
pm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Pmarinus.step.reann.rds")
#make plot
p18 <- metaProfilePauseEvolution(pm.genes, LAM_3p_merge_bam)
pi18 <- calcPI(pm.genes,LAM_3p_merge_bam, output = "individual")
figure.data$PI[18] <- calcPI(pm.genes,LAM_3p_merge_bam,output = "mean")
#calculate motif enrichment, add to protected vector
pm.pauses <- promoters(pm.genes, upstream = 0, downstream = 100)
pm.background <- random.flank(pm.pauses)
pm.pause.seqs <- extractSequenceFromFasta(pm.pauses, "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Lampetra/fasta/default/Lampetra.fa", rename.chrs = T)
pm.background.seqs <- extractSequenceFromFasta(pm.background, "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Lampetra/fasta/default/Lampetra.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[18] <- motifEnrichment(pm.pause.seqs, pm.background.seqs, pause.button)
pm.tss.seqs <- extractSequenceFromFasta(promoters(pm.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/bowtie2Genomes/Lampetra/fasta/default/Lampetra.fa", rename.chrs = T)
figure.data$InitiatorScore[18] <- motifEnrichment(pm.tss.seqs, pm.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# M. musculus
#load proseq ## PE2 and 3p are the same dataset with and without normalization
MOUSE_PE2  <- import_bigWig("/local/workdir/James/PauseEvolution/data/mouse_mESC/Abood123.192021.RE.F.norm_plus.bw","/local/workdir/James/PauseEvolution/data/mouse_mESC/Abood123.192021.RE.F.norm_minus.bw")
#load new TSS
mm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Mmusculus.step.reann.rds")
#make plot
p19 <- metaProfilePauseEvolution(mm.genes, MOUSE_PE2)
pi19 <- calcPI(mm.genes,MOUSE_PE2, output = "individual")
#clear envt
figure.data$PI[19] <- calcPI(mm.genes,MOUSE_PE2,output = "mean")
#calculate motif enrichment, add to protected vector
mm.pauses <- promoters(mm.genes, upstream = 0, downstream = 100)
mm.background <- random.flank(mm.pauses)
mm.pause.seqs <- extractSequenceFromFasta(mm.pauses, "/fs/cbsubscb17/storage/data/short_read_index/mm10/mm10.fa.gz", rename.chrs = T)
mm.background.seqs <- extractSequenceFromFasta(mm.background, "/fs/cbsubscb17/storage/data/short_read_index/mm10/mm10.fa.gz", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[19] <- motifEnrichment(mm.pause.seqs, mm.background.seqs, pause.button)
mm.tss.seqs <- extractSequenceFromFasta(promoters(mm.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/mm10/mm10.fa.gz", rename.chrs = T)
figure.data$InitiatorScore[19] <- motifEnrichment(mm.tss.seqs, mm.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

# H. sapiens (using ChROseq data, not Core/Martins data)
#load proseq
HUM_PE <- import_bam("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/H.sapiens/ChROseq_merged_0h.bam",revcomp = F, trim.to = "3p", paired_end = T)
HUM_cap <- import_bigWig("/local/workdir/bab364/tmp/procap/K562_grocap_plus.bigWig","/local/workdir/bab364/tmp/procap/K562_grocap_minus.bigWig", genome = "hg19")
#load new TSS
hs.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Hsapiens.step.reann.rds")
#HUM.reannotated.metaplot <- MetaProfilePlotMulti_sample_supplement(HUM_PE, HUM_cap, ANNOTATIONS = promoters(hs.genes, upstream = 500, downstream = 500), SPECIES = "Human Final")
#HUM.consensus.metaplot <- MetaProfilePlotMulti_sample_supplement(HUM_PE, HUM_cap, ANNOTATIONS = promoters(hum.genes.filter, upstream = 500, downstream = 500), SPECIES = "Human OG")
#make plot
p20 <- metaProfilePauseEvolution(hs.genes, HUM_PE, last.plot=T)
pi20 <- calcPI(hs.genes,HUM_PE, output = "individual")
#clear envt
figure.data$PI[20] <- calcPI(hs.genes,HUM_PE,output = "mean")
#calculate motif enrichment, add to protected vector
hs.pauses <- promoters(hs.genes, upstream = 0, downstream = 100)
hs.background <- random.flank(hs.pauses)
hs.pause.seqs <- extractSequenceFromFasta(hs.pauses, "/fs/cbsubscb17/storage/data/short_read_index/hg19/hg19.fa", rename.chrs = T)
hs.background.seqs <- extractSequenceFromFasta(hs.background, "/fs/cbsubscb17/storage/data/short_read_index/hg19/hg19.fa", rename.chrs = T, BACKGROUND = T)
figure.data$pauseButtonScore[20] <- motifEnrichment(hs.pause.seqs, hs.background.seqs, pause.button)
hs.tss.seqs <- extractSequenceFromFasta(promoters(hs.genes, upstream = 50, downstream = 50), "/fs/cbsubscb17/storage/data/short_read_index/hg19/hg19.fa", rename.chrs = T)
figure.data$InitiatorScore[20] <- motifEnrichment(hs.tss.seqs, hs.background.seqs, initiator)
#clear envt, run after each block to keep workspace clean
rm(list = setdiff(ls(), c(ls(pattern = "^p\\d+$"),ls(pattern = "^pi\\d+$"),"figure.data", "maxMatch","initiator", "pause.button", lsf.str())))

####################################################################################################################################
# make plots
####################################################################################################################################

# Fig 1c
figure1c <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20, align = "v", ncol=1, rel_heights = c(rep(1,19),2))
ggsave2("/local/workdir/bab364/PauseEvolution/15may2024/fig1c.pdf", plot = figure1c, device = "pdf", width = 1.75, height = 8.2, units = "in") #export to match size of OG figures, allow direct substitution

# Fig 1e
figure.data$NELF <- factor(figure.data$NELF, levels = c("None", "NELF-B and C/D", "All, or lacking NELF-E"))
figure1f <- ggplot(figure.data, aes(x = NELF, y = PI)) +
  geom_boxplot(aes(fill = NELF), alpha = .5) +
  geom_point(aes(color = NELF), alpha = 0.7) +
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),
        legend.position = "right", 
        legend.text = element_text( color="black", size=10),
        axis.text.y = element_text( color="black", size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +  
  labs(y = "PI values")+
  geom_signif(comparisons = list(c("None", "NELF-B and C/D"), c("NELF-B and C/D", "All, or lacking NELF-E"), c("None", "All, or lacking NELF-E")), y_position = c(55,65,75), tip_length = 0)+
  coord_cartesian(ylim=c(0, 85))
ggsave2("/local/workdir/bab364/PauseEvolution/15may2024/fig1f.pdf", plot = figure1f, device = "pdf", width = 4.4, height = 2.4, units = "in") #export to match size of OG figures, allow direct substitution

# Fig 2d
figure2f <- ggplot(figure.data, aes(x=log10(PI), y=pauseButtonScore)) +
  theme_classic() + 
  geom_point(aes(color=NELF)) +
  geom_smooth(method=lm , color="#E7298A", fill="#666666", se=TRUE, alpha = .1) +
  scale_color_brewer(palette = "Dark2")
summary(lm(formula = log(PI) ~ pauseButtonScore, data = figure.data))
cor.test(figure.data$pauseButtonScore, log(figure.data$PI), method = "pearson")
pdf("/local/workdir/bab364/PauseEvolution/15may2024/fig2f_5nov2024.pdf",width = 4.5, height = 2); figure2f; dev.off()

# Fig S4
piPlots <- list()
for (i in c(1:20)) {
  obj_name <- paste0("pi", i)
  # Check if the object exists in the workspace
  if (exists(obj_name)) {
    # Create the histogram for the current object
    plt_data <- data.frame(PI = get(obj_name))
    quartiles <- quantile(log10(plt_data$PI), probs = c(0.25, 0.5, 0.75))
    density_pi <- density(log10(plt_data$PI))
    dens.df <- data.frame(x = density_pi$x,y = density_pi$y,quartile = factor(findInterval(density_pi$x, quartiles)))
    if (i<20) {
      piPlots[[i]] <- ggplot(dens.df)+
        geom_line(aes(x=x,y=y))+
        geom_ribbon(aes(x=x,ymin=0,ymax=y,fill=quartile),alpha=.2)+
        scale_fill_manual(name = "",values = c("green","yellow","orange","red"))+
        geom_vline(xintercept =0, linetype="dashed")+
        theme_classic()+
        coord_cartesian(xlim=c(-5, 6))+ ##################################### originally -3 to 6 ############################################
      theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank(),plot.margin = unit(c(0, .5, 0, 4), "cm"))
    } else {
      piPlots[[i]] <- ggplot(dens.df)+
        geom_line(aes(x=x,y=y))+
        geom_ribbon(aes(x=x,ymin=0,ymax=y,fill=quartile),alpha=.2)+
        scale_fill_manual(name = "",values = c("green","yellow","orange","red"))+
        geom_vline(xintercept =0, linetype="dashed")+
        theme_classic()+
        coord_cartesian(xlim=c(-5, 6))+
        theme(legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y = element_blank(),plot.margin = unit(c(0, .5, 0, 4), "cm"))+
        labs(x="log10(Pause Index)")
    }
  } else {
    warning(paste("Object", obj_name, "does not exist"))
  }
}
labels <- c("E. coli","H. mediterranei","A. thaliana","Z. mays","O. sativa","D. discoideum", "S. pombe", "S. cerevisiae","S. arctica","C. fragrantissima", "C. owczarzaki","N. vectensis","C. elegans","D. pulex","D. iulia","D. melanogaster", "S. purpuratus", "P. marinus", "M. musculus", "H. sapiens")
FigS1 <- plot_grid(plotlist = piPlots, align = "v",axis = "1", ncol=1, rel_heights = c(rep(1,19),2), labels = as.list(labels), hjust = 0, vjust = 3)
FigS1
pdf("/local/workdir/bab364/PauseEvolution/PI.curves.06nov2024_wide.pdf",width = 5, height = 10); FigS1; dev.off()


# Fig S9
figure.data$no_NELF<- factor(figure.data$no_NELF, levels = c(0,2,3,4))
s4b <- ggplot(figure.data, aes(x = HEXIM, y = PI)) +
  geom_boxplot(aes(fill = HEXIM), alpha = .5) +
  geom_point(aes(color = HEXIM), alpha = 0.7) +
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),
        legend.position = "right", 
        legend.text = element_text( color="black", size=10),
        axis.text.y = element_text( color="black", size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +  
  labs(y = "PI values")+
  geom_signif(comparisons = list(c("Absent","Present")), y_position = 95, tip_length = 0)+
  coord_cartesian(ylim=c(0, 100))

s4a <- ggplot(figure.data, aes(x = any_NELF, y = PI)) +
  geom_boxplot(aes(fill = any_NELF), alpha = .5) +
  geom_point(aes(color = any_NELF), alpha = 0.7) +
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),
        legend.position = "right", 
        legend.text = element_text( color="black", size=10),
        axis.text.y = element_text( color="black", size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +  
  labs(y = "PI values")+
  geom_signif(comparisons = list(c("Absent","Present")), y_position = 95, tip_length = 0)+
  coord_cartesian(ylim=c(0, 100))
FigS4 <- plot_grid(plotlist = list(s4a, s4b), align = "h", ncol=2, labels = "AUTO")
ggsave2("PauseEvolution/15may2024/figS4.pdf", plot = FigS4, device = "pdf", width = 7.5, height = 3, units = "in")

# Fig s10
figS5a <- ggplot(figure.data, aes(x = NELF, y = pauseButtonScore)) +
  geom_boxplot(aes(fill = NELF), alpha = .5) +
  geom_point(aes(color = NELF), alpha = 0.7) +
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),
        legend.position = "none", 
        legend.text = element_text( color="black", size=10),
        axis.text.y = element_text( color="black", size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = unit(c(0, .5, 0, 1), "cm")) +  
  labs(y = "Pause Button Motif Enrichment")+
  geom_signif(comparisons = list(c("None", "NELF-B and C/D"), c("NELF-B and C/D", "All, or lacking NELF-E"), c("None", "All, or lacking NELF-E")), y_position = c(1.33,1.36,1.39), tip_length = 0)+
  coord_cartesian(ylim=c(1.17, 1.41))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")

figS5c <- ggplot(figure.data, aes(x=log10(PI), y=InitiatorScore)) +
  theme_classic() + 
  geom_point(aes(color=NELF)) +  
  labs(y = "Initiator Motif Enrichment")+
  geom_smooth(method=lm , color="#E7298A", fill="#666666", se=TRUE, alpha = .1) +
  scale_color_brewer(palette = "Dark2") +
  theme(plot.margin = unit(c(.5, .5, 0, 1), "cm"))
summary(lm(formula = log(PI) ~ InitiatorScore, data = figure.data))
cor.test(figure.data$InitiatorScore, figure.data$PI, method = "pearson")

figS5b <- ggplot(figure.data, aes(x = NELF, y = InitiatorScore)) +
  geom_boxplot(aes(fill = NELF), alpha = .5) +
  geom_point(aes(color = NELF), alpha = 0.7) +
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),
        legend.position = "none", 
        legend.text = element_text( color="black", size=10),
        axis.text.y = element_text( color="black", size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), plot.margin = unit(c(0, .5, 0, 1), "cm")) +  
  labs(y = "Initiator Motif Enrichment")+
  geom_signif(comparisons = list(c("None", "NELF-B and C/D"), c("NELF-B and C/D", "All, or lacking NELF-E"), c("None", "All, or lacking NELF-E")), y_position = c(1.75,1.77,1.79), tip_length = 0)+
  coord_cartesian(ylim=c(1.55, 1.8))+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")
top.row <- plot_grid(figS5a,figS5b, align = "h", labels = c("A","B"))
FigS5 <- plot_grid(top.row, figS5c, ncol = 1, align = "v", axis = 1, labels = c("","C"))
ggsave2("/local/workdir/bab364/PauseEvolution/15may2024/figS5_05nov2024.pdf", plot = FigS5, device = "pdf", width = 6, height = 6, units = "in")

# assay type supplemental figures
figure.data
figure.data$Assay <- "A"
figure.data$Assay[c(1,2,4,5,7,8,11,12,15,16,19)] <- "PRO-seq"
figure.data$Assay[c(3,13)] <- "GRO-seq"
figure.data$Assay[c(6,9,10,14,17,18,20)] <- "ChRO-seq"

# Now figure s2
figure.data$NELF <- factor(figure.data$NELF, levels = c("None", "NELF-B and C/D", "All, or lacking NELF-E"))
FigS16_2 <- ggplot(figure.data, aes(x = NELF, y = PI, fill = Assay)) +
  geom_boxplot(aes(colour = Assay),alpha = .5, outlier.alpha = .7) + 
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),legend.position = "right", legend.text = element_text( color="black", size=10),axis.text.y = element_text( color="black", size=10),axis.text.x = element_text( color="black", size=10),axis.ticks.x = element_blank()) +  
  labs(y = "Pausing Index", x = "NELF subunits")+
  coord_cartesian(ylim=c(0, 100))+ geom_signif(comparisons = list(c("ChRO-seq","Pro-seq")), y_position = 95, tip_length = 0)
FigS16_2
ggsave2("/local/workdir/bab364/PauseEvolution/FigureS2_31mar2025.pdf", plot = FigS16_2, device = "pdf", width = 6, height = 3, units = "in")
####################################################################################################################################
# create needed functions
####################################################################################################################################

metaProfilePauseEvolution <- function(ANNOTATIONS, PROseq, binsize = 10, ncores = 10, last.plot = FALSE){ #smoothness improves with binsize = 5, but runs slower
  ANNOTATIONS <- promoters(ANNOTATIONS, upstream =  150, downstream = 500)
  means <- metaSubsample(PROseq, ANNOTATIONS, binsize = binsize, prop.sample = 0.2, lower = 0.25, upper = 0.75, ncores = ncores)
  scalefactor <- round(max(means$upper),6) # this establishes an upper bound for plotting
  displaymax <- ifelse(scalefactor >= 10, round(max(means$upper), 1), round(max(means$upper), 2)) #this establishes a label to indicate scale
  print(paste0("actual max =", scalefactor))
  print(paste0("labeled max =", displaymax))
  if (last.plot == F){ #all plots except last
    PLOT <- ggplot(means, aes(x, mean, color = sample.name)) +
      geom_line(lwd=0.7, alpha=0.8) + theme_classic() + theme(panel.border = element_blank(),  
                                                              # Remove panel background
                                                              panel.background = element_blank(),
                                                              legend.background = element_blank(),
                                                              # Add axis line
                                                              axis.line = element_line(colour = "black", size=0.5),
                                                              panel.grid.major = element_line(colour = "white"),
                                                              axis.text.x = element_blank(),
                                                              axis.text.y = element_text( color="black", size=10, vjust = 1),
                                                              axis.title.x = element_blank(),
                                                              axis.title.y = element_blank(),
                                                              plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      scale_color_manual(values = c("red"))+
      geom_ribbon(aes(x, ymin = lower, ymax = upper, fill = sample.name, color = NULL), alpha = 0.2)+ # add/remove aes(color = NULL), remove size = 0.2 to change plot look if desired.
      labs(x = "Distance from TSS") + theme(legend.position = "none")+ 
      scale_x_continuous(labels = function(x) (x)-150, breaks = c(0, 150, 650))+
      scale_y_continuous(breaks = c(displaymax))+
      coord_cartesian(ylim=c(0, max(c(scalefactor,displaymax))))
  }else{ #set for last plot only, includes x axis markers
    PLOT <- ggplot(means, aes(x, mean, color = sample.name)) +
      geom_line(lwd=0.7, alpha=0.8) + theme_classic() + theme(panel.border = element_blank(),  
                                                              # Remove panel background
                                                              panel.background = element_blank(),
                                                              legend.background = element_blank(),
                                                              # Add axis line
                                                              axis.line = element_line(colour = "black", size=0.5),
                                                              panel.grid.major = element_line(colour = "white"),
                                                              axis.text.x = element_text( color="black", size=10),
                                                              axis.text.y = element_text( color="black", size=10, vjust = 1),
                                                              axis.title.x = element_text(color="black", size=10),
                                                              axis.title.y = element_blank(),
                                                              plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      scale_color_manual(values = c("red"))+
      geom_ribbon(aes(x, ymin = lower, ymax = upper, fill = sample.name, color = NULL), alpha = 0.2)+
      labs(x = "Distance from TSS") + theme(legend.position = "none")+ 
      scale_x_continuous(labels = function(x) (x)-150, breaks = c(0, 150, 650))+
      scale_y_continuous(breaks = c(displaymax))+
      coord_cartesian(ylim=c(0, max(c(scalefactor,displaymax))))
  }
  return(PLOT)
}

# Use of PRO/GRO/CHRO-seq:
figure.data
figure.data$Assay <- "Assay"
figure.data$Assay[c(4:7,19)] <- "PRO-seq"
figure.data$Assay[c(3,8,16)] <- "GRO-seq"
figure.data$Assay[c(1,2,11,12,10,14,15,13,17,18,20)] <- "ChRO-seq"
figure.data$species[figure.data$Assay == "PRO-seq"]
figure.data$species[figure.data$Assay == "GRO-seq"]
figure.data$species[figure.data$Assay == "ChRO-seq"]

proseq.only <- figure.data[c(4:7,19),] #proseq only
groseq.only <- figure.data[c(3,8,16),] #groseq only
chroseq.only <- figure.data[c(1,2,11,12,10,14,15,13,17,18,20),] #chroseq only
proseq.subunits <- ggplot(proseq.only, aes(x = any_NELF, y = PI)) +
  geom_boxplot(aes(fill = any_NELF), alpha = .5) + geom_point(aes(color = any_NELF), alpha = 0.7) +
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),legend.position = "none", legend.text = element_text( color="black", size=10),axis.text.y = element_text( color="black", size=10),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank()) +  
  labs(y = "PI values")+geom_signif(comparisons = list(c("Absent","Present")), y_position = 10, tip_length = 0)+
  coord_cartesian(ylim=c(0, 11))
groseq.subunits <- ggplot(groseq.only, aes(x = any_NELF, y = PI)) +
  geom_boxplot(aes(fill = any_NELF), alpha = .5) + geom_point(aes(color = any_NELF), alpha = 0.7) +
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),legend.position = "none", legend.text = element_text( color="black", size=10),axis.text.y = element_text( color="black", size=10),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank()) +  
  labs(y = "PI values")+ geom_signif(comparisons = list(c("Absent","Present")), y_position = 70, tip_length = 0)+
  coord_cartesian(ylim=c(0, 75))
chroseq.subunits <- ggplot(chroseq.only, aes(x = any_NELF, y = PI)) +
  geom_boxplot(aes(fill = any_NELF), alpha = .5) + geom_point(aes(color = any_NELF), alpha = 0.7) +
  theme_classic() +                                                  
  theme(axis.line = element_line(color = "black"),legend.position = "right", legend.text = element_text( color="black", size=10),axis.text.y = element_text( color="black", size=10),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank()) +  
  labs(y = "PI values")+ geom_signif(comparisons = list(c("Absent","Present")), y_position = 95, tip_length = 0)+
  coord_cartesian(ylim=c(0, 100))
FigS16 <- plot_grid(proseq.subunits, groseq.subunits,chroseq.subunits, nrow = 1, align = "h", axis = 1, rel_widths = c(1,1,1.5), labels = "AUTO")
ggsave2("/local/workdir/bab364/PauseEvolution/FigS16_S2.final.pdf", plot = FigS16, device = "pdf", width = 6, height = 3, units = "in")

########################################################################################################################################
# calculate new pause index:

# this function will take (counts(newTSS to +150)/150) / (counts(+200 to TES)/(width-200))
calcPI <- function(annots, signal, output = c("mean", "individual")) {
  annots <- annots %>% filter(width(annots)>=400)
  signal <- signal %>% filter(signal$score>0)
  PC <- min(abs(signal$score))
  annots <- annots %>% filter(getCountsByRegions(signal, annots)>=5*PC)
  pause <- promoters(annots, upstream = 0, downstream = 150)
  genebody <- resize(annots, width = width(annots)-200, fix = "end")
  pause.count <- getCountsByRegions(signal, pause)+PC
  gb.counts <- getCountsByRegions(signal, genebody)+PC
  pause.norm <- pause.count/150
  gb.norm <- gb.counts/width(genebody)
  PI <- pause.norm/gb.norm
  if (output == "mean"){
    return(mean(PI))
  } else if (output == "individual"){
    return(PI)
  } else { print("output must be set to either mean or individual") }
}

########################################################################################################################################
# motif enrichment:

# find the pause:
findPause <- function(genes, proseq){
  search.window <- promoters(genes, upstream = 0, downstream = 100)
  proseq.counts <- getCountsByPositions(proseq, search.window, binsize = 1, field="score")
  print(dim(proseq.counts))
  index.max.proseq.counts <- numeric(length = nrow(proseq.counts))
  for (i in 1:length(search.window)) {
    max_index <- which.max(proseq.counts[i,])
    index.max.proseq.counts[i] <- max_index[1]
  }
  new.annots <- resizeGRanges(gr = search.window, resizeVector = index.max.proseq.counts, windowSize = 1) #output should be single base max pause
  pauses <- resize(new.annots, width = 1, fix = "center")
  pauses <- pauses %>% filter(start(pauses) > 50)
  return(pauses)
}
# calculate match scores for an individual sequence
calculate_match_score <- function(sequence, pwm, maxMatch) {
  n_bases <- ncol(pwm)
  n_positions <- nchar(sequence) - n_bases + 1
  match_scores <- numeric(n_positions)
  
  for (i in 1:n_positions) {
    subsequence <- substr(sequence, i, i + n_bases - 1)
    subsequence <- toupper(subsequence) #convert to upper case
    score <- 0
    for (j in 1:n_bases) {
      base <- substr(subsequence, j, j)
      score <- score + pwm[base, j]/maxMatch[j]
    }
    match_scores[i] <- score / n_bases
  }
  
  return(match_scores)
}
# load motif:
pause.button <- read.csv('/fs/cbsubscb17/storage/projects/AC_Projects/Lab_notebooks/Pause.evolution_project/PauseMotifAnalysis_Sep22/Gilad_use.defined.motif/Vivian_pause_motif_nospaces.csv', header = F)
rownames(pause.button)<- pause.button[,1]
pause.button<- pause.button[, colnames(pause.button) != "V1"]
colnames(pause.button)<- 1:ncol(pause.button)
# must be defined to run function
#maxMatch <- apply(pause.button, 2, max) #built in now
# load INR
initiator <- read.csv('/fs/cbsubscb17/storage/projects/AC_Projects/Lab_notebooks/Pause.evolution_project/PauseMotifAnalysis_Sep22/Gilad_use.defined.motif/Mammalian.Inr.csv', header = F, fill = T)[2:5,]
rownames(initiator)<- initiator[,1]
initiator<- initiator[, colnames(initiator) != "V1"]
colnames(initiator)<- 1:ncol(initiator)

# parallelize:
process_sequence <- function(k, sequences, pwm, maxMatch) {
  score.k <- calculate_match_score(as.character(sequences[[k]]), pwm, maxMatch)
  return(score.k)
}
# extract flanking sequences (take a 1kb region, located 500bp upstream or downstream of each pause): 
random.flank <- function(pause, seed = 42) {
  set.seed(seed)
  pause.1k <- resize(pause, width = 1000, fix = "center")
  upstream <- promoters(pause.1k, upstream = 1000, downstream = 0)
  downstream <- flank(pause.1k, width = 1000, start = F, both = F)
  # Create random indices
  random_indices <- runif(length(pause)) > 0.5
  # Initialize a GRanges object for the randomized regions
  randomized <- pause.1k
  # Assign regions based on random indices
  randomized[which(random_indices)] <- upstream[which(random_indices)]
  randomized[which(!random_indices)] <- downstream[which(!random_indices)]
  randomized[which(start(pause.1k) < 1001 & strand(pause.1k) == "+")] <- downstream[which(start(pause.1k) < 1001 & strand(pause.1k) == "+")]
  randomized[which(start(pause.1k) < 1001 & strand(pause.1k) == "-")] <- upstream[which(start(pause.1k) < 1001 & strand(pause.1k) == "-")]
  return(randomized)
}
# calculate fold-enrichment: 
motifEnrichment <- function(target, background, motif, cores = 20, method = "max") {
  maxMatch <- apply(motif, 2, max)
  scores <- mclapply(1:length(target), function(k) {
    process_sequence(k, target, motif, maxMatch)
  }, mc.cores = cores)
  scores.background <- mclapply(1:length(background), function(k) {
    process_sequence(k, background, motif, maxMatch)
  }, mc.cores = cores)
  if (method == "max") {
    enrichment <- mean(sapply(scores, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))), na.rm = T)/mean(sapply(scores.background, mean, na.rm = T), na.rm = T)
  } else {
    enrichment <- mean(sapply(scores, mean, na.rm = T), na.rm = T)/mean(sapply(scores.background, mean, na.rm = T), na.rm = T)
  }
  return(enrichment)
} #modified to take max by default, can disable and take "mean" if desired

#################################################################
# adjust background ranges based on sequence size, filters short contigs
extractSequenceFromFasta <- function(gr, fasta_file, rename.chrs = FALSE, BACKGROUND = FALSE) {
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file, use.names = FALSE)
  seq.with.names <- readDNAStringSet(fasta_file)
  
  # Rename chromosomes if set
  if (rename.chrs) {
    sequence.names <- numeric(length(sequences))
    for (j in 1:length(sequences)) {
      sequence.names[j] <- strsplit(names(seq.with.names[j]), " ")[[1]][1]
    } 
    names(sequences) <- sequence.names
  } else {
    names(sequences) <- names(seq.with.names)
  }
  
  sequences_at_site <- character()
  excluded_contigs <- c()
  
  for (i in seq_along(gr)) {
    seq_range <- gr[i]
    seq_name <- as.character(seqnames(gr)[i])
    seq_start <- start(seq_range)
    seq_end <- end(seq_range)
    seq_strand <- as.character(strand(seq_range))
    
    # Check if the sequence is longer than 2kb
    if (seqlengths(sequences)[seq_name] < 2000) {
      excluded_contigs <- c(excluded_contigs, seq_name)  # Store the name of the excluded contig
      next  # Skip to the next iteration if the sequence is shorter than 2kb
    }
    
    # Check if the coordinates are within the bounds of the sequence
    if (seq_start > 0 && seq_end <= seqlengths(sequences)[seq_name]) {
      # Extract the sequence using the valid coordinates
      if (seq_strand == "+") {
        sequences_at_site <- c(sequences_at_site, as.character(subseq(sequences[[seq_name]], seq_start, seq_end)))
      } else {
        sequences_at_site <- c(sequences_at_site, as.character(reverseComplement(subseq(sequences[[seq_name]], seq_start, seq_end))))
      }
    } else {
      # If the coordinates are out of bounds and BACKGROUND is enabled, shift the range upstream by 2kb
      if (BACKGROUND) {
        new_start <- max(1, seq_start - 2000)  # Ensure the start position is not negative
        new_end <- new_start + width(gr[i]) - 1  # Maintain the width of the original range
        sequences_at_site <- c(sequences_at_site, as.character(subseq(sequences[[seq_name]], new_start, new_end)))
        
        # Print the coordinates where BACKGROUND was applied
        cat("Range", i, "shifted upstream to:", new_start, "-", new_end, "\n")
      } else {
        next  # Skip to the next iteration if the range is out of bounds and BACKGROUND is not enabled
      }
    }
  }
  
  # Print the names of excluded contigs or chromosomes
  if (length(excluded_contigs) > 0) {
    cat("Excluded contigs/chromosomes due to being less than 2kb:", paste(excluded_contigs, collapse = ", "), "\n")
  }
  
  return(sequences_at_site)
}

#for multiple data input over same ranges
MetaProfilePlotMulti_sample_supplement <- function(..., ANNOTATIONS, plus.minus.ratio = 0.3, DATATYPE = "PRO", SPECIES, xmin=0, xmax=1000,aspect.ratio = 1/3, binsize = 5, norm=T){
  samples <- list(...)
  for (i in 1:length(samples)) {
    sample <- samples[[i]]
    assign("meansTMP", metaSubsample(sample, ANNOTATIONS, binsize = binsize,lower = 0.25, upper = 0.75, ncores = 4))
    if (i == 1){
      means <- meansTMP
      means[means == "plus"] <- paste0("plus",i)
      means$sample.color <- rep(sapply(substitute(list(...)), deparse)[-1][i], length(meansTMP))
      if(norm==T){
        print("signal normalization on")
        NF <- max(abs(means$mean))
        print(paste0("base normalization factor (sample 1 absolute maximum) = ",NF))
      }
    } else {
      meansTMP$sample.color <- rep(sapply(substitute(list(...)), deparse)[-1][i], length(meansTMP))
      meansTMP[meansTMP == "plus"] <- paste0("plus",i)
      if(norm==T){
        ND <- max(abs(meansTMP$mean))
        print(paste0("sample ",i , " absolute maximum = ",ND))
        print(paste0("sample ",i , " normalization factor = ",NF/ND))
        meansTMP$mean <- meansTMP$mean*(NF/ND)
        meansTMP$upper <- meansTMP$upper*(NF/ND)
        meansTMP$lower <- meansTMP$lower*(NF/ND)
      }
      means <- rbind(means, meansTMP)
    }
  }
  means$sample.color <- as.character(means$sample.color)
  means$sample.name <- as.character(means$sample.name)
  scalefactor <- round(max(abs(means$upper)), 6) # this establishes an upper bound for plotting. setting plus.minus.ratio sets the lower bound (-scalefactor*plus.minus.ratio)
  displaymax <- ifelse(scalefactor >= 10, round(max(means$upper), 1), round(max(means$upper), 2)) #this establishes a label to indicate scale
  PLOT <- ggplot(means, aes(x, mean, color = sample.color)) +
    geom_line() + theme_classic() + theme(aspect.ratio=aspect.ratio,
                                          axis.text.y = element_text(),legend.position = "none", axis.title.y = element_text(size=8), axis.title.x = element_text(size=10)) +
    geom_ribbon(aes(x, ymin = lower, ymax = upper, 
                    color = NULL, fill = sample.color),
                alpha = 0.2) + 
    labs(title = NULL,
         x = "Distance from TSS" ,
         y = "PRO-seq Signal + 50% CI") +
    scale_x_continuous(labels = function(x) x-(500))+
    scale_y_continuous(breaks = c(0,displaymax), sec.axis = sec_axis(trans= ~./(NF/ND), name = "PRO-cap signal + 50% CI", breaks = round(c(0,displaymax/(NF/ND)), 2)))+
    geom_vline(xintercept = 500, linetype = "dashed")+
    coord_cartesian(xlim = c(xmin,xmax), ylim=c(0, max(c(scalefactor,displaymax))))
  return(PLOT)
}

########################### Correlation between replicates ####################################

ec.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Ecoli.step.reann.rds")
ec.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ecoli/Ec-C1_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
ec.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ecoli/Ec-C2_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
ec.r3 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ecoli/Ec-C3_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(ec.r1, ec.genes),getCountsByRegions(ec.r2, ec.genes), method = "pearson")
cor.test(getCountsByRegions(ec.r2, ec.genes),getCountsByRegions(ec.r3, ec.genes), method = "pearson")
cor.test(getCountsByRegions(ec.r1, ec.genes),getCountsByRegions(ec.r3, ec.genes), method = "pearson")

hm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Hmediterranei.step.reann.rds")
hm.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Hmediterranei/HMB_5_CKDL210003158-1a-22_HFHGLCCX2_L4_dedup_QC_end.sort.bam",revcomp = T, trim.to = "3p", paired_end = T)
hm.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Hmediterranei/HMB_6_CKDL210003158-1a-AK1923_HFHGLCCX2_L4_dedup_QC_end.sort.bam",revcomp = T, trim.to = "3p", paired_end = T)
hm.r3 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Hmediterranei/HMB_7_CKDL210003158-1a-AK1932_HFHGLCCX2_L4_dedup_QC_end.sort.bam",revcomp = T, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(hm.r1, hm.genes),getCountsByRegions(hm.r2, hm.genes), method = "pearson")
cor.test(getCountsByRegions(hm.r2, hm.genes),getCountsByRegions(hm.r3, hm.genes), method = "pearson")
cor.test(getCountsByRegions(hm.r1, hm.genes),getCountsByRegions(hm.r3, hm.genes), method = "pearson")

dd.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Ddiscoideum.step.reann.rds")
dd.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ddiscoideum/dicty_r1_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
dd.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Ddiscoideum/dicty_r2_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(dd.r1, dd.genes),getCountsByRegions(dd.r2, dd.genes), method = "pearson")

sa.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Sarctica.step.reann.rds")
sa.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Sarctica/sphaeroforma_r1_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
sa.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Sarctica/sphaeroforma_r2_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(sa.r1, sa.genes),getCountsByRegions(sa.r2, sa.genes), method = "pearson")

cf.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Cfragrantissima.step.reann.rds")
cf.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Cfragrantissima/creolimax_r1_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cf.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Cfragrantissima/creolimax_r2_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(cf.r1, cf.genes),getCountsByRegions(cf.r2, cf.genes), method = "pearson")

co.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Cowczarzaki.step.reann.rds")
co.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Cowczarzaki/capsaspora_r1_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
co.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Cowczarzaki/capsaspora_r2_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(co.r1, co.genes),getCountsByRegions(co.r2, co.genes), method = "pearson")

nv.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Nvectensis.step.reann.rds")
nv.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Nvectensis/Cnidarian_JW1_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
nv.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Nvectensis/Cnidarian_JW2_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(nv.r1, nv.genes),getCountsByRegions(nv.r2, nv.genes), method = "pearson")

dp.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Dpulex.step.reann.rds")
dp.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Dpulex/Daphnia_Part1_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
dp.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Dpulex/Daphnia_Part2_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(dp.r1, dp.genes),getCountsByRegions(dp.r2, dp.genes), method = "pearson")

di.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Diulia.step.reann.rds")
di.r1 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Diulia/Dryas_iulia_F1_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
di.r2 <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Diulia/Dryas_iulia_H1_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(di.r1, di.genes),getCountsByRegions(di.r2, di.genes), method = "pearson")

sp.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Spurpuratus5.step.reann.rds")
sp.r1 <- import_bam("/fs/cbsubscb17/storage/data/short_read_index/sp5.0/mapped_PROseq/My_proseq_output_dir-05_21_2024_12hr/Sea_Urchin_Embryos_12_hours_dedup_QC_end.sort.bam",revcomp = T, trim.to = "3p", paired_end = T)
sp.r2 <- import_bam("/fs/cbsubscb17/storage/data/short_read_index/sp5.0/mapped_PROseq/My_proseq_output_dir-05_21_2024_20hr/Sea_Urchin_Embryos_20_hours_dedup_QC_end.sort.bam",revcomp = T, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(sp.r1, sp.genes),getCountsByRegions(sp.r2, sp.genes), method = "pearson")

pm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Pmarinus.step.reann.rds")
pm.f <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Pmarinus/Lamprey_Muscle_Female_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
pm.m <- import_bam("/local/workdir/bab364/PauseEvolution/PauseEvoGeoCorrected/Pmarinus/Lamprey_Muscle_Male_dedup_QC_end.sort.bam",revcomp = F, trim.to = "3p", paired_end = T)
cor.test(getCountsByRegions(pm.f, pm.genes),getCountsByRegions(pm.m, pm.genes), method = "pearson")

