# motif analysis from Shao-Pei via Alex:


library(seqLogo)
#define function that divides the frequency by the row sum i.e. proportions
proportion <- function(x){
  rs <- sum(x);
  return(x / rs);
}

seq_upperCase <- function(seq){
  seq[seq=="a"]<- "A"
  seq[seq=="t"]<- "T"
  seq[seq=="c"]<- "C"
  seq[seq=="g"]<- "G"
  
  a <- NULL
  t <- NULL
  c <- NULL
  g <- NULL
  for (i in 1:NCOL(seq)){
    a <- c(a, sum(seq[,i]=="A"))
    t <- c(t, sum(seq[,i]=="T"))
    c <- c(c, sum(seq[,i]=="C"))
    g <- c(g, sum(seq[,i]=="G"))
  }
  return (data.frame(a,c,g,t))
}

SeqLogo <- function(seq, output, range=NULL) {
  #seq=m$HighAlleleSeq
  seq<- data.frame(do.call(rbind, strsplit(as.character(seq), "")))
  df <- seq_upperCase(seq)
  #create position weight matrix
  pwm <- apply(df, 1, proportion)
  if (!is.null(range)){
    p = makePWM((pwm[,range]))    
  }else{
    p = makePWM((pwm))  
  }
  print(p)
  #p <- makePWM(pwm)
  # slotNames(p)
  # p@consensus
  # p@ic
  # p@width
  # p@alphabet
  pdf(output, width = 3, height = 4.5)
  seqLogo(p)
  dev.off()
  print(pwm) #originally return(pwm), leaves seqlogo to be output to pdf
  return(seqLogo(p))
} # note: this function has no special instructions for handling vectors of different lengths. It will loop over shorter sequences.

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
  endPositions <- startPositions # +1 for even-number-sized ranges
  # Construct and return the resized GRanges
  output <- GRanges(seqnames = seqnames(gr), ranges = IRanges(start = startPositions, end = endPositions), strand = strand(gr))
  return(resize(output, width = 1, fix = "start"))
} #note: this function has been modified for motif analysis in windows of odd-number size! Do not use for other applications! 
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
  pauses <- resize(new.annots, width = 17, fix = "center")
  return(pauses)
}
extractSequenceFromFasta <- function(gr, fasta_file, rename.chrs = F) {
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file, use.names = F)
  seq.with.names <- readDNAStringSet(fasta_file)
  # rename weird chromosomes if set 
  if (rename.chrs) {
    sequence.names <- numeric(length(sequences))
    for (j in 1:length(sequences)) {
      sequence.names[j] <- strsplit(names(seq.with.names[j])," ")[[1]][1]
    } 
    names(sequences)<- sequence.names
  } else {
    names(sequences) <- names(seq.with.names)
  }
  sequences_at_site <- character(length(gr))
  for (i in seq_along(gr)) {
    seq_range <- gr[i]
    seq_name <- as.character(seqnames(gr)[i])
    seq_start <- start(seq_range)
    seq_end <- end(seq_range)
    seq_strand <- as.character(strand(seq_range))
    sequences_at_site[i] <- ifelse(seq_strand == "+", as.character(subseq(sequences[[seq_name]], seq_start, seq_end)), as.character(reverseComplement(subseq(sequences[[seq_name]], seq_start, seq_end))))
  }
  return(sequences_at_site)
}


############################################################################################
# test motif analysis with pause site extraction, generate plot for human
HUM_3p <- import_bam("/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/PE_data/H.sapiens/ChROseq_merged_0h.bam",revcomp = F, trim.to = "3p", paired_end = T)
HUM_3p <- tidyChromosomes(HUM_3p)
hs.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Hsapiens.step.reann.rds")
hs.genes <- tidyChromosomes(hs.genes)
hs.genes <- hs.genes %>% filter(getCountsByRegions(HUM_3p,promoters(hs.genes, upstream = 0,downstream = 100))!=0)
hs.pause <- findPause(hs.genes, HUM_3p)
hs.sequences <- extractSequenceFromFasta(hs.pause, "/fs/cbsubscb17/storage/data/short_read_index/hg19/hg19.fa.gz")
#SeqLogo(hs.sequences, "/local/workdir/bab364/tmp/human_logo.pdf")
SeqLogo(hs.sequences, "/local/workdir/bab364/PauseEvolution/human_logo.pdf")

# pipeline works, repeat for lamprey using direct extraction of bases from fasta
LAM_PE_m <- import_bigWig("/local/workdir/James/PauseEvolution/data/Lampetra_fish/Lamprey_Muscle_Male_dedup_QC_end_plus.rpm.bw", "/local/workdir/James/PauseEvolution/data/Lampetra_fish/Lamprey_Muscle_Male_dedup_QC_end_minus.rpm.bw")
LAM_PE_f <- import_bigWig("/fs/cbsubscb17/storage/projects/NASA_2020/Lamprey_Muscle_Female_3Map/Lamprey_Muscle_Female_dedup_QC_end_plus.rpm.bw","/fs/cbsubscb17/storage/projects/NASA_2020/Lamprey_Muscle_Female_3Map/Lamprey_Muscle_Female_dedup_QC_end_minus.rpm.bw")
LAM_PE <- mergeGRangesData(LAM_PE_m,LAM_PE_f)
rm(LAM_PE_m,LAM_PE_f)
pm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Pmarinus.step.reann.rds")
pm.genes <- pm.genes %>% filter(getCountsByRegions(LAM_PE,promoters(pm.genes, upstream = 0,downstream = 100))!=0)
pm.pause <- findPause(pm.genes, LAM_PE)
pm.sequences <- extractSequenceFromFasta(pm.pause, "/local/storage/data/short_read_index/petMar2/petMar2.fa")
SeqLogo(pm.sequences, "/local/workdir/bab364/PauseEvolution/lamprey_logo.pdf")

# function to extract sequences directly from fasta works, now proceed to proto-paused organisms
# S.pombe
SPOM_PE <- import_bigWig("/local/workdir/James/PauseEvolution/data/S.pombe/Pombe_WT_merged_plus.rpm.bw","/local/workdir/James/PauseEvolution/data/S.pombe/Pombe_WT_merged_minus.rpm.bw")
spom.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Spombe.step.reann.rds")
length(spom.genes)
spom.genes <- spom.genes %>% filter(getCountsByRegions(SPOM_PE,promoters(spom.genes, upstream = 0,downstream = 100))!=0)
length(spom.genes)
spom.pause <- findPause(spom.genes, SPOM_PE)
spom.sequences <- extractSequenceFromFasta(spom.pause, "/local/storage/data/short_read_index/S_pombe/Schizosaccharomyces_pombe.ASM294v2.simple.chrs.fa")
SeqLogo(spom.sequences, "/local/workdir/bab364/PauseEvolution/pombe_logo.pdf")


# O.sativa
OSAT_PE <- import_bigWig("/fs/cbsubscb17/storage/data/short_read_index/O.sativa/O.sativa_3p.end_plus.bw","/fs/cbsubscb17/storage/data/short_read_index/O.sativa/O.sativa_3p.end_minus.bw")
os.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Osativa.step.reann.rds")
length(os.genes)
os.genes <- os.genes %>% filter(getCountsByRegions(OSAT_PE,promoters(os.genes, upstream = 0,downstream = 100))!=0)
length(os.genes)
os.pause <- findPause(os.genes, OSAT_PE)
os.sequences <- extractSequenceFromFasta(os.pause, "/local/storage/data/short_read_index/O.sativa/O.sativa_GCF_001433935.1_IRGSP-1.0_genomic.fna", rename.chrs = T)
SeqLogo(os.sequences, "/local/workdir/bab364/PauseEvolution/sativa_logo.pdf")

# Z. mays
ZMAYS_PE2 <- import_bigWig("/fs/cbsubscb17/storage/data/short_read_index/Zea.mays/Z.mays_3p.end_plus.bw","/fs/cbsubscb17/storage/data/short_read_index/Zea.mays/Z.mays_3p.end_minus.bw")
zm.genes <- readRDS(file = "/fs/cbsubscb17/storage/data/Brent_TSS.annotations_test_10May23/Trim_PROseq_to.5prime/Reann.TSS.stepfunc.BB/Zmays.step.reann.rds")
length(zm.genes)
zm.genes <- zm.genes %>% filter(getCountsByRegions(ZMAYS_PE2,promoters(zm.genes, upstream = 0,downstream = 100))!=0)
length(zm.genes)
zm.pause <- findPause(zm.genes, ZMAYS_PE2)
zm.sequences <- extractSequenceFromFasta(zm.pause, "/local/storage/data/short_read_index/Zea.mays/ncbi-genomes-2022-09-20/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna", rename.chrs = T)
zm.logo <- SeqLogo(zm.sequences, "/local/workdir/bab364/PauseEvolution/mays_logo.pdf")


############################################################################################

