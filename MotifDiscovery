#R --vanilla --slave --args $(pwd) HT_allReads_TSS_maxTSNs+-10.txt HT_allReads_TSS_maxTSNs+-10_SeqLogo.pdf < getSeqLogo.R
# only use bed regions with at least 5 VALID mat AND 5 pat VALID reads

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input=args[2]
output=args[3]



#s<- read.table("/Volumes/SPC_SD/KD_IGV/TSS/identifyTSS/HT_allReads_TSS_maxTSNs+-10.txt", sep="")
s<-read.table(input)
seq<- data.frame(do.call(rbind, strsplit(as.character(s$V1), "")))
head(s)
head(seq)
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

library(seqLogo)

#create data frame using the four vectors
df <- data.frame(a,c,g,t)
df

#define function that divides the frequency by the row sum i.e. proportions
proportion <- function(x){
  rs <- sum(x);
  return(x / rs);
}

#create position weight matrix
pwm <- apply(df, 1, proportion)
p <- makePWM(pwm)
# slotNames(p)
# p@consensus
# p@ic
# p@width
# p@alphabet
pdf(output)
seqLogo(p)
dev.off()

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
    b <- c(a, sum(seq[,i]=="A"))
    t <- c(t, sum(seq[,i]=="T"))
    c <- c(c, sum(seq[,i]=="C"))
    g <- c(g, sum(seq[,i]=="G"))
  }
  return (data.frame(a,c,g,t))
}


SeqLogo <- function(seq, output) {
  #seq=m$HighAlleleSeq
  seq<- data.frame(do.call(rbind, strsplit(as.character(seq), "")))
  df <- seq_upperCase(seq)
  #create position weight matrix
  pwm <- apply(df, 1, proportion)
  p = makePWM((pwm[,190:210]))
  #p <- makePWM(pwm)
  # slotNames(p)
  # p@consensus
  # p@ic
  # p@width
  # p@alphabet
  pdf(output)
  seqLogo(p)
  dev.off()
  return (pwm)
}
