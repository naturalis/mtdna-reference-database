#load in the needed pacakges
library(pegas)
library(adegenet)

#read in the data, x is the fasta created by the two_seq_matrix.py script and 
#r_msa is created by the for_r_msa.py script. For reading in the r_msa the 
#right columns have to be specified, the col.pop is the most important and has
#to be the column containing the country.
#r_msa is only needed for the country column, so any file containing 
r_msa <- read.loci("../results/nuclear/msa_P78406_S4_for_r.csv", loci.sep = "\t"
                   , col.loci = 2:8, col.pop = 11, row.names = 1)
fasta_reads <- read.dna("../results/nuclear/reads_P78406.fas", "f")

#check whether its all fine and create haplotypes.
alview(fasta_reads)
h <- haplotype(fasta_reads)

#calculate distances and rmst
d <- dist.dna(h, "N")
nt <- rmst(d, quiet = TRUE)

#plot without the colour country codes
(sz <- summary(h))
(nt.labs <- attr(nt, "labels"))
sz <- sz[nt.labs]
plot(nt, size = sz)

#plot with the colour country codes
pop <- r_msa$population
(P <- haploFreq(fasta_reads, fac = pop, haplo = h))
P <- P[nt.labs, ]
plot(nt, size = sz, pie = P, legend = TRUE)
xy <- replot()