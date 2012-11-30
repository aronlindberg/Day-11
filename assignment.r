# Load libraries
library(TraMineR)
library(TraMineRextras)

# Load and format data
data(biofam)
biofam$cohort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960),
    labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"),
    right=FALSE)
bf.states <- c("Parent", "Left", "Married", "Left/Married", "Child",
    "Left/Child", "Left/Married/Child", "Divorced")
bf.shortlab <- c("P","L","M","LM","C","LC", "LMC", "D")
biofam.seq <- seqdef(biofam[,10:25], states=bf.shortlab,
    labels=bf.states, weights=biofam$wp00tbgs)
weight <- attr(biofam.seq, "weight")

#1. The biofam data was created from time stamped event data. Table 1 describes how the states were defined in terms of event occurrences.
## First we define the assumed events behind each state-transition.
tm <- c(
   ## "P",  "L",   "M",   "LM",  "C",     "LC",    "LMC",   "D"
   ##------------------------------------------------------------
    "P",   "L",   "M",   "L,M", "C",     "L,C",   "L,M,C", "M,D",
    "P",   "L",   "P,M", "M",   "P,C",   "C",     "M,C",   "M,D",
    "",    "D,L", "M",   "L",   "D,C",   "D,L,C", "L,C",   "D",
    "P,D", "D",   "P",   "L,M", "D,P,C", "D,C",   "C",     "D",
    "P",   "L",   "M",   "L,M", "C",     "L",     "L,M",   "M,D",
    "P",   "",    "M",   "M",   "P",     "L,C",   "M",     "D",
    "P,D", "L",   "",    "",    "D,P",   "D",     "L,M,C", "D",
    "P",   "L",   "P,M", "M",   "P,C",   "C",     "M,C",   "D"
    )
tm <- matrix(tm, nrow=8, ncol=8, byrow=TRUE)

## Assigning row and column names from an automatic transition matrix
dimnames(tm) <- dimnames(seqetm(biofam.seq, method="transition",
use.labels=FALSE))
tm
## Transforming the state sequences into event sequences
bf.seqe <- seqecreate(biofam.seq, tevent=tm, use.labels=FALSE)

#2. Plot the event sequences by sex. What is the main difference between men and women.
seqpcplot(bf.seqe, group=biofam$sex, alpabet=c("P", "L", "M", "C", "D"),
	filter = list(type = "function",
		value = "minfreq",
		level = 0.15),
	ltype="non-embeddable",
	cex=1.5, lwd=1
)


#3. Plot the event sequences by the birth cohorts defined in Assignment 10 (Before end of word war II versus after end of WW-II). Represent non-embeddable sequence patterns and color only those with a support of at least 15%. Comment differences between the two cohorts.
biofam$cohort <- cut(biofam$birthyr, c(1900,1946,1960),
	labels-c("Before WWII", "After WWII"), right=FALSE)

seqpcplot(bf.seqe, group=biofam$cohort, alphabet=c("P", "L", "M", "C", "D"),
	filter = list(type = "function",
		value = "minfreq",
		level = 0.15),
	ltype="non-embeddable",
	cex=1.5, lwd=1
)

#4. Find the most frequent subsequences (minimum support of 10%) with at least 2 events. Among those who left home and got married, what is the proportion who did it the same year? Plot the 10 most frequent subsequences.
bf.sub <- seqefsub(bf.seqe, pMinSupport=.1)

## computing the number of events by subsequence
bf.fsubsn <- seqentrans(bf.fsub)

## selecting subsequences with 2 or more events
bf.fsb <- bf.fsubn[bf.fsubn$data$nevent >1]
bf.fsb
plot(bf.fsb[1:10,])

#5. Display and plot the 10 subsequences which best discriminate (a) women from men, and (b) birth cohorts.

## between sex
bf.discr.sex <- seqecmprgroup(bf.fsub, group=biofam$sex)
bf.discr.sex[1:10,])

## between birth cohorts
bf.discr.bthc <- seqecmpgroup(bf.fsub, group=biofam$cohort)
bf.discr.bthc[1:10,]
plot(bf.discr.bthc[1:10,])
