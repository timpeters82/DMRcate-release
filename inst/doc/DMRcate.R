
## ----bioconductor, message=FALSE, warning=FALSE, eval=FALSE--------------
## source("http://bioconductor.org/biocLite.R")
## biocLite("DMRcate")


## ----libr, message=FALSE, warning=FALSE----------------------------------
library(DMRcate)


## ----loaddata------------------------------------------------------------
data(dmrcatedata)
myMs <- logit2(myBetas)


## ----filter--------------------------------------------------------------
nrow(illuminaSNPs)
nrow(myMs)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
nrow(myMs.noSNPs)


## ----annotate------------------------------------------------------------
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- annotate(myMs.noSNPs, analysis.type="differential", design=design, coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01)


## ----length--------------------------------------------------------------
length(myannotation$ID)


## ----dmrcate-------------------------------------------------------------
dmrcoutput <- dmrcate(myannotation, bw=1000)


## ----plotting------------------------------------------------------------
dmrcoutput$results[1,]
DMR.plot(dmrcoutput=dmrcoutput, dmr=1, betas=myBetas, phen.col=c(rep("orange", 38), rep("blue", 38)), pch=16, toscale=TRUE)


## ----sessionInfo---------------------------------------------------------
sessionInfo()


