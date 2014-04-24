\name{DMRcate-package}
\alias{DMRcate-package}
\alias{DMRcate}
\docType{package}
\title{
Illumina 450K methylation array spatial analysis
}
\description{
  \emph{De novo} identification and extraction of differentially
  methylated regions (DMR) in the human genome using 450k array
  data. \code{DMRcate} extracts and annotates differentially methylated regions
  (DMRs) using an array-bias corrected density estimate. Functions are
  provided for filtering probes possibly confounded by SNPs and
  cross-hybridisation. Includes bedGraph and plotting functions.
}
\author{
Tim J. Peters <Tim.Peters@csiro.au>
}
\references{
Peters T.J., Buckley M.J., Pidsley R., Clark S.J., Molloy P.L. \emph{De
Novo} identification of differentially methylated regions in the human
genome. To submit
}

\examples{
data(dmrcatedata)
myMs <- logit2(myBetas)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- annotate(myMs.noSNPs, analysis.type="differential",
    design=design, coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01)
dmrcoutput <- dmrcate(myannotation, bw=1000)
makeBedgraphs(dmrcoutput=dmrcoutput, betas=myBetas, samps=c(1,39))
phen.col <- c(rep("orange", 38), rep("blue", 38))
DMR.plot(dmrcoutput=dmrcoutput, dmr=1, betas=myBetas, phen.col=phen.col,
    pch=16, toscale=TRUE)
}