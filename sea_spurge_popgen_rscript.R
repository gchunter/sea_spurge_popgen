#Load dartR
library(dartR)
#Assign dartseq file and ind file to a genlight object
gl_all_samples <- gl.read.dart(filename = "data/Report_DEup19-4268_1_moreOrders_SNP_1.csv",
                   ind.metafile = "data/sea_spurge_ind2_metrics.csv")
#structure of the genlight object
gl_all_samples
#number of individuals in the genlight object
nInd(gl_all_samples)
#number of loci
nLoc(gl_all_samples)
#number of populations
nPop(gl_all_samples)
#populaitons to which the individuals belong
pop(gl_all_samples)
#table of individuals per population
table(pop(gl_all_samples))
#population names
levels(pop(gl_all_samples))
#calculate and report the overall call rate across the loci
#Use a missing threshold value < 10% (2202 loci retained)
gl.report.callrate(gl_all_samples)
#filter on call rate and retain only those loci with less than 10% missing data
gl_filter1_callrate <- gl.filter.callrate(gl_all_samples, method = "loc",
                                          threshold = 0.90, recalc = TRUE, plot = TRUE)
#filter on departure of Hardy-Weinberg equilibrium for every loci
gl_filter2_hwe <- gl.filter.hwe(gl_filter1_callrate,
                                alpha = 0.05, basis = "any",
                                bon = TRUE, v = 5)
#report the minor allele frequencies for filter2_HWE
gl.report.maf(gl.filter2_hwe, v = 5)
#filter on MAF at 0.01
gl_filter3_maf <- gl.filter.maf(gl_filter2_hwe, threshold = 0.01, v = 5)
#report the reproducibility for the loci
gl.report.repavg(gl_filter3_maf)
#filter out loci for thich the repeateability is less than 0.98
gl_filter4_repavg <- gl.filter.repavg(gl_filter3_maf, threshold = 0.98, v = 3)
#report the number of monomorphs
gl.report.monomorphs(gl_filter4_repavg)
#report secondaries
gl.report.secondaries(gl_filter4_repavg)
#filter out secondaries
gl_filter5_secondaries <- gl.filter.secondaries(gl_filter4_repavg, method = "random", v = 5)
#produce a smearplot
glPlot(gl_filter5_secondaries)



