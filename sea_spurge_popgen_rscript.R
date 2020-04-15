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
#report call rate for individuals
gl.report.callrate(gl_filter5_secondaries, method = "ind", plot = TRUE, v = 2)
#filter call rate for individuals at 95%
gl_filter6_ind_cr <- gl.filter.callrate(gl_filter5_secondaries,
                                        method = "ind", threshold = 0.95,
                                        recalc = TRUE, plot = TRUE)
gl_filter7_loc_cr <- gl.filter.callrate(gl_filter6_ind_cr,
                                        method = "loc", threshold = 0.95, recalc = TRUE, 
                                        mono.rm = TRUE, plot = TRUE, v = 2)

#population names in gl_filter7_loc_cr
popNames(gl_filter7_loc_cr)
#convert genlight to treemix input file format
gl2treemix(gl_filter7_loc_cr, outpath = "data")
#change from genlight to genind
gind_filter7_loc_cr <- gl2gi(gl_filter7_loc_cr, v =1)
#Load Hierfstat
library(hierfstat)
#convert genind to hierfstat format
hf_filter7 <- genind2hierfstat(gind_filter7_loc_cr, pop = NULL)
attach(hf_filter7)
hf_filter7_loci <- hf_filter7[,c(2:752)]
hf_filter7_pop <- hf_filter7$pop

hf_filter7_basic_stats <- basic.stats(data.frame(hf_filter7_pop, 
                                                 hf_filter7_loci),
            diploid = TRUE, digits = 4)
#calculate allelic richness
allelic.richness(hf_filter7)
#calculate basic statistics
hf_filter7_basic_stats <- basic.stats(hf_filter7)

hf_filter7

?genind2hierfstat

write.csv(hf_filter7, file = "hf_filter7.csv")





head(hf_filter7)
table(diploid[,1])
table(hf_filter7[,1])
data(gind_filter7_loc_cr)
head(genind2hierfstat(gind_filter7_loc_cr)[,1:10])
dim(hf_filter7)
gind_filter7_loc_cr
allele.count(hf_filter7, diploid = TRUE)
allelic.richness(hf_filter7, diploid = TRUE)
hf_ar <- allelic.richness(hf_filter7, diploid = TRUE)
rm(hf_ar)
hf_boot_fis <- boot.ppfis(hf_filter7,
                          nboot = 10,
                          quant = c(0.025, 0.975),
                          diploid = TRUE,
                          dig = 3)
hf_boot_fis
ho <- hf_filter7_basic_stats[["Ho"]]
library(poppr)
genind2genalex(gind_filter7_loc_cr,
               filename = "gind_filter7_loc_cr.csv",
               overwrite = FALSE,
               quiet = FALSE,
               pop = NULL,
               allstrata = TRUE,
               geo = FALSE,
               sep = ",",
               sequence = TRUE)
gind_filter7_loc_cr
genind2genalex(gind_filter7_loc_cr,
               filename = "gind2.csv")
library(adegenet)
