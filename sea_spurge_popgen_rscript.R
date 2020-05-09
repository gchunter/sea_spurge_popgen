#Load packages
library(dartR)
library(hierfstat)
library(adegenet)
library(tidyverse)
library(ggplot2)
library(strataG)
library(poppr)
library(ggplot2)
library(pegas)
library(StAMPP)

#1. SAMPLE SNP DATA

#Assign dartseq file and ind file to a genlight object
gl_all_samples <- gl.read.dart(filename = "data/Report_DEup19-4268_1_moreOrders_SNP_1.csv",
                   ind.metafile = "data/sea_spurge_ind4_metrics.csv")

#2. SNP FILTERING

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
ss_gl_f1 <- gl.filter.callrate(gl_all_samples, method = "loc",
                                          threshold = 0.90, recalc = TRUE, plot = TRUE)

#filter on departure of Hardy-Weinberg equilibrium for every loci
ss_gl_f2 <- gl.filter.hwe(ss_gl_f1,
                                alpha = 0.05, basis = "any",
                                bon = TRUE, v = 5)

#report the minor allele frequencies for ss_gl_f2
gl.report.maf(ss_gl_f2, v = 5)

#filter on MAF at 0.01
ss_gl_f3 <- gl.filter.maf(ss_gl_f2, threshold = 0.01, v = 5)

#report the reproducibility for the loci
gl.report.repavg(ss_gl_f3)

#filter out loci for thich the repeateability is less than 0.98
ss_gl_f4 <- gl.filter.repavg(ss_gl_f3, threshold = 0.98, v = 3)

#report the number of monomorphs
gl.report.monomorphs(ss_gl_f4)

#report secondaries
gl.report.secondaries(ss_gl_f4)

#filter out secondaries
ss_gl_f5 <- gl.filter.secondaries(ss_gl_f4, method = "random", v = 5)

#produce a smearplot
glPlot(ss_gl_f5)

#report call rate for individuals
gl.report.callrate(ss_gl_f5, method = "ind", plot = TRUE, v = 2)

#filter call rate for individuals at 95%
ss_gl_f6 <- gl.filter.callrate(ss_gl_f5,method = "ind", threshold = 0.95,
                                        recalc = TRUE, plot = TRUE)

#filte on call rate for loci at 95%
ss_gl_f7 <- gl.filter.callrate(ss_gl_f6, method = "loc", threshold = 0.95, recalc = TRUE, 
                                        mono.rm = TRUE, plot = TRUE, v = 2)

#Recode the population names to reflect the states they were collected in
ss_gl_f7_state <- gl.edit.recode.pop(ss_gl_f7)

#Recode the population names to reflect the bioregion where they were collected
ss_gl_f7_bioreg <- gl.edit.recode.pop(ss_gl_f7)

#table of number of samples at the state level
table(pop(ss_gl_f7_state))

#population names in ss_gl_f7
popNames(ss_gl_f7)

#number of individuals in each population in ss_gl_f7
table(pop(ss_gl_f7))

#3. GENETIC DIVERSITY

#change from genlight to genind

#population level
ss_gl_f7_gi <- gl2gi(ss_gl_f7, v =1)

#state level
ss_gl_f7_gi_state <- gl2gi(ss_gl_f7_state, v = 1)

#bioreg level
ss_gl_f7_gi_bioreg <- gl2gi(ss_gl_f7_bioreg, v = 1)

#convert genind(gi) to genpop(gp)
#population level
ss_gl_f7_gi_gp <- genind2genpop(ss_gl_f7_gi, pop = NULL,
                                quiet = FALSE, process.other = TRUE)

#state level
ss_gl_f7_gi_gp_state <- genind2genpop(ss_gl_f7_gi_state, pop = NULL,
                                quiet = FALSE, process.other = TRUE)

#bioreg level
ss_gl_f7_gi_gp_bioreg <- genind2genpop(ss_gl_f7_gi_bioreg, pop = NULL,
                                      quiet = FALSE, process.other = TRUE)


#convert from genind (gi) to gtypes (gt) format for StrataG package
#population level
ss_gl_f7_gt <- genind2gtypes(ss_gl_f7_gi)

#state level
ss_gl_f7_gt_state <- genind2gtypes(ss_gl_f7_gi_state)

#bioreg level
ss_gl_f7_gt_bioreg <- genind2gtypes(ss_gl_f7_gi_bioreg)

#Calculate summary stats for each stratum (population) seperately
#population level
ss_gl_f7_gt_allsum <- summarizeLoci(ss_gl_f7_gt, by.strata = TRUE)

#state level
ss_gl_f7_gt_state_allsum <- summarizeLoci(ss_gl_f7_gt_state, by.strata = TRUE)

#bioreg level
ss_gl_f7_gt_bioreg_allsum <- summarizeLoci(ss_gl_f7_gt_bioreg, by.strata = TRUE)

#calculate the number of private alleles in each strata and locus with StrataG
#population level
ss_gl_f7_gt_priv_allele <- privateAlleles(ss_gl_f7_gt)

pop_pri_alle_df <- data.frame(ss_gl_f7_gt_priv_allele) #change to data.frame
colSums(pop_pri_alle_df)

#state level
ss_gl_f7_gt_state_priv_allele <- privateAlleles(ss_gl_f7_gt_state)

state_pri_alle_df <- data.frame(ss_gl_f7_gt_state_priv_allele) #change to data.frame
colSums(state_pri_alle_df) #sum each column

#bioreg level
ss_gl_f7_gt_bioreg_priv_allele <- privateAlleles(ss_gl_f7_gt_bioreg)

bioreg_pri_alle_df <- data.frame(ss_gl_f7_gt_bioreg_priv_allele) #change to data.frame
colSums(bioreg_pri_alle_df) #sum each column

#convert genind to hierfstat
ss_gl_f7_hf <- genind2hierfstat(ss_gl_f7_gi) #population level
ss_gl_f7_hf_state <- genind2hierfstat(ss_gl_f7_gi_state) #state level
ss_gl_f7_hf_bioreg <- genind2hierfstat(ss_gl_f7_gi_bioreg) #bioregion level

ss_gl_f7_hf_loci <- ss_gl_f7_hf[,c(2:752)] #specify the loci
ss_gl_f7_hf_pop <- ss_gl_f7_hf$pop #specify a vector of populaitons
attach(ss_gl_f7_hf) #attach the hierfstat object
ss_gl_f7_hf_basic_stats <- basic.stats(data.frame(ss_gl_f7_hf_pop,
                                                  ss_gl_f7_hf_loci),
                                       diploid = TRUE, digits = 4)

ss_gl_f7_hf_state_basic_stats <- basic.stats(ss_gl_f7_hf_state,
                                             diploid = TRUE, digits = 4)

ss_gl_f7_hf_bioreg_basic_stats <- basic.stats(ss_gl_f7_hf_bioreg,
                                              diploid = TRUE, digits = 4)

ss_gl_f7_hf_basic_stats_Ho <- ss_gl_f7_hf_basic_stats[["Ho"]] #isolate Ho
hf_basic_stats_Ho_df <- data.frame(ss_gl_f7_hf_basic_stats_Ho) # Ho to data frame
hf_basic_stats_Ho_df[is.na(hf_basic_stats_Ho_df)] = 0 #change NA to 0
colMeans(hf_basic_stats_Ho_df) #means of all columns in Ho

ss_gl_f7_hf_basic_stats_Hs <- ss_gl_f7_hf_basic_stats[["Hs"]] #isolate Hs
hf_basic_stats_Hs_df <- data.frame(ss_gl_f7_hf_basic_stats_Hs) # Hs to data frame
hf_basic_stats_Hs_df[is.na(hf_basic_stats_Hs_df)] = 0 #change NA's to 0
colMeans(hf_basic_stats_Hs_df) #means of all columns in Hs

ss_gl_f7_hf_basic_stats_Fis <- ss_gl_f7_hf_basic_stats[["Fis"]] #isoalte Fis
ss_gl_f7_hf_basic_stats_Fis_df <- data.frame(ss_gl_f7_hf_basic_stats_Fis)#Fis to data frame
ss_gl_f7_hf_basic_stats_Fis_df[is.na(ss_gl_f7_hf_basic_stats_Fis_df)] = 0 #change NA to 0
colMeans(ss_gl_f7_hf_basic_stats_Fis_df) #means of all columns in Fis

ss_gl_f7_hf_state_basic_stats_Ho <- ss_gl_f7_hf_state_basic_stats[["Ho"]]
hf_state_basic_stats_Ho_df <- data.frame(ss_gl_f7_hf_state_basic_stats_Ho)
hf_state_basic_stats_Ho_df[is.na(hf_state_basic_stats_Ho_df)] = 0
colMeans(hf_state_basic_stats_Ho_df)

ss_gl_f7_hf_state_basic_stats_Hs <- ss_gl_f7_hf_state_basic_stats[["Hs"]]
hf_state_basic_stats_Hs_df <- data.frame(ss_gl_f7_hf_state_basic_stats_Hs)
hf_state_basic_stats_Hs_df[is.na(hf_state_basic_stats_Hs_df)] = 0
colMeans(hf_state_basic_stats_Hs_df)

ss_gl_f7_hf_state_basic_stats_Fis <- ss_gl_f7_hf_state_basic_stats[["Fis"]]
hf_state_basic_stats_Fis_df <- data.frame(ss_gl_f7_hf_state_basic_stats_Fis)
hf_state_basic_stats_Fis_df[is.na(hf_state_basic_stats_Fis_df)] = 0
colMeans(hf_state_basic_stats_Fis_df)

ss_gl_f7_hf_bioreg_basic_stats_Ho <- ss_gl_f7_hf_bioreg_basic_stats[["Ho"]]
ss_gl_f7_hf_bioreg_basic_stats_Ho_df <- data.frame(ss_gl_f7_hf_bioreg_basic_stats_Ho)
ss_gl_f7_hf_bioreg_basic_stats_Ho_df[is.na(ss_gl_f7_hf_bioreg_basic_stats_Ho_df)] = 0
colMeans(ss_gl_f7_hf_bioreg_basic_stats_Ho_df)

ss_gl_f7_hf_bioreg_basic_stats_Hs <- ss_gl_f7_hf_bioreg_basic_stats[["Hs"]]
ss_gl_f7_hf_bioreg_basic_stats_Hs_df <- data.frame(ss_gl_f7_hf_bioreg_basic_stats_Hs)
ss_gl_f7_hf_bioreg_basic_stats_Hs_df[is.na(ss_gl_f7_hf_bioreg_basic_stats_Hs_df)] = 0
colMeans(ss_gl_f7_hf_bioreg_basic_stats_Hs_df)

ss_gl_f7_hf_bioreg_basic_stats_Fis <- ss_gl_f7_hf_bioreg_basic_stats[["Fis"]]
ss_gl_f7_hf_bioreg_basic_stats_Fis_df <- data.frame(ss_gl_f7_hf_bioreg_basic_stats_Fis)
ss_gl_f7_hf_bioreg_basic_stats_Fis_df[is.na(ss_gl_f7_hf_bioreg_basic_stats_Fis_df)] = 0
colMeans(ss_gl_f7_hf_bioreg_basic_stats_Fis_df)

#Allelic Richness
#population level
f7_hf_ar <- allelic.richness(ss_gl_f7_hf, min.n = NULL, diploid = TRUE)
f7_hf_ar_df <- data.frame(f7_hf_ar)
f7_hf_ar_df[is.na(f7_hf_ar_df)] = 0
colMeans(f7_hf_ar_df)

#state level
f7_hf_state_ar <- allelic.richness(ss_gl_f7_hf_state, min.n = NULL, diploid = TRUE)
f7_hf_state_ar_df <- data.frame(f7_hf_state_ar)
f7_hf_state_ar_df[is.na(f7_hf_state_ar_df)] = 0
colMeans(f7_hf_state_ar_df)

#bioregion level
f7_hf_bioreg_ar <- allelic.richness(ss_gl_f7_hf_bioreg, min.n = NULL, diploid = TRUE)
f7_hf_bioreg_ar_df <- data.frame(f7_hf_bioreg_ar)
f7_hf_bioreg_ar_df[is.na(f7_hf_bioreg_ar_df)] = 0
colMeans(f7_hf_bioreg_ar_df)

rm(parameters)

#4. POPULATION STRUCTURE

#pairwise Fst




---------------------------------------------------------------------------

#5. POPULATION DIFFERENTIATION

#install and load radiator genomic converter
library(radiator)

?genomic_converter

genomic_converter(data = gl_f7_loc_cr,
                  strata = NULL,
                  filename = "gl_f7_arlequin.csv",
                  verbose = TRUE)

genomic_converter(data = gi_ss_gl_f7, strata = NULL, output = "structure")

#convert genlight object to arlequin format using radiator genomic converter
test <- genomic_converter(data = gi_f7_loc_cr, strata = NULL, output = "arlequin")
genomic_converter(data = gi_ss_gl_f7, strata = NULL, output = "arlequin")

#5.POPULATION SPLITS

#convert genlight to treemix input file format
gl2treemix(gl_7_loc_cr, outpath = "data")

#plot tree of m3 bootstrap 100 for gl_filter7_loc_cr
source("treemix/plotting_funcs.R")
plot_tree("treemix/m2/gl7m2b100_stem")
plot_resid("treemix/m4/gl7m4b100_stem", "treemix/poporder2")


gl2structure(ss_gl_f7, indNames = NULL, addcolumns = NULL,
             ploidy = 2, exportMarkerNames = TRUE,
             outfile = "gl.str", outpath = ".", v = 5)
gl.write.csv(ss_gl_f7, outfile = "ss_gl_f7.csv", outpath = "interim/")
gl2structure(ss_gl_f7, indNames = NULL, ploidy = 0, exportMarkerNames = TRUE,
             outfile = "interim/ss_gl_f7.str", outpath = "interim/")
gi_ss_gl_f7_df <- genind2df(gi_ss_gl_f7, pop = NULL, usepop = FALSE, oneColPerAll = FALSE)
View(gi_ss_gl_f7_df)
rm(gi_ss_gl_f7_df)
----------------------------------------------------------------------
View(ss_gl_f7)

gi_ss_gl_f7_gi2df <- genind2df(gi_ss_gl_f7, pop = TRUE,
                               sep = " ", usepop = TRUE,
                               oneColPerAll = TRUE)
rm(gi_ss_gl_f7_gi2df)

pca1 <- glPca(ss_gl_f7)
pca1
#glPca object can be displayed using scatter (a scatterplot of principal components)
scatter(pca1, posi = "bottomright")
title("PCA of sea spurge populaitons in Australia")
library(ape)
tre <- nj(dist(as.matrix(ss_gl_f7)))
tre
plot(tre, typ = "fan", cex = 0.7)
title("NJ tree of the Aus data")

myCol <- colorplot(pca1$scores, pca1$scores, transp = TRUE, cex = 4)
abline(h = 0, v = 0, col = "grey")
add.scatter.eig(pca1$eig[1:40], 2, 1, 2, posi = "topright", inset = .05, ratio = .3)

plot(tre, typ = "fan", show.tip = FALSE)
tiplabels(pch = 20, col = myCol, cex = 4)
title("NJ tree of the Aus sea spurge data")

#DAPC PCA
dapc1 <- dapc(ss_gl_f7, n.pca = 10, n.da = 1)
scatter(dapc1, scree.da = FALSE, bg = "white", posi.pca = "topright",
        legend = TRUE, txt.leg = paste("group", 1:2), col = c("red", "blue"))


#Write summarizeLoci data to .csv for all levels
#population level
write.csv(gt_ss_gl_f7_allsum,
          file = "results/gt_ss_f7_allsum.csv",
          row.names = TRUE, col.names = TRUE)
#state level
write.csv(gt_ss_gl_f7_state_allsum,
          file = "interim/gt_ss_gl_f7_state_allsum.csv",
          row.names = TRUE, col.names = TRUE)


