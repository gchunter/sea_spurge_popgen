#Load packages
library(dartR)
library(hierfstat)
library(adegenet)
library(tidyverse)
library(ggplot2)
library(strataG)

#1. SAMPLE SNP DATA

#Assign dartseq file and ind file to a genlight object
gl_all_samples <- gl.read.dart(filename = "data/Report_DEup19-4268_1_moreOrders_SNP_1.csv",
                   ind.metafile = "data/sea_spurge_ind2_metrics.csv")

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
gl_f1_callrate <- gl.filter.callrate(gl_all_samples, method = "loc",
                                          threshold = 0.90, recalc = TRUE, plot = TRUE)

#filter on departure of Hardy-Weinberg equilibrium for every loci
gl_f2_hwe <- gl.filter.hwe(gl_f1_callrate,
                                alpha = 0.05, basis = "any",
                                bon = TRUE, v = 5)

#report the minor allele frequencies for filter2_HWE
gl.report.maf(gl_f2_hwe, v = 5)

#filter on MAF at 0.01
gl_f3_maf <- gl.filter.maf(gl_f2_hwe, threshold = 0.01, v = 5)

#report the reproducibility for the loci
gl.report.repavg(gl_f3_maf)

#filter out loci for thich the repeateability is less than 0.98
gl_f4_repavg <- gl.filter.repavg(gl_f3_maf, threshold = 0.98, v = 3)

#report the number of monomorphs
gl.report.monomorphs(gl_f4_repavg)

#report secondaries
gl.report.secondaries(gl_f4_repavg)

#filter out secondaries
gl_f5_second <- gl.filter.secondaries(gl_f4_repavg, method = "random", v = 5)

#produce a smearplot
glPlot(gl_f5_second)

#report call rate for individuals
gl.report.callrate(gl_f5_second, method = "ind", plot = TRUE, v = 2)

#filter call rate for individuals at 95%
gl_f6_ind_cr <- gl.filter.callrate(gl_f5_second,
                                        method = "ind", threshold = 0.95,
                                        recalc = TRUE, plot = TRUE)
gl_f7_loc_cr <- gl.filter.callrate(gl_f6_ind_cr,
                                        method = "loc", threshold = 0.95, recalc = TRUE, 
                                        mono.rm = TRUE, plot = TRUE, v = 2)

#population names in gl_f7_loc_cr
popNames(gl_f7_loc_cr)

#3. GENETIC DIVERSITY

#change from genlight to genind
gi_f7_loc_cr <- gl2gi(gl_f7_loc_cr, v =1)

#convert genind to gtypes
gi_f7_g <- genind2gtypes(gi_f7_loc_cr)

#calculate all summary stats for each stratum(population) seperately
allsum_gi_f7_g <- summarizeLoci(gi_f7_g, by.strata = TRUE)
write.csv(allsum_gi_f7_g,
          file = "results/allsum_gi_f7_g.csv",
          row.names = TRUE, col.names = TRUE)

#Calculate the number of private alleles in each strata and locus
pa_gi_f7_g <- privateAlleles(gi_f7_g)

#write a csv of the private alleles object
write.csv(pa_gi_f7_g, file = "results/pa_gi_f7_g.csv",
          row.names = TRUE, col.names = TRUE)

#convert genind to hierfstat format
hf_f7 <- genind2hierfstat(gi_f7_loc_cr, pop = NULL)
attach(hf_filter7)
hf_filter7_loci <- hf_filter7[,c(2:752)]
hf_filter7_pop <- hf_filter7$pop
hf_filter7_basic_stats <- basic.stats(data.frame(hf_filter7_pop, 
                                                 hf_filter7_loci),
                                      diploid = TRUE, digits = 4)

#Calculate basic.stats for hf_filter7 and assign to object
hf_f7_bs <- basic.stats(hf_filter7, diploid = TRUE, digits = 4)

#allelic richness
hf_f7_ar <- allelic.richness(hf_filter7, min.n = NULL, diploid = TRUE)

#View allelic richness results
View(hf_f7_ar[["Ar"]])

#make a table of Ar values from hf_filter7
table_hf_f7_ar <- hf_f7_ar[["Ar"]]

#write a .csv file of table_hf_f7_ar
write.csv(table_hf_f7_ar,
          file = "results/table_hf_f7_ar.csv",
          append = FALSE, quote = FALSE, row.names = TRUE,
          col.names = TRUE)

#Make a table of Hs values only from hf_f7_bs
table_hf_f7_bs_Hs <- hf_f7_bs[["Hs"]]

#write a .csv file of table_hf_f7_bs_Hs
write.csv(table_hf_f7_bs_Hs,
          file = "results/table_hf_f7_bs_Hs.csv",
          append = FALSE, quote = FALSE, row.names = TRUE,
          col.names = TRUE)

#create a dataframe of Hs from the .csv file
df_hf_f7_bs_Hs <- read_csv("results/table_hf_f7_bs_Hs.csv")

#remove column X1
df2_hf_f7_bs_Hs <- select(df_hf_f7_bs_Hs, -X1)

#View Fis
view(hf_f7_bs[["Fis"]])

#make a table of Fis values only from hf_f7_bs
table_hf_f7_bs_Fis <- hf_f7_bs[["Fis"]]

#write a .csv of the table_hf_f7_bs_Fis
write.csv(table_hf_f7_bs_Fis, file = "results/table_hf_f7_bs_Fis.csv",
          append = FALSE, quote = FALSE, row.names = TRUE,
          col.names = TRUE)

#View Ho
view(hf_f7_bs[["Ho"]])

#Make a table of Ho values only from hf_f7_bs
table_hf_f7_bs_Ho <- hf_f7_bs[["Ho"]]

#Write a .csv of the object table_hf_f7_bs_Ho
write.csv(table_hf_f7_bs_Ho,
          file = "results/table_hf_f7_bs_Ho.csv",
          append = FALSE, quote = FALSE, row.names = TRUE,
          col.names = TRUE)

#calculate allelic richness and print output to csv
hf_filter7_ar <- allelic.richness(data.frame(hf_filter7_pop,
                                             hf_filter7_loci),
                                  min.n = NULL,
                                  diploid = TRUE)

table_ar <- hf_filter7_ar[["Ar"]]

write.csv(table_ar, file = "hierfstat/table_ar.csv", row.names = TRUE,
          col.names = TRUE)

table(hf_filter7[,1])

#calculate basic statistics
hf_filter7_basic_stats <- basic.stats(hf_filter7)

#4. POPULATION DIFFERENTIATION

#install and load radiator genomic converter
library(radiator)
?genomic_converter
genomic_converter(data = gl_f7_loc_cr,
                  strata = NULL,
                  filename = "gl_f7_arlequin.csv",
                  verbose = TRUE)

#convert genlight object to arlequin format using radiator genomic converter
test <- genomic_converter(data = gi_f7_loc_cr, strata = NULL, output = "arlequin")


#5.POPULATION SPLITS

#convert genlight to treemix input file format
gl2treemix(gl_7_loc_cr, outpath = "data")

#plot tree of m3 bootstrap 100 for gl_filter7_loc_cr
source("treemix/plotting_funcs.R")
plot_tree("treemix/m2/gl7m2b100_stem")
plot_resid("treemix/m4/gl7m4b100_stem", "treemix/poporder2")

poporder <- list("Anxious Bay beach Elliston", "Aslings Beach Eden", "Balnarring Beach",
                 "Bawley Point", "Beach near Lake Bunga", "Bennetts Beach Yacaaba Peninsula",
                 "Bernie Central", "Birthday Bay Southwest", "Blossom Beach", "Bmerwerre Beach",                   
                 "Boydtown beach", "Broulee South", "Brown Bay Port MacDonnell",
                 "Camel Rock Beach", "Cape Borks Lighthouse", "Dunes Beach south west",            
                 "Endeovon Bay SouthWest", "Esperence Harbour Beach", 
                 "Greens Beach", "Hankerchief beach","Hardwick Bay", "Hilkier Bay", 
                 "Kilcunda", "Lagoon Beach Lord Howe Island", "Loaders", 
                 "Lone Pine Sleaford", "Lowlands Beach", "Lucky Bay", 
                 "Maria Island", "Myponga beach", "North Bernie", "Pambula Beach", 
                 "Pardoe Northdown", "Pardoe Northdown Conservation Area", 
                 "Peppermint Grove Beach", "Port Latta", "Port Neill Foreshore Beach", 
                 "Princetown Beach", "Pt Lesueur Maria Island", "Riedle Beach",
                 "Rifle Butts Beach Port Victoria", "Silverleaves Phillip Island",
                 "Ski Beach Tumby Bay", "Squeaky Beach", "Stansbury", "Steamers Beach", 
                 "Tea Tree Lane", "Venus Bay", "West Beach Hopetoun", 
                 "Woolamai Beach Phillip Island", "Yanchep")






