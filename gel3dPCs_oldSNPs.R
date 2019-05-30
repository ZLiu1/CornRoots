setwd("~/GEL3D")
gelPC<-read.csv('remainingTraits_PCs_VIF8.csv',header = TRUE)
# genotype
###---------------- Genot matrix file creation--------------------------------------------
ipSNP<-read.table("GAPIT.Genotype.Numerical.txt",header=TRUE) #Choose the GAPIT.Genotype.Numerical.txt,  
#which is transformed from HapMap format and imputed with Major allele in GAPIT
str(ipSNP)
colnames(ipSNP)
#rename the SNPs to SNP1.....1000, et al 
names<-gsub("recgenoMx", "SNP", colnames(ipSNP))
colnames(ipSNP)<-names
#add one column to the data.frame
Ind_id<-rownames(ipSNP)
ipSNP$Ind_id<-Ind_id
ipSNP$Ind_id
rownames(ipSNP)<-NULL #remove the row names
ipSNPf<-ipSNP[,c(ncol(ipSNP),1:(ncol(ipSNP)-1))] #move the Ind_id column to the start
genot<-ipSNPf
# order by 'Ind_id' to make sure the genot file has the same order with phenot file!!!!
genot <- genot[order(genot$Ind_id),] 
#transform data.frame to matrix
genot_mat <- as.matrix(genot[, 2:ncol(genot)])
rownames(genot_mat) <- genot$Ind_id
str(genot_mat)

### ---------------------map file creation-----------------------
map <- read.table("ipMap.txt", sep = "\t", header = TRUE)#choose the ipMap.txt
head(map)
#rename the SNPs to SNP1.....1000, et al, make sure the SNPs have the same name 
names<-gsub("ip", "", map$SNP)
map$SNP<-names
head(map)
###----common genotypes between genot and phenot--------------------------------------------------------------
phenot<-gelPC #option 1
intersect(genot$Ind_id,phenot$Genotype)-> A #common genotypes/Ind_id shared by genot and phenot
setdiff(genot$Ind_id,phenot$Genotype) #Ind_id not shared 
genot_com<-genot_mat[rownames(genot_mat) %in% A,] #genotypic data only with shared Ind_id
class(genot_com)
# rename genot_com to genot_mat to be convenient
genot_mat<-genot_com
# the common phenotypic data
phenot_f<-phenot[which(phenot$Genotype %in% A =='TRUE'),]
phenot_f<-phenot_f[order(phenot_f$Genotype),] 
# just in case there were NAs, replace NAs with column mean
for(i in 2:ncol(phenot_f)){
  phenot_f[is.na(phenot_f[,i]), i] <- mean(phenot_f[,i], na.rm = TRUE)
}
write.csv(phenot_f,'phenotypic_genotypedPCwithBiomass.csv',row.names = FALSE) #only for the first time

##########--------could read in the 
##-------A relationship (K) matrix should be estimated from the genotype data,
##-------using e.g. the following code:
average <- colMeans(genot_mat, na.rm = T)
stdev <- apply(genot_mat, 2, sd)
genot_stand <- sweep(sweep(genot_mat, 2, average, "-"), 2, stdev, "/")
K_mat <- (genot_stand %*% t(genot_stand)) / ncol(genot_stand)
str(K_mat)
K_mat[1:5, 1:5]
ncol(K_mat)
length(phenot_f)

####---------RUNNIING MLMM------------------------------################
source("mlmm.r")
source("emma.r")
source("plot_mlmm.r")
#---- perform mlmm (20 steps), it can take few minutes...
#---thresh = user define threshold, delete this if no user defined threshold

###--- final loop for all traits including rPCs and PCs and also output the p value and TAS effect--------------------------------------
for (i in 2:length(phenot_f)) {
  mygwas <- mlmm(Y = phenot_f[,i], X = genot_mat, K = K_mat, nbchunks = 2, maxsteps = 21,thresh = 1.0 * 10^-5) # Not 2.0*10-6 this time!!!
  #optimal step according to user defined threshold
  setwd('~/GEL3D/plot')
  pdf(paste("mPlot", colnames(phenot_f)[i], ".pdf", sep = ""),height = 4,width = 10)
  plot_opt_GWAS(mygwas, opt = "thresh", snp_info = map, pval_filt = 0.1)#dotted line correspond to the user defined threshold
  qqplot_opt_GWAS(mygwas, opt = "thresh") 
  #% variance plot
  plot_step_RSS(mygwas)
  dev.off()
  #opt thresh--user difined threshold
  setwd('~/GEL3D/optimal')
  head(mygwas$opt_thresh$out)
  opt_thresh_p<-subset(merge(map,mygwas$opt_thresh$out,by='SNP'),pval<=1e-3)
  write.csv(opt_thresh_p,paste("opt_thresh_p_", colnames(phenot_f)[i], ".csv", sep = ""),row.names=FALSE)
  #including SNP effects
  setwd('~/GEL3D/opt_eff')
  mygwas$opt_thresh$coef
  effects_SNP<-mygwas$opt_thresh$coef
  write.csv(effects_SNP,paste("Effects_SNP_", colnames(phenot_f)[i], ".csv", sep = ""),row.names = FALSE)
  # maximum Cofactor model
  setwd('~/GEL3D/maxModel')
  maxCofModel <- subset(merge(map,mygwas$pval_step[[length(mygwas$pval_step)]]$out,by='SNP'),pval<=1e-3)
  write.csv(maxCofModel,paste("EmaxCofModel", colnames(phenot_f)[i], ".csv", sep = ""),row.names = FALSE)
}

##############----$$$$$$$$$$___merge GWAS results into one file--------------++++++++++#############
setwd ("~/GEL3D/maxModel/")
##---Loading multiple .csv files as separate data frames
folder <-"~/GEL3D/maxModel/"  # path to folder that holds multiple .csv files
file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder

trait<- gsub('.csv', '',file_list) #simplify the file name
trait<-gsub('EmaxCofModel','',trait) # simplify the file name to be only variable names
##---read in each .csv file in file_list and create a data frame with the triat name/variable as the .csv file
for (i in 1:length(file_list)){
  assign(paste('df',trait[i],sep = '_'), #change to file_list[i] if you need the same file name as dataframe name
         read.csv(paste(folder, file_list[i], sep=''))
  )}


##----create a list of all the data frames/files
# Loads all data.frames in current environment that begins with df 
# followed by 1 to any amount of numbers
mylist<-lapply(ls(pattern="df_"), function(x) get(x))
##----add "Trait" column to each data frame
for (i in (1:length(mylist))){
  mylist[[i]]$Trait<-rep(trait[i],nlevels(mylist[[i]]$SNP))
}

##----rbind all dataframes into one
library(dplyr)
data<-bind_rows(mylist)
extr0.05<-data[which(data$pval <= 1.0*10^-5),]
#extr_r0.05<-data[which(data$pval <= 2.0*10^-6),]
# add physical positions-----no need in this case
#extr_r0.05<-merge(extr_r0.05, map, by='SNP', all.x = TRUE)
##------write the extracted .csv files
write.csv(extr0.05,'extr_e-5maxModel.csv',row.names = FALSE)
#write.csv(extr_r0.05,'extracted results_adjustedBfe-6maxModel.csv',row.names = FALSE)

###------optimal model--------------########
setwd ("~/GEL3D/optimal/")
##---Loading multiple .csv files as separate data frames
folder <-"~/GEL3D/optimal/"  # path to folder that holds multiple .csv files
file_list <- list.files(path=folder, pattern="*.csv") # create list of all .csv files in folder

trait<- gsub('.csv', '',file_list) #simplify the file name
trait<-gsub('opt_thresh_p_','',trait) # simplify the file name to be only variable names
##---read in each .csv file in file_list and create a data frame with the triat name/variable as the .csv file
for (i in 1:length(file_list)){
  assign(paste('df',trait[i],sep = '_'), #change to file_list[i] if you need the same file name as dataframe name
         read.csv(paste(folder, file_list[i], sep=''))
  )}


##----create a list of all the data frames/files
# load Loads all data.frames in current environment that begins with df 
# followed by 1 to any amount of numbers
mylist<-lapply(ls(pattern="df_"), function(x) get(x))

##----add "Trait" column to each data frame
for (i in (1:length(mylist))){
  mylist[[i]]$Trait<-rep(trait[i],nlevels(mylist[[i]]$SNP))
}

##----rbind all dataframes into one
library(dplyr)
data<-bind_rows(mylist)
extr0.05<-data[which(data$pval < 1.0*10^-5),]
#extr_r0.05<-data[which(data$pval <= 2.0*10^-6),]
# add physical positions-----no need in this case
#extr_r0.05<-merge(extr_r0.05, map, by='SNP', all.x = TRUE)
##------write the extracted .csv files
write.csv(extr0.05,'extr_e-5optimal.csv',row.names = FALSE)
#write.csv(extr_r0.05,'extracted results_adjustedBfe-6optimal.csv',row.names = FALSE)
file_list
trait
