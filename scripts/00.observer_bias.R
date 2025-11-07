
# Read in cleaned up anonymized phenotype data FULL
data=read.csv(file="rawdata/phenotype-data-anon_full.csv",stringsAsFactors = F)

### Define haplotype groups
source("scripts/01.haplotype-groups.R")

# Cleanup backcross form data (F1 cross info)
source("scripts/02.backcross-cleanup.R")

# Merge phenotype data with backcross/treatment data
source("scripts/03.merge-data.R")

merged$num_males=merged$X0000 + merged$X0001 + merged$X0010 + merged$X0011 + merged$X0100 + merged$X0101 + merged$X0110 + merged$X0111 + merged$X1000 + merged$X1001 + merged$X1010 + merged$X1011 + merged$X1100 + merged$X1101 + merged$X1110 + merged$X1111

#Check observer triple crossover rates
obs_check=aggregate(merged[, c(haplotypes, "num_males")], by= list(merged$observer), sum)
obs_check$TCO_rate=(obs_check$X1010+obs_check$X0101)/obs_check$num_males
hist(obs_check$TCO_rate)

#TCO rate estimates
mean_tco=mean(obs_check$TCO_rate)
sd_tco=sd(obs_check$TCO_rate)
md=mad(obs_check$TCO_rate)
iqr=IQR(obs_check$TCO_rate,na.rm = TRUE)
md

# Expected TCO rate across dataset using Comeron Map Distances
exp_tco=(0.12 * 0.21 * 0.24)

# Observers removed from dataset were:
## Individuals #7 and #18 who did not observe >50 males in the total dataset
## Individuals #5 and # 19 who has an excess of TCO events as compared to the rest of the dataset and far beyond expectation

# difference from exp TCO
obs_check$TCO_rate[obs_check$Group.1=="19"]/exp_tco
obs_check$TCO_rate[obs_check$Group.1=="5"]/exp_tco
#Both observers were over nearly an order of magnitude above the Expected TCO rate

# Mean TCO without removed observers:
exp_tco
mean_tco
adj_mean=mean(obs_check$TCO_rate[obs_check$Group.1!="5" & obs_check$Group.1!="19" & obs_check$Group.1!="7" & obs_check$Group.1!="18"])
adj_mean

# Make plot for supplement: 
tco_check=ggplot(data=obs_check,aes(x=Group.1,y=TCO_rate))+geom_point()+geom_hline(yintercept = mean_tco ,col="red")+geom_hline(yintercept = exp_tco,col="purple")+geom_hline(yintercept = mean_tco+md,col="grey50")+geom_hline(yintercept = mean_tco+iqr,col="grey30")+geom_hline(yintercept = adj_mean,col="orange")+xlab("Anonymized observer ID")+ylab("Triple Crossover Rate")+geom_text(aes(label=ifelse(TCO_rate>mean_tco+md,as.character(Group.1),ifelse(num_males<50,as.character(Group.1),''))),hjust=0,vjust=0)
tco_check
ggsave("images/Observer_Bias_TCO.png")


# Check for happlotype skews
source("scripts/04.haplotype-bias.R")

# Recombination rate comparison (not doing Kosambi correction due to LARGE number of crossovers that leads to NaN values produced)

recomb <- by_vialday

recomb$nco_inds <- rowSums(recomb[names(recomb) %in% haps_nco], )
recomb$sco_inds <- rowSums(recomb[names(recomb) %in% haps_sco], )
recomb$dco_inds <- rowSums(recomb[names(recomb) %in% haps_dco], )
recomb$tco_inds <- rowSums(recomb[names(recomb) %in% haps_tco], )
recomb$co_inds <- recomb$sco_inds + 2*recomb$dco_inds + 3*recomb$tco_inds
recomb$num_M <- rowSums(recomb[names(recomb) %in% haplotypes], )

recomb$recomb_rate <- recomb$co_inds / recomb$num_M
recomb <- na.omit(recomb)

## Rough recombination rate by interval for validation purposes
recomb$ycv_count <-rowSums(recomb[names(recomb) %in% ycv_haps], )
recomb$cvv_count <-rowSums(recomb[names(recomb) %in% cvv_haps], )
recomb$vf_count <-rowSums(recomb[names(recomb) %in% vf_haps], )

recomb$ycv_rr <- recomb$ycv_count / recomb$num_M
recomb$cvv_rr <-recomb$cvv_count / recomb$num_M
recomb$vf_rr <-recomb$vf_count / recomb$num_M

total_ycv_count <- sum(recomb[names(recomb) %in% ycv_haps])
total_cvv_count <- sum(recomb[names(recomb) %in% cvv_haps])
total_vf_count <- sum(recomb[names(recomb) %in% vf_haps])
total_num_M <- sum(recomb["num_M"])

total_ycv_rr <- total_ycv_count / total_num_M
total_cvv_rr <-total_cvv_count / total_num_M
total_vf_rr <-total_vf_count / total_num_M

## Remove haplotype columns
recomb <- recomb[,!names(recomb) %in% haplotypes]

recomb_rate=summaryBy(co_inds+ycv_count+cvv_count+vf_count+num_M~MaternalVial+vial_letter+PaternalStock+Treatment,data=recomb,FUN=sum,na.rm=T)
recomb_rate=recomb_rate[recomb_rate$num_M.sum>=5,]
recomb_rate$ycv_rr=recomb_rate$ycv_count.sum/recomb_rate$num_M.sum
recomb_rate$cvv_rr=recomb_rate$cvv_count.sum/recomb_rate$num_M.sum
recomb_rate$vf_rr=recomb_rate$vf_count.sum/recomb_rate$num_M.sum
recomb_rate$recomb_rate=rowSums(recomb_rate[,c("ycv_rr","cvv_rr","vf_rr")],na.rm=TRUE)
recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$PaternalStock), mean)
names(recomb_summary)[1:2] <- c("Treatment", "PaternalStock")

recomb_summary$treat=ifelse(recomb_summary$Treatment=="2x",2,ifelse(recomb_summary$Treatment=="1x",1,0.5))

write.csv(recomb_summary, "output/recombination_summary_FULL.csv", row.names = FALSE)


# Levene's Test for Homogeneity of Variance across dataset:
leveneTest(num_offspring~interaction(PaternalStock, Treatment, vial_letter),recomb)
# p=0.03903
# Significant heterogeneity of variance detected!!!


#Similar but this time across observers:
leveneTest(num_males~as.factor(observer),merged)
#p< 2.2e-16




