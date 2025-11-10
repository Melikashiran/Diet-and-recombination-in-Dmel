
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
#recomb <- recomb[,!names(recomb) %in% haplotypes]

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


### Crossover interference calculation

#COI <- recomb
COI=summaryBy(co_inds+ycv_count+cvv_count+vf_count+num_M+X0100 + X1011+X0010 + X1101+X0110 + X1001~PaternalStock+Treatment+vial_letter,data=recomb,FUN=sum,na.rm=T)
COI$COI=COI$co_inds.sum/COI$num_M.sum
COI$ycv_rr=COI$ycv_count.sum/COI$num_M.sum
COI$cvv_rr=COI$cvv_count.sum/COI$num_M.sum
COI$vf_rr=COI$vf_count.sum/COI$num_M.sum
COI=COI[COI$num_M.sum>=5,]

# expected vs observed
COI$Exp_DCO_ycv_cvv <-  round((COI$ycv_rr * COI$cvv_rr)*COI$num_M,0)
COI$Obs_DCO_ycv_cvv <-  COI$X0100 + COI$X1011

COI$Exp_DCO_cvv_vf <-  round((COI$cvv_rr * COI$vf_rr)*COI$num_M,0)
COI$Obs_DCO_cvv_vf <-  COI$X0010 + COI$X1101

COI$Exp_DCO_ycv_vf <-  round((COI$ycv_rr * COI$vf_rr)*COI$num_M,0)
COI$Obs_DCO_ycv_vf <-  COI$X0110 + COI$X1001


# coefficient of coincidence and interference
COI$COC_ycv_cvv <- COI$Obs_DCO_ycv_cvv / COI$Exp_DCO_ycv_cvv
COI$Interference_ycv_cvv <- 1 - COI$COC_ycv_cvv

COI$COC_cvv_vf <- COI$Obs_DCO_cvv_vf / COI$Exp_DCO_cvv_vf
COI$Interference_cvv_vf <- 1 - COI$COC_cvv_vf

COI$COC_ycv_vf <- COI$Obs_DCO_ycv_vf / COI$Exp_DCO_ycv_vf
COI$Interference_ycv_vf <- 1 - COI$COC_ycv_vf


# Plotting interference against genetic map distance

# Need to reorganize the data frame: Interference measure and Genetic Distance Measure
wide_df=COI[,c("PaternalStock","Treatment","vial_letter","Interference_cvv_vf","Interference_ycv_cvv","Interference_ycv_vf")]
colnames(wide_df)=c("PaternalStock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df=melt(wide_df,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Interference")

wide_df2=COI[,c("PaternalStock","Treatment","vial_letter","ycv_count.sum","cvv_count.sum","vf_count.sum","co_inds.sum")]
wide_df2$ycv_cvv=rowSums(wide_df2[,c("ycv_count.sum","cvv_count.sum")],na.rm=TRUE)
wide_df2$cvv_vf=rowSums(wide_df2[,c("cvv_count.sum","vf_count.sum")],na.rm=TRUE)
wide_df2$ycv_vf=wide_df2$co_inds.sum
wide_df2=wide_df2[,c("PaternalStock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")]
long_df2=melt(wide_df2,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Num_COs")

wide_df3=COI[,c("PaternalStock","Treatment","vial_letter","num_M.sum")]
wide_df3$sub1=wide_df3$num_M.sum
wide_df3$sub2=wide_df3$num_M.sum
long_df3=melt(wide_df3,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Total")

wide_df5=COI[,c("PaternalStock","Treatment","vial_letter","Exp_DCO_cvv_vf","Exp_DCO_ycv_cvv","Exp_DCO_ycv_vf")]
colnames(wide_df5)=c("PaternalStock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df5=melt(wide_df5,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Exp_DCO")

wide_df6=COI[,c("PaternalStock","Treatment","vial_letter","Obs_DCO_cvv_vf","Obs_DCO_ycv_cvv","Obs_DCO_ycv_vf")]
colnames(wide_df6)=c("PaternalStock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df6=melt(wide_df6,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Obs_DCO")

long_df4=cbind(long_df,long_df2[,5],long_df3[,5],long_df5[,5],long_df6[,5])
colnames(long_df4)=c("PaternalStock","Treatment","vial_letter","Interval","Interference","Num_COs","Total","Exp_DCO","Obs_DCO")
long_df4$Gendist=long_df4$Num_COs/long_df4$Total
long_df4$Gendist=ifelse(long_df4$Num_COs/long_df4$Total<0.5,long_df4$Num_COs/long_df4$Total,0.499+(((long_df4$Num_COs/long_df4$Total)-0.5)/100))
#long_df4$kosambi_distances <- (100/4)*(log((1+(2*long_df4$Gendist))/(1-(2*long_df4$Gendist))))
long_df4$PaternalStock=as.factor(long_df4$PaternalStock)

ggplot(data=long_df4,aes(x=Gendist,y=Interference,col=Treatment))+geom_line()+facet_wrap(~PaternalStock)



COI_final=summaryBy(Num_COs+ Total+ Exp_DCO +Obs_DCO~PaternalStock+Treatment+Interval,data=long_df4,FUN=sum,na.rm=T)

COI_final$COC <- COI_final$Obs_DCO.sum / COI_final$Exp_DCO.sum
COI_final$Interference <- 1 - COI_final$COC
COI_final$Gendist=ifelse(COI_final$Num_COs/COI_final$Total<0.5,COI_final$Num_COs/COI_final$Total,0.499)
#COI_final$kosambi_distances <- (100/4)*(log((1+(2*COI_final$Gendist))/(1-(2*COI_final$Gendist))))

COI_final=summaryBy(Interference+Gendist~PaternalStock+Treatment+Interval,data=long_df4,FUN=mean,na.rm=T)
col_pal2=c("#fdd0a2","#fd8d3c","#a63603","#c6dbef","#6baed6","#08519c")

diet_gxe_int=ggplot(data=COI_final,aes(x=Gendist.mean,y=Interference.mean,col=Treatment))+geom_line()+facet_wrap(~PaternalStock) + xlab("Recombination Rate (Percent)") + ylab("Crossover Interference") +   #ggtitle("Interference versus Recombination Rate") +
  theme_base() + scale_color_manual(values = col_pal2[1:3])
#ylim(0.1,1) + xlim(0.15,0.52)
diet_gxe_int
ggsave("images/COI-vs-RR_42_FULL.png", height=5)

diet_gxe_int=ggplot(data=COI_final,aes(x=Gendist.mean,y=Interference.mean,col=Treatment))+geom_line()+facet_wrap(~PaternalStock) + xlab("Recombination Rate (Percent)") + ylab("Crossover Interference") +   #ggtitle("Interference versus Recombination Rate") +
  theme_base() + scale_color_manual(values = col_pal2[4:6])
#ylim(0.1,1) + xlim(0.15,0.52)
diet_gxe_int
ggsave("images/COI-vs-RR_217_FULL.png", height=5)

