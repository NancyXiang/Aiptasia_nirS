
############################################################
##################### Alpha-diversity ######################
############################################################

setwd("~/WP3_Aiptasia_nirS/nirS_analysis")

install.packages("gridExtra")
install.packages("plyr")
install.packages("vegan")
install.packages('ggpubr')

library(vegan)
library(ggplot2)
library(plyr)
library(gridExtra)
library(ggpubr)

col_Symb=c("#CECACC", "#E0A800", "#006C67")
asv=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta.txt", header = TRUE, sep ='\t')
map=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_input/Aip_nirS_metadata.txt", header = TRUE, row.names = 1, sep ='\t')

### rarefying ###########
cnts=t(asv[, 1:43]) # define samples
cnts=subset(cnts, rownames(cnts) != "CA3") # remove CA3,  samples without DNA input
cnts=subset(cnts, rownames(cnts) != "CA5") # remove CA5,  samples without DNA input
rownames(cnts)

min(rowSums(cnts)) #12521 reads
asv.rar=data.frame(rarefy(cnts, 12521)) # min reads for one sample
colnames(asv.rar)[1] <- "otu.tab.rff"
asv.rar=merge(asv.rar, cnts, by="row.names")
rownames(asv.rar)= asv.rar[,1]
asv.rar=asv.rar[,-1]

asv.grp=merge(asv.rar, map, by="row.names")
rownames(asv.grp)= asv.grp[,1]
asv.grp=asv.grp[,-1]

#############################
########  Prep Plot  ########
#############################

alpha=as.data.frame(t(estimateR(asv.grp[, 2:428],  smallsample = TRUE))) # define the location of ASV
alpha$Shannon=diversity(asv.grp[, 2:428], index = "shannon") # shannon
alpha$simpson=diversity(asv.grp[, 2:428], index = "simpson") # simpson
alpha$invsimpson=diversity(asv.grp[, 2:428], index = "invsimpson") # invsimpson

alpha$host=map$host[match(rownames(alpha),rownames(map))]
alpha$symbiont=map$symbiont[match(rownames(alpha),rownames(map))]
alpha$other=map$other[match(rownames(alpha),rownames(map))]
alpha[is.na(alpha)] <- "None"

alpha$host=factor(alpha$host, levels = c( "None","CC7", "H2"))
alpha$symbiont=factor(alpha$symbiont, levels = c("None","SSA01", "SSB01"))
alpha$symbiont=gsub("None", "APO", alpha$symbiont)
alpha$combine=paste(alpha$host, alpha$symbiont, sep = "_")

#####################################################
##################### boxplot  ######################
#####################################################

# Chao1 richness, Shannon diversity, Simpson evenness
#Chao1 richness
chao=ggplot(alpha, aes(x=symbiont, y=S.chao1, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=col_Symb)  +
  facet_grid(~host) +
  theme_classic() + labs( y= "Chao1 richness", x="")

#Shannon diversity
shannon=ggplot(alpha, aes(x=symbiont, y=Shannon, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=col_Symb)  +
  facet_grid(~host) +
  ylim(0, 2.0)+
  theme_classic() + labs( y= "Shannon diversity", x="")

# Simpson evenness
simpson=ggplot(alpha, aes(x=symbiont, y=simpson, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=col_Symb)  +#change color question
  ylim(0, 0.8)+
  facet_grid(~host) +
  theme_classic() + labs( y= "Simpson evenness", x="")

# output six combs
pdf("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_output/new_alphaDiversity_six.pdf", width=6,height=8, pointsize = 12)
ggarrange(chao, shannon, simpson + rremove("x.text"),
          #legend = "none",
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()


################### subset for six combinations #####################
alpha_six=subset(alpha, host == "CC7"| host == "H2")
alpha_six$symbiont=gsub("None", "APO", alpha_six$symbiont)
#Chao1 richness
# You need to create a new variable that has host + symbiont and use that one there.
alpha_six$combine=paste(alpha_six$host, alpha_six$symbiont, sep = "_")
ggplot(alpha_six, aes(x=host, y=S.chao1, fill=combine)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=c("#AF81C9", "#F89A7E", "#F2CA85","#53D1F1","#6E87BB","#4EB398"))  +
  #facet_grid(~host) +
  theme_classic() + labs( y= "Chao1 richness", x="")

sixChao=ggplot(alpha_six, aes(x=host, y=S.chao1, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=c("#AF81C9", "#F89A7E", "#F2CA85"))  +
  #facet_grid(~host) +
  theme_classic() + labs( y= "Chao1 richness", x="")

#Shannon diversity
sixShan=ggplot(alpha_six, aes(x=host, y=Shannon, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=c("#AF81C9", "#F89A7E", "#F2CA85"))  +
  #facet_grid(~symbiont) +
  ylim(0, 2.0)+
  theme_classic() + labs( y= "Shannon diversity", x="")

# Simpson evenness
sixSimp=ggplot(alpha_six, aes(x=host, y=simpson, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=c("#53D1F1","#6E87BB","#4EB398"))  +#change color question
  ylim(0, 0.8)+
  #facet_grid(~symbiont) +
  theme_classic() + labs( y= "Simpson evenness", x="")

# output six combs
pdf("WP3_Aiptasia_nirS/nirS_analysis/R_output/alphaDiversity_six.pdf", width=3,height=8.5, pointsize = 12)
ggarrange(sixChao, sixShan, sixSimp + rremove("x.text"),
          legend = "none",
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()


#################### subset for algae ##########################
# Chao1
Sym_Chao=ggplot(alpha_Sym, aes(x=symbiont, y=S.chao1, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=c("#BCABD0", "#A8CE9F"))  +
  ylim(0, 40)+
  # facet_grid(~symbiont) +
  theme_classic() + labs( y= "Chao1 richness", x="")
# Shannon
Sym_Shan=ggplot(alpha_Sym, aes(x=host, y=Shannon, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  scale_fill_manual(values=c("#BCABD0", "#A8CE9F"))  +
  #facet_grid(~symbiont) +
  ylim(0, 2.0)+
  theme_classic() + labs( y= "Shannon diversity", x="")
# Simpson evenness
Sym_Simp=ggplot(alpha_Sym, aes(x=host, y=simpson, fill=symbiont)) +
  stat_boxplot(geom = "errorbar")  +
  geom_boxplot(alpha = 1) +
  ylim(0, 0.8)+
  scale_fill_manual(values=c("#BCABD0", "#A8CE9F"))  +
  #facet_grid(~symbiont) +
  theme_classic() + labs( y= "Simpson evenness", x="")

# output algae
pdf("WP3_Aiptasia_nirS/nirS_analysis/R_output/alphaDiversity_algae.pdf", width=2,height=8.5, pointsize = 12)
ggarrange(Sym_Chao, Sym_Shan, Sym_Simp + rremove("x.text"),
          legend = "none",
          #labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()


##################################################
##################### Stats ######################
##################################################

########### six combinations using two-way ANOVA ###########
shapiro.test(alpha$S.chao1) # P=0.6392, p-value > 0.05 implying we can assume the normality.
shapiro.test(alpha$Shannon) # P=0.598,  p-value > 0.05 implying we can assume the normality.
shapiro.test(alpha$simpson) # P=0.1423,  p-value > 0.05 implying we can assume the normality.

### chao1
chao_AN <- aov( S.chao1 ~ host * symbiont, data=alpha_six)
#plot(fitted(chao_AN), residuals(chao_AN))
library(nortest)
ad.test(residuals(chao_AN)) # P=0.5345, should >0,05, Anderson-Darling test for normality
cvm.test(residuals(chao_AN)) # P=0,6758, should >0,05, Cramer-Von Mises Test for normal distribution
shapiro.test(residuals(chao_AN)) # P=0,3243, should > 0,05, Normality Test
summary(chao_AN)
TukeyHSD(chao_AN)

### Shannon
Shannon_AN <- aov( Shannon ~ host * symbiont, data=alpha_six)
ad.test(residuals(Shannon_AN)) # P=0.08
cvm.test(residuals(Shannon_AN)) # P=0,09
shapiro.test(residuals(Shannon_AN)) # P=0,08
summary(Shannon_AN)
TukeyHSD(Shannon_AN)

### simpson
simpson_AN <- aov(log10(simpson) ~ host * symbiont, data=alpha_six) #log transformation to meet normality
ad.test(residuals(simpson_AN)) # P=0.2332
cvm.test(residuals(simpson_AN)) # P=0,227
shapiro.test(residuals(simpson_AN)) # P=0,2131
summary(simpson_AN)
TukeyHSD(simpson_AN)


########### Symbiont strains using t-test ###########
### chao1
t.test(S.chao1~ symbiont, alpha_Sym)# p-value = 0.296
### Shannon
t.test(Shannon~ symbiont, alpha_Sym)#p-value = 0.8011
### simpson
t.test(simpson~ symbiont, alpha_Sym)#p-value = 0.7664


########## finish 02/28/2022 ################################
