##########
# for plotting the collated alpha diversity from the rarefactions
# did rarefactions from 100 to 1000, by steps of 25, each run 100 times

library(ggplot2)


#################
#shannon
#################
## read the data for rarefied shannon diversity
alph.s<-read.delim("~/shannon.txt",
                   na.strings = 'n/a', header = TRUE)
names(alph.s)

# calculate the means and SEs for each sample at each sampling depth:
m<-c() # makes an empty vector for the results
for(i in 4:81){#column numbers with samples (check this for your file)
  mm<-tapply(alph.s[,i], factor(alph.s$sequences.per.sample), mean, na.rm=TRUE) 
  m<-c(m, mm)
}

ses<-c()
for(i in 4:81){
  ss<-tapply(alph.s[,i], factor(alph.s$sequences.per.sample), se, na.rm=TRUE)
  ses<-c(ses, ss)
}

# put the means and SEs into a data frame for easier plotting
# figure out the number of sampling depths:
length(levels(factor(alph.s$sequences.per.sample))) 
# repeat the names of the samples the number of times list above (in this case 37)
# include the number of sequences sampled
pdata<-data.frame(id=rep(names(alph.s[,c(4:81)]), each=37), sequences=seq(100,1000,25))
# add the means and the SEs:
pdata2<-cbind(pdata, m, ses)

# calculate 95% confidence intervals for each mean estimate of diversity:
pdata2 <- within(pdata2, {
  LL <- m - 1.96 * ses
  UL <- m + 1.96 * ses
})

# check the data:
head(pdata2)
tail(pdata2)

#make 2 figures each with half the samples
dim(pdata2) # 2886    6
# to split the file in half
2886/2 # 1443

# these can be plotted with varying y-axes (1) or all on the same axis (2)
# 1:
quartz(height=8, width=14)
ggplot(data=pdata2[1:1443, ], aes(x=factor(sequences), y=m)) +  
  geom_errorbar(aes(x=factor(sequences), ymin=LL, ymax=UL), col="black", width=0.1) + 
  geom_point(aes(x=factor(sequences), y=m), size=2) + 
  geom_line() + labs(x="number of sequences", y="shannon diversity") + theme_bw() + 
  facet_wrap(~id, nrow = 7, ncol = 6, scales = "free_y") + 
  theme(axis.text.x = element_text(size = 5.5, angle = 90))

# 2:
quartz(height=8, width=14)
ggplot(data=pdata2[1:1443, ], aes(x=factor(sequences), y=m)) +  
  geom_errorbar(aes(x=factor(sequences), ymin=LL, ymax=UL), col="black", width=0.1) + 
  geom_point(aes(x=factor(sequences), y=m), size=2) + 
  geom_line() + labs(x="number of sequences", y="shannon diversity") + theme_bw() + 
  facet_wrap(~id, nrow = 7, ncol = 6, scales = "free_x") + 
  theme(axis.text.x = element_text(size = 5.5, angle = 90))

quartz(height=8, width=14)
ggplot(data=pdata2[1444:2886, ], aes(x=factor(sequences), y=m)) +  
  geom_errorbar(aes(x=factor(sequences), ymin=LL, ymax=UL), col="black", width=0.1) + 
  geom_point(aes(x=factor(sequences), y=m), size=2) + 
  geom_line() + labs(x="number of sequences", y="shannon diversity") + theme_bw() +
  facet_wrap(~id, nrow = 7, ncol = 6, scales = "free_y") + 
  theme(axis.text.x = element_text(size = 5.5, angle = 90))

##
# make a dataframe with the ids and the means of all the diversity measures rarefied to 1000
# reads for those with over 1000, and up to max for those under 1000

# this just looks at the list of mean shannon for all samples with at least 1000 reads:
pdata2$m[pdata2$sequences == 1000]

# create a dataframe with the sample ids and the means at 1000 reads sampling depth
alph.divs <- data.frame(id = levels(pdata2$id), sh.mean = pdata2$m[pdata2$sequences == 1000])
# identify the samples which are missing mean shannon for 1000 reads (because fewer total):
alph.divs$id[which(is.na(alph.divs$sh.mean) == TRUE)]

# could automate this, but done by hand instead.
# identify the row numbers of missing means:
alph.divs[which(is.na(alph.divs$sh.mean) == TRUE)]

# replace them with the numbers from the highest number of reads for each sample
# print the data file for the ids which are in the list of those missing and make note of the row numbers
pdata2[which(pdata2$id %in% alph.divs$id[which(is.na(alph.divs$sh.mean) == TRUE)]), ]

alph.divs$sh.mean[which(is.na(alph.divs$sh.mean) == TRUE)] <- c(pdata2$m[c(1,213,289,394,466,904,
                                                                    989,1948,2367,2518,2679)])
alph.divs

#################
#chao1
#################
alph.c<-read.delim("~/Research/WHRI/Catalyst/Catalyst_tutorial_Sask/qiime_runthrough/nn_alpha_collated/chao1.txt",
                   na.strings = 'n/a', header = TRUE)
names(alph.c)

m<-c()
for(i in 4:81){#column numbers with samples
  mm<-tapply(alph.c[,i], factor(alph.c$sequences.per.sample), mean, na.rm=T)
  m<-c(m, mm)
}

ses<-c()
for(i in 4:81){
  ss<-tapply(alph.c[,i], factor(alph.c$sequences.per.sample), se, na.rm=T)
  ses<-c(ses, ss)
}

pdata<-data.frame(id=rep(names(alph.c[,c(4:81)]), each=37), sequences=rep(seq(100,1000,25)))
pdata2<-cbind(pdata, m, ses)

pdata2 <- within(pdata2, {
  LL <- m - 1.96 * ses
  UL <- m + 1.96 * ses
})

head(pdata2)

#make 2 figures 
quartz(height=8, width=14)
ggplot(data=pdata2[1:1443, ], aes(x=factor(sequences), y=m)) +  
  geom_errorbar(aes(x=factor(sequences), ymin=LL, ymax=UL), col="black", width=0.1) + 
  geom_point(aes(x=factor(sequences), y=m), size=2) + 
  geom_line() + labs(x="number of sequences", y="Chao1 richness") + theme_bw() + 
  facet_wrap(~id, nrow = 7, ncol = 6, scales = "free_y") + 
  theme(axis.text.x = element_text(size = 5.5, angle = 90))

quartz(height=8, width=14)
ggplot(data=pdata2[1444:2886, ], aes(x=factor(sequences), y=m)) +  
  geom_errorbar(aes(x=factor(sequences), ymin=LL, ymax=UL), col="black", width=0.1) + 
  geom_point(aes(x=factor(sequences), y=m), size=2) + 
  geom_line() + labs(x="number of sequences", y="Chao1 richness") + theme_bw() + 
  facet_wrap(~id, nrow = 7, ncol = 6, scales = "free_y") + 
  theme(axis.text.x = element_text(size = 5.5, angle = 90))


##
# make a dataframe with the ids and the means of all the diversity measures rarefied to 1000
# reads for those with over 1000, and up to max for those under 1000
# add chao1 to the existing dataframe with shannon:

alph.divs <- data.frame(alph.divs, chao1.mean = pdata2$m[pdata2$sequences == 1000])

alph.divs$chao1.mean[which(is.na(alph.divs$chao1.mean) == TRUE)] <- c(pdata2$m[c(1,213,289,394,466,904,
                                                                    989,1948,2367,2518,2679)])
alph.divs

# write a .csv file with the sample ids and the means of the alpha diversity measures:
write.csv(alph.divs, "~/Research/WHRI/Catalyst/Catalyst_tutorial_Sask/R stuff/mean divs from qiime nn.csv",
          row.names = FALSE)