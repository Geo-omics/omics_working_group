
setwd('~/Desktop/del15k/MAF5%/GWAS/')
gwas.result=read.csv('size.PC2/GAPIT.Blink.sclerotia.size.GWAS.Results.csv',header=T)

bon.threshold=-log10(0.05/nrow(gwas.result))
last.sig.p=	gwas.result[which(gwas.result$FDR_Adjusted_P.values<0.05),4] #The last significant SNP in result
FDR.threshold=-log10(last.sig.p[length(last.sig.p)]) 
gwasResults=data.frame(SNP=gwas.result$SNP,CHR=gwas.result$Chromosome,BP=gwas.result$Position,P=gwas.result$P.value)
library(gaston)
pdf('gapit.blink.sc.size.PC2.qqplot.pdf',width = 6,height = 6)
qqplot.pvalues(gwasResults$P,col.abline = "orange3",CB=T,col.CB = "lightgrey",pch=16,col='indianred4',
               main='BLINK.sc.size',xlab=expression(Expected -log[10](italic(P))), ylab=expression(Observed -log[10](italic(P))))
dev.off()
library(dplyr)
don <- gwasResults %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
#prepare the X axis
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# manhattan
library(ggplot2)
pdf('gapit.blink.sc.size.PC2.mahatton.pdf',width = 12,height = 4)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # threshold
  geom_hline(yintercept = bon.threshold,col='orange3',lwd=0.5,lty=2)+
  #geom_hline(yintercept = FDR.threshold,col='grey70',lwd=0.5,lty=2)+
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("indianred4","pink"), 47 )) +
  #scale_color_manual(values = rep(c("grey10","grey60"), 47 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ,na.value = 'NA') +
  scale_y_continuous(expand = c(0.01, 0.1) ) +     # remove space between plot area and x axis
  
  # highlight:
  geom_point(data=subset(don, P<=last.sig.p[length(last.sig.p)]), color="salmon3", size=3) +
  
  # Descriptions:
  labs(x='Contigs',y=expression(-log[10](italic(p))),caption='GWAS model: BLINK',title = 'Sclerotia size')+
  # Custom the theme:
  theme_classic() +
  theme( 
    plot.background = element_rect(),
    legend.position="none",
    plot.title= element_text(size=24,color='black',hjust=0.5),
    axis.title = element_text(size=18,color='black'),
    axis.text=element_text(size=14,color='black'),
    axis.title.x = element_text(vjust=-3),
    axis.text.x =element_text(vjust=-1),
    plot.caption=element_text(size=10,face='plain',color='#606F7B',margin=margin(t=15)),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
dev.off()
