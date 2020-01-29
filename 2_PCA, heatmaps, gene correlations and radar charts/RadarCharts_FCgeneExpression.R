
###############################################################################################################
############################################### Spider plots ##################################################
###############################################################################################################

# this is the code to perform the radar plots of the manuscript

rm(list=ls())

library(reshape)
library(ggplot2)
library(fmsb)


###############################################################################
##################################### CPS #####################################
###############################################################################

# load differential gene expression PfRBC vs uRBC data obatined with the code "CPS_analysis3.R"

# pre-immunization
i7 <- read.csv2("data/genelist3_new.csv", dec=".")
# post-immunization
c1 <- read.csv2("data/genelist4_new.csv", dec=".")
# protected post-immunization
c1p<- read.csv2("data/genelist7_new.csv", dec=".")
# unprotected pre-immunization
c1u<- read.csv2("data/genelist8_new.csv", dec=".")

i7 <- i7[,c("logFC", "Gene.Symbol")]
c1 <- c1[,c("logFC", "Gene.Symbol")]
c1p<- c1p[,c("logFC", "Gene.Symbol")]
c1u<-c1u[,c("logFC", "Gene.Symbol")]



####################   Enriched in monocytes M11.0 iRBC ####################   

# Load genes of BTM in the leading edge in the GSEA analysis of CPS protection (PfRBC/uRBC)
edge <- read.csv2("data/Leading edge btm m11_0 fc genes protection C1.csv")

monos1 <- merge(edge,i7, by="Gene.Symbol", all.x=TRUE)
monos1 <- monos1[!duplicated(monos1$Gene.Symbol), ]
monos2 <- merge(monos1,c1, by="Gene.Symbol", all.x=TRUE)
monos2 <- monos2[!duplicated(monos2$Gene.Symbol), ]
monos3 <- merge(monos2,c1p, by="Gene.Symbol", all.x=TRUE)
monos3 <- monos3[!duplicated(monos3$Gene.Symbol), ]
monos4 <- merge(monos3,c1u, by="Gene.Symbol", all.x=TRUE)
monos4 <- monos4[!duplicated(monos4$Gene.Symbol), ]

colnames(monos4) <- c("Gene.Symbol", "logFC_i7", "logFC_c1", "logFC_c1p", "logFC_c1u")

# remove genes with NA in FC
mo<- monos4[complete.cases(monos4), ]
mo <- mo[order(mo$logFC_c1p),]

# add "max" and "min" variables in 2nd and 3rd row (these are values for the plots)
mo$max <- 1
mo$min <- -3
mo <- mo[,c(1,6,7,2,4,5)]
m <- mo[,-1]

# transpose
m <- as.data.frame(t(m))

colnames(m) <- as.character(unlist(mo[,1]))

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.11,0.67,0.85,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.11,0.67,0.85,0.4) )


pdf("radar chart CPS enriched in monocytes M11_0.pdf", width = 12, height = 10)

radarchart( m  , axistype=1 , seg=4, 
            #5 segments
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, caxislabels=c("-3", "-2", "-1", "0", "1"),
            #custom labels
            vlcex=0.5, title="Enriched in monocytes (M11.0)")
legend(x=1.0, y=1.2, legend = rownames(m[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=0.9, pt.cex=3)

dev.off()




############### TLR and inflammatory signaling (M16) ##################

edge <- read.csv2("data/Leading edge btm m16 fc genes protection C1.csv", dec=".")

tlr1 <- merge(edge,i7, by="Gene.Symbol", all.x=TRUE)
tlr1 <- tlr1[!duplicated(tlr1$Gene.Symbol), ]
tlr2 <- merge(tlr1,c1, by="Gene.Symbol", all.x=TRUE)
tlr2 <- tlr2[!duplicated(tlr2$Gene.Symbol), ]
tlr3 <- merge(tlr2,c1p, by="Gene.Symbol", all.x=TRUE)
tlr3 <- tlr3[!duplicated(tlr3$Gene.Symbol), ]
tlr4 <- merge(tlr3,c1u, by="Gene.Symbol", all.x=TRUE)
tlr4 <- tlr4[!duplicated(tlr4$Gene.Symbol), ]

colnames(tlr4) <- c("Gene.Symbol", "logFC_i7", "logFC_c1", "logFC_c1p", "logFC_c1u")

# remove genes with NA in FC
t<- tlr4[complete.cases(tlr4), ]

t <- t[order(t$logFC_c1p),]

# add "max" and "min" variables in 2nd and 3rd row (these are values for the plots)
t$max <- 2
t$min <- -2
t <- t[,c(1,6,7,2,4,5)]
tm <- t[,-1]

# transpose
tm <- as.data.frame(t(tm))

colnames(tm) <- as.character(unlist(t[,1]))

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.11,0.67,0.85,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.11,0.67,0.85,0.4) )


pdf("radar chart CPS enriched in TLR and Inflammatory signaling M16.pdf", width = 12, height = 10)

radarchart( tm  , axistype=1 , seg=4, 
            #5 segments
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, caxislabels=c( "-2", "-1", "0", "1", "2"),
            #custom labels
            vlcex=0.7, title="TLR and Inflammatory signaling (M16)")
legend(x=1.0, y=1.2, legend = rownames(m[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=0.9, pt.cex=3)

dev.off()





############### Antigen presentation and immune response (M5) #############


edge <- read.csv2("data/Leading edge btm m4_3 M5 S10 fc genes protection C1.csv", dec=".")
colnames(edge)[2] <- "Gene.Symbol"

a<- merge(edge,i7, by="Gene.Symbol", all.x=TRUE)
a <- a[!duplicated(a$Gene.Symbol), ]
b <- merge(a,c1, by="Gene.Symbol", all.x=TRUE)
b <- b[!duplicated(b$Gene.Symbol), ]
c <- merge(b,c1p, by="Gene.Symbol", all.x=TRUE)
c <- c[!duplicated(c$Gene.Symbol), ]
d <- merge(c,c1u, by="Gene.Symbol", all.x=TRUE)
d <- d[!duplicated(d$Gene.Symbol), ]

d <- d[,c(1,5:8)]
colnames(d) <- c("Gene.Symbol", "logFC_i7", "logFC_c1", "logFC_c1p", "logFC_c1u")

# remove genes with NA in FC
t<- d[complete.cases(d), ]

t <- t[order(t$logFC_c1p),]

# add "max" and "min" variables in 2nd and 3rd row (these are values for the plots)
t$max <- 1
t$min <- -1
t <- t[,c(1,6,7,2,4,5)]
tm <- t[,-1]

# transpose
tm <- as.data.frame(t(tm))

colnames(tm) <- as.character(unlist(t[,1]))

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.11,0.67,0.85,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.11,0.67,0.85,0.4) )



pdf("radar chart Regulation of antigen presentation and immune response M5.pdf", width = 12, height = 10)

radarchart( tm  , axistype=1 , seg=4, 
            #5 segments
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, caxislabels=c( "-1", "-0.5", "0", "0.5", "1"),
            #custom labels
            vlcex=0.7, title="Regulation antigen presentation & immune response (M5)")
legend(x=1.0, y=1.2, legend = rownames(m[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=0.9, pt.cex=3)

dev.off()





############### Myeloid cell enriched receptors & transporters (M4.3) #############

edge <- read.csv2("data/Leading edge btm m4_3 M5 S10 fc genes protection C1.csv", dec=".")
colnames(edge)[1] <- "Gene.Symbol"

a<- merge(edge,i7, by="Gene.Symbol", all.x=TRUE)
a <- a[!duplicated(a$Gene.Symbol), ]
b <- merge(a,c1, by="Gene.Symbol", all.x=TRUE)
b <- b[!duplicated(b$Gene.Symbol), ]
c <- merge(b,c1p, by="Gene.Symbol", all.x=TRUE)
c <- c[!duplicated(c$Gene.Symbol), ]
d <- merge(c,c1u, by="Gene.Symbol", all.x=TRUE)
d <- d[!duplicated(d$Gene.Symbol), ]

d <- d[,c(1,5:8)]
colnames(d) <- c("Gene.Symbol", "logFC_i7", "logFC_c1", "logFC_c1p", "logFC_c1u")

# remove genes with NA in FC
t<- d[complete.cases(d), ]

t <- t[order(t$logFC_c1p),]

# add "max" and "min" variables in 2nd and 3rd row (these are values for the plots)
t$max <- 0
t$min <- -2
t <- t[,c(1,6,7,2,4,5)]
tm <- t[,-1]

# transpose
tm <- as.data.frame(t(tm))

colnames(tm) <- as.character(unlist(t[,1]))

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.11,0.67,0.85,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.11,0.67,0.85,0.4) )



pdf("radar chart Myeloid cell enriched receptors & transporters M4.3.pdf", width = 12, height = 10)

radarchart( tm  , axistype=1 , seg=4, 
            #5 segments
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, caxislabels=c( "-2", "-1.5", "-1", "-0.5", "0"),
            #custom labels
            vlcex=0.7, title="Myeloid cell enriched receptors & transporters (M4.3)")
legend(x=1.0, y=1.2, legend = rownames(m[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=0.9, pt.cex=3)

dev.off()




###############################################################################
##################################### RTS,S #####################################
###############################################################################

# load differential gene expression CSP7 vs DMSO data obatined with the code "RTSS_analysis2.R"

# comparators
comp <- read.csv2("data/genelist4_March16.csv", dec=".")
# RTS,S protected
rp <- read.csv2("data/genelist1_March16.csv", dec=".")
# RTS,S unprotected
ru<- read.csv2("data/genelist2_March16.csv", dec=".")


comp <- comp[,c("logFC", "Gene.Symbol")]
rp <- rp[,c("logFC", "Gene.Symbol")]
ru<- ru[,c("logFC", "Gene.Symbol")]




####################### IFN signature (M67, M68, M75, M111.1, M127, M150) ##############################

edge <- read.csv2("data/Leading edge btm IFN and monos fc genes CSP protection RTSS.csv", dec=".")
colnames(edge)[1] <- "Gene.Symbol" #select only genes from IFN signature
edge <- edge[!duplicated(edge$Gene.Symbol), ]


ifn1 <- merge(edge,comp, by="Gene.Symbol", all.x=TRUE)
ifn1 <- ifn1[!duplicated(ifn1$Gene.Symbol), ]
ifn2 <- merge(ifn1,rp, by="Gene.Symbol", all.x=TRUE)
ifn2 <- ifn2[!duplicated(ifn2$Gene.Symbol), ]
ifn3 <- merge(ifn2,ru, by="Gene.Symbol", all.x=TRUE)
ifn3 <- ifn3[!duplicated(ifn3$Gene.Symbol), ]



d <- ifn3[,-2]
colnames(d) <- c("Gene.Symbol", "logFC_comp", "logFC_rp", "logFC_ru")

# remove genes with NA in FC
t<- d[complete.cases(d), ]

t <- t[order(t$logFC_rp),]

# add "max" and "min" variables in 2nd and 3rd row (these are values for the plots)
t$max <- 0.5
t$min <- -0.5
t <- t[,c(1,5,6,2,3,4)]
tm <- t[,-1]

# transpose
tm <- as.data.frame(t(tm))

colnames(tm) <- as.character(unlist(t[,1]))

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.11,0.67,0.85,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.11,0.67,0.85,0.4) )



pdf("radar chart IFN signatures protection CSP RTSS 2.pdf", width = 12, height = 10)

radarchart( tm  , axistype=1 , seg=4, 
            #5 segments
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, caxislabels=c( "-0.5", "-0.25", "0", "0.25", "0.5"),
            #custom labels
            vlcex=0.7, title="IFN signatures (M67, M75, M111.1, M127 and M150)")
legend(x=1.0, y=1.2, legend = rownames(tm[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=0.9, pt.cex=3)

dev.off()





####################### Enriched in monocytes M11.0 ##############################


edge <- read.csv2("data/Leading edge btm IFN and monos fc genes CSP protection RTSS.csv", dec=".")
colnames(edge)[2] <- "Gene.Symbol" #select only genes from M11.0 BTM
edge <- edge[!duplicated(edge$Gene.Symbol), ]

mo1 <- merge(edge,comp, by="Gene.Symbol", all.x=TRUE)
mo1 <- mo1[!duplicated(ifn1$Gene.Symbol), ]
mo2 <- merge(mo1,rp, by="Gene.Symbol", all.x=TRUE)
mo2 <- mo2[!duplicated(ifn2$Gene.Symbol), ]
mo3 <- merge(mo2,ru, by="Gene.Symbol", all.x=TRUE)
mo3 <- mo3[!duplicated(mo3$Gene.Symbol), ]

d <- mo3[,-2]
colnames(d) <- c("Gene.Symbol", "logFC_comp", "logFC_rp", "logFC_ru")

# remove genes with NA in FC
t<- d[complete.cases(d), ]

t <- t[order(t$logFC_rp),]

# add "max" and "min" variables in 2nd and 3rd row (these are values for the plots)
t$max <- 0.5
t$min <- -0.5
t <- t[,c(1,5,6,2,3,4)]
tm <- t[,-1]

# transpose
tm <- as.data.frame(t(tm))

colnames(tm) <- as.character(unlist(t[,1]))

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.11,0.67,0.85,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.11,0.67,0.85,0.4) )



pdf("radar chart BTM monos M11_0 protection CSP RTSS.pdf", width = 12, height = 10)

radarchart( tm  , axistype=1 , seg=4, 
            #5 segments
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, caxislabels=c( "-0.5", "-0.25", "0", "0.25", "0.5"),
            #custom labels
            vlcex=0.7, title="Enriched in monocytes (II) (M11.0)")
legend(x=1.0, y=1.2, legend = rownames(tm[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=0.9, pt.cex=3)

dev.off()

