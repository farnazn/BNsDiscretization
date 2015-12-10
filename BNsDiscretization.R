#################################################################
# Load Packages
#################################################################
library(bnlearn)
library(moments)
library(ROCR)
library(caTools)
library(pROC)
library(xtable)
library(ggplot2)
library(lattice)
library(GGally)
#################################################################
# Set the Working Directory
#################################################################
# Mac
base <- "~/Dropbox/Research/Discretization/Data"
setwd(base)
# dataDIR <- paste(base, "Data", sep="/")
# Desktop
# setwd("~/Discretization/Data")
# getwd()
#################################################################
# Remove Everything
#################################################################
rm(list = ls(all = TRUE)) 
#################################################################
# Finnish Lake Data-set - Data
#################################################################
Lake.Data <- read.csv("summerAll.csv", header=TRUE, sep=",")
View(Lake.Data)
dim(Lake.Data)
colnames(Lake.Data)
# Change column names
colnames(Lake.Data) <- c('P', 'Chla', 'Type', 'Lake', 'Year', 'N', 'Month', 'Depth', 'Surface Area', 'Color')
#################################################################
# log-transform the data
Lake.Data[,'P'] = log(Lake.Data[,'P'])
Lake.Data[,'Chla'] = log(Lake.Data[,'Chla'])
Lake.Data[,'N'] = log(Lake.Data[,'N'])
# Remove one data point with log(Chla) = -23
Lake.Data <- Lake.Data[!Lake.Data[,'Chla']==min(Lake.Data[,'Chla']),]
# Remove one data point with log(N) > 9
Lake.Data <- Lake.Data[!Lake.Data[,'N']>9,]
Lake.Data <- Lake.Data[!Lake.Data[,'N']<4,]
#################################################################
pdf("Plot-Shape.pdf", width=8, height=10)
# par(mfrow=c(2,2))
# Break Points Three
Range.Chla <- range(Lake.Data[,'Chla'])[2]-range(Lake.Data[,'Chla'])[1]
Breaks.I.Chla = c(min(Lake.Data[,'Chla']), min(Lake.Data[,'Chla'])+Range.Chla/3,  max(Lake.Data[,'Chla'])-Range.Chla/3, max(Lake.Data[,'Chla']))
hist(Lake.Data[,'Chla'], breaks = Breaks.I.Chla, col='black', main='Three Intervals', xlab='log(Chlorophyll a)',border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
# Break Points Four
Range.Chla <- range(Lake.Data[,'Chla'])[2]-range(Lake.Data[,'Chla'])[1]
Breaks.I.Chla = c(min(Lake.Data[,'Chla']), min(Lake.Data[,'Chla'])+Range.Chla/4,  min(Lake.Data[,'Chla'])+2*Range.Chla/4, min(Lake.Data[,'Chla'])+3*Range.Chla/4, max(Lake.Data[,'Chla']))
hist(Lake.Data[,'Chla'], breaks = Breaks.I.Chla, col='black', main='Four Intervals', xlab='log(Chlorophyll a)',border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
# Break Points Five
Range.Chla <- range(Lake.Data[,'Chla'])[2]-range(Lake.Data[,'Chla'])[1]
Breaks.I.Chla = c(min(Lake.Data[,'Chla']), min(Lake.Data[,'Chla'])+Range.Chla/5,  min(Lake.Data[,'Chla'])+2*Range.Chla/5, min(Lake.Data[,'Chla'])+3*Range.Chla/5, min(Lake.Data[,'Chla'])+4*Range.Chla/5, max(Lake.Data[,'Chla']))
hist(Lake.Data[,'Chla'], breaks = Breaks.I.Chla, col='black', main='Five Intervals', xlab='log(Chlorophyll a)',border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
invisible(dev.off())
#################################################################
# Pairs Plot
Data <- Lake.Data
colnames(Data) = c("Phosphorus", "Chlorophylla", "Type", "Lake", "Year", "Nitrogen", "Month", "Depth", "Surface Area", "Color")

pdf("PairsPlot.pdf", width=12, height=10)
ggpairs(Data[, c("Phosphorus","Nitrogen","Chlorophylla")],
        lower=list(continuous="smooth", params=c(colour="darkgray")),
        diag=list(continuous="density", params=c(colour="black")), 
        upper=list(params=c(Size=20, colour="black")),
        axisLabels="show") + theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "dashed", colour = "black", fill = NA),
        axis.title = element_text(size = 200),
        axis.text.x = element_text(size=20, colour = "black"),
        axis.text.y = element_text(size=20, colour = "black"))
invisible(dev.off())
#################################################################
# Sample 10%/90% of eash lake type
# We will use 90% of data for model development & 10% of data for model evaluation
for (i in 1:9){
assign(paste("Lake.",i,sep=""), Lake.Data[which(Lake.Data$Type==i),])  
}
Type.Size <- c(nrow(Lake.1), nrow(Lake.2), nrow(Lake.3), nrow(Lake.4)
             , nrow(Lake.5), nrow(Lake.6), nrow(Lake.7), nrow(Lake.8)
             , nrow(Lake.9))

for (j in 1:10){
  set.seed(j+10)
  assign(paste("Model.",j,sep=""), 0)
  assign(paste("Evaluation.",j,sep=""), 0)

  for(i in 1:9){
      assign(paste("L.",i,".M.",j,sep=""), get(noquote(paste('Lake.', i, sep="")))[sample(1:Type.Size[i], size= round(0.9*Type.Size[i]), replace=FALSE),]
    )
      assign(paste("L.", i,".M.E.",j, sep=""), get(noquote(paste('Lake.', i, sep="")))[-sample(1:Type.Size[i], size= round(0.9*Type.Size[i]), replace=FALSE),]
    )
    # Keeping variables: N, P, and Chla
      Vars <- c('P', 'Chla', 'N')
      assign(paste("L.",i,".M.",j,sep=""), get(noquote(paste('L.',i,".M.",j, sep="")))[,Vars])
    # Keeping variables: N, P, and Chla
      Vars <- c('P', 'Chla', 'N')
      assign(paste("L.",i,".M.",j,sep=""), get(noquote(paste('L.',i,".M.",j, sep="")))[,Vars])
      assign(paste("L.",i,".M.E.",j,sep=""), get(noquote(paste('L.',i,".M.E.",j, sep="")))[,Vars])
      
      assign(paste("Model.",j,sep=""), rbind(get(noquote(paste("Model.",j,sep=""))), get(noquote(paste("L.",i,".M.",j,sep="")))))
      assign(paste("Evaluation.",j,sep=""), rbind(get(noquote(paste("Evaluation.",j,sep=""))), get(noquote(paste("L.",i,".M.E.",j,sep="")))))
  }
}
#################################################################
#BN discretized with equal interval
#################################################################
# Break Points
Range.Chla <- range(Lake.Data[,'Chla'])[2]-range(Lake.Data[,'Chla'])[1]
Breaks.I.Chla = c(min(Lake.Data[,'Chla']), min(Lake.Data[,'Chla'])+Range.Chla/3,  max(Lake.Data[,'Chla'])-Range.Chla/3, max(Lake.Data[,'Chla']))

Range.P <- range(Lake.Data[,'P'])[2]-range(Lake.Data[,'P'])[1]
Breaks.I.P = c(min(Lake.Data[,'P']), min(Lake.Data[,'P'])+Range.P/3,  max(Lake.Data[,'P'])-Range.P/3, max(Lake.Data[,'P']))

Range.N <- range(Lake.Data[,'N'])[2]-range(Lake.Data[,'N'])[1]
Breaks.I.N = c(min(Lake.Data[,'N']), min(Lake.Data[,'N'])+Range.N/3,  max(Lake.Data[,'N'])-Range.N/3, max(Lake.Data[,'N']))

for (j in 1:10){
  assign(paste("Interval.",j,sep=""), data.frame(P = cut(get(noquote(paste("Model.",j,sep="")))$P, breaks = Breaks.I.P, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                               , Chla = cut(get(noquote(paste("Model.",j,sep="")))$Chla, breaks = Breaks.I.Chla, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                               , N = cut(get(noquote(paste("Model.",j,sep="")))$N, breaks = Breaks.I.N, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)))
  assign(paste("Interval.",j,sep=""), get(noquote(paste("Interval.",j,sep="")))[complete.cases(get(noquote(paste("Interval.",j,sep="")))),])
  
  levels(Interval.1$N) <- c("Low Nitrogen","Medium Nitrogen","High Nitrogen")
  levels(Interval.1$P) <- c("Low Phosphorus","Medium Phosphorus","High Phosphorus")
  levels(Interval.1$Chla) <- c("Low","Medium","High")

  assign(paste("Interval.E.",j,sep=""), data.frame(P = cut(get(noquote(paste("Evaluation.",j,sep="")))$P, breaks = Breaks.I.P, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                               , Chla = cut(get(noquote(paste("Evaluation.",j,sep="")))$Chla, breaks = Breaks.I.Chla, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                               , N = cut(get(noquote(paste("Evaluation.",j,sep="")))$N, breaks = Breaks.I.N, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)))
  assign(paste("Interval.E.",j,sep=""), get(noquote(paste("Interval.E.",j,sep="")))[complete.cases(get(noquote(paste("Interval.E.",j,sep="")))),])

  levels(Interval.E.1$N) <- c("Low Nitrogen","Medium Nitrogen","High Nitrogen")
  levels(Interval.E.1$P) <- c("Low Phosphorus","Medium Phosphorus","High Phosphorus")
  levels(Interval.E.1$Chla) <- c("Low","Medium","High")

  Model.Interval <- empty.graph(names(get(noquote(paste("Interval.",j,sep="")))))
  modelstring(Model.Interval) <- "[P],[N],[Chla|N:P]"
  assign(paste("fit.Interval.",j,sep=""), bn.fit(Model.Interval, get(noquote(paste("Interval.",j,sep="")))))

#   Model evaluation for equal interval
  assign(paste("logLik.Interval.",j,sep=""), logLik(get(noquote(paste("fit.Interval.",j,sep=""))), get(noquote(paste("Interval.",j,sep="")))))
#   AIC = -2 ( ln ( likelihood )) + 2 K
  assign(paste("AIC.Interval.",j,sep=""), AIC(get(noquote(paste("fit.Interval.",j,sep=""))), get(noquote(paste("Interval.",j,sep=""))), k = 2))
  assign(paste("BIC.Interval.",j,sep=""), BIC(get(noquote(paste("fit.Interval.",j,sep=""))), get(noquote(paste("Interval.",j,sep="")))))
#   predict class on training sample
  assign(paste("Interval.Pred.ChlaT.",j,sep=""), predict(get(noquote(paste("fit.Interval.",j,sep=""))), "Chla", get(noquote(paste("Interval.",j,sep="")))))
#   AUC for evaluation sample  
  assign(paste("Interval.Pred.ChlaE.",j,sep=""), predict(get(noquote(paste("fit.Interval.",j,sep=""))), "Chla", get(noquote(paste("Interval.E.",j,sep="")))))
#   make confusion matrix for training sample
  assign(paste("Interval.CM.ChlaE.",j,sep=""), table(get(noquote(paste("Interval.E.",j,sep="")))[,'Chla'],get(noquote(paste("Interval.Pred.ChlaE.",j,sep="")))))
#   compute accuracy for training sample
  assign(paste("Interval.Accuracy.",j,sep=""), sum(diag(table(get(noquote(paste("Interval.E.",j,sep="")))[,'Chla'],get(noquote(paste("Interval.Pred.ChlaE.",j,sep="")))))/length(get(noquote(paste("Interval.Pred.ChlaE.",j,sep=""))))))

  assign(paste("Interval.Pred.ChlaE.",j,sep=""), data.frame(get(noquote(paste("Interval.Pred.ChlaE.",j,sep="")))))
  assign(paste("AUC.Interval.",j,sep=""), auc(as.matrix(get(noquote(paste("Interval.Pred.ChlaE.",j,sep="")))), get(noquote(paste("Interval.E.",j,sep="")))[ ,'Chla']))
}
#################################################################
# Sum of Square Errors
#################################################################
SSE= NULL
for (i in 1:1943){
  if(Interval.Pred.ChlaE.1[i,1]=="(3.51,6.41]") {SSE[i]=4.96}
  if(Interval.Pred.ChlaE.1[i,1]=="(0.603,3.51]") {SSE[i]=2.0565}
  if(Interval.Pred.ChlaE.1[i,1]=="(-2.3,0.603]") {SSE[i]=-0.8485}
}

Square= NULL
 for (i in 1:1943){Square[i]=(SSE[i]-Evaluation.1[i+1,"Chla"])^2
                   }

sum(Square, na.rm=TRUE)
#################################################################
# LATEX tables and figures
# Calculating logLik, AIC, BIC, Accuracy, and AUC
#################################################################
# LATEX OUTPUT by xtable 1.7-3 package
Interval <- matrix(NA, 0, 5)

for(j in 1:10){
  assign(paste("Interval",sep="")
         , rbind(get(noquote(paste("Interval",sep="")))
         , cbind(get(noquote(paste("logLik.Interval.",j,sep="")))
               , get(noquote(paste("AIC.Interval.",j,sep="")))
               , get(noquote(paste("BIC.Interval.",j,sep="")))
               , get(noquote(paste("Interval.Accuracy.",j,sep="")))
               , get(noquote(paste("AUC.Interval.",j,sep=""))))))
}

colnames(Interval) = c("logLik", "AIC", "BIC","Accuracy", "AUC")

Interval <- rbind(Interval, 
      cbind(mean(Interval[,"logLik"]), mean(Interval[,"AIC"]), mean(Interval[,"BIC"]), mean(Interval[,"Accuracy"]), mean(Interval[,"AUC"])))


Interval.Table <- xtable(Interval, digits=3)
print(Interval.Table)
#################################################################
# Convert CPTs to LATEX output
#################################################################
# CPTs to LATEX N
print(xtable(fit.Interval.1$N$prob), floating=FALSE)
# CPTs to LATEX P
print(xtable(fit.Interval.1$P$prob), floating=FALSE)
# CPTs to LATEX Chla
print(xtable(fit.Interval.1$Chla$prob[1:3, 1:3, 1]), floating=FALSE)
print(xtable(fit.Interval.1$Chla$prob[1:3, 1:3, 2]), floating=FALSE)
print(xtable(fit.Interval.1$Chla$prob[1:3, 1:3, 3]), floating=FALSE)
#################################################################
# DAG model
#################################################################
pdf("Plot-Interval.pdf", width=8, height=10)
plot(Model.Interval)
#################################################################
# Conditional Probability Tables
#################################################################
bn.fit.barchart(fit.Interval.1$P, ylab="Phosphorus Levels")
bn.fit.barchart(fit.Interval.1$N, ylab="Nitrogen Levels")
bn.fit.barchart(fit.Interval.1$Chla, ylab="Chlorophyll a Levels")
#################################################################
# Histogram of data for model building: continuous vs discretized
hist(Model.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
hist(Model.1[ ,"Chla"], breaks=c(-2.31, 0.6027629, 3.5081109, 6.42), main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Interval.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Interval", 1,line=2,at=c(-3),col="darkgreen")
#################################################################
# Histogram of data for model evaluation: continuous vs discretized
hist(Evaluation.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
hist(Evaluation.1[ ,"Chla"], breaks=c(-2.31, 0.6027629, 3.5081109, 6.42), main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Interval.E.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Interval", 1,line=2,at=c(-3),col="darkgreen")
#################################################################
# Histogram of continuous data for model evaluation vs discretized model prediction
hist(Evaluation.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
pred.1 <- runif(length(Interval.Pred.ChlaE.1[Interval.Pred.ChlaE.1=="[-2.3,0.603]"]), min = -2.3, max = 0.6)
pred.2 <- runif(length(Interval.Pred.ChlaE.1[Interval.Pred.ChlaE.1=="(0.603,3.51]"]), min = 0.604, max = 3.51)
pred.3 <- runif(length(Interval.Pred.ChlaE.1[Interval.Pred.ChlaE.1=="(3.51,6.41]"]), min = 3.52, max = 6.41)
pred <- c(pred.1, pred.2, pred.3)
hist(pred, breaks=c(-2.31, 0.6027629, 3.5081109, 6.42), , main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Interval.E.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Interval", 1,line=2,at=c(-3),col="darkgreen")
invisible(dev.off())
#################################################################
pdf("Plot-Interval-3.pdf", width=8, height=10)
temp <- trellis.par.get()
temp$axis.text$cex=2.4
temp$plot.polygon$col <- "darkgray" 
trellis.par.set(temp)
bn.fit.barchart(fit.Interval.1$P, xlab= list(label="Probabilities", cex=3), ylab="", main="")
bn.fit.barchart(fit.Interval.1$N, xlab= list(label="Probabilities", cex=3), ylab="",main="")
temp$axis.text$cex=1.2
trellis.par.set(temp)
bn.fit.barchart(fit.Interval.1$Chla, ylab=list(label="Chlorophyll a Levels", cex=1.8), xlab= list(label="Probabilities", cex=1.8), main= list(label="Conditional Probabilities", cex=1.8))
invisible(dev.off())
#################################################################
# palette <- palette(gray(seq(0, 1,len=10))) 
# trellis.par.set(list(par.xlab.text=list(cex=0.85) 
# temp$superpose.polygon$col=list(col=palette()) 
#################################################################
# BN discretized with equal quantile
#################################################################
# Break Points
Breaks.Q.Chla=c(quantile(Lake.Data[,"Chla"], probs = seq(0, 1, by = 1/3)))
Breaks.Q.P=c(quantile(Lake.Data[,"P"], probs = seq(0, 1, by = 1/3)))
Breaks.Q.N=c(quantile(Lake.Data[,"N"], probs = seq(0, 1, by = 1/3)))

for (j in 1:10){
  assign(paste("Quantile.",j,sep=""), data.frame(P = cut(get(noquote(paste("Model.",j,sep="")))$P, breaks = Breaks.Q.P, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                                 , Chla = cut(get(noquote(paste("Model.",j,sep="")))$Chla, breaks = Breaks.Q.Chla, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                                 , N = cut(get(noquote(paste("Model.",j,sep="")))$N, breaks = Breaks.Q.N, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)))
  assign(paste("Quantile.",j,sep=""), get(noquote(paste("Quantile.",j,sep="")))[complete.cases(get(noquote(paste("Quantile.",j,sep="")))),])
  
  
  levels(Quantile.1$N) <- c("Low Nitrogen","Medium Nitrogen","High Nitrogen")
  levels(Quantile.1$P) <- c("Low Phosphorus","Medium Phosphorus","High Phosphorus")
  levels(Quantile.1$Chla) <- c("Low","Medium","High")
  assign(paste("Quantile.E.",j,sep=""), data.frame(P = cut(get(noquote(paste("Evaluation.",j,sep="")))$P, breaks = Breaks.Q.P, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                                   , Chla = cut(get(noquote(paste("Evaluation.",j,sep="")))$Chla, breaks = Breaks.Q.Chla, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                                   , N = cut(get(noquote(paste("Evaluation.",j,sep="")))$N, breaks = Breaks.Q.N, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)))
  assign(paste("Quantile.E.",j,sep=""), get(noquote(paste("Quantile.E.",j,sep="")))[complete.cases(get(noquote(paste("Quantile.E.",j,sep="")))),])

  levels(Quantile.E.1$N) <- c("Low Nitrogen","Medium Nitrogen","High Nitrogen")
  levels(Quantile.E.1$P) <- c("Low Phosphorus","Medium Phosphorus","High Phosphorus")
  levels(Quantile.E.1$Chla) <- c("Low","Medium","High")

  Model.Quantile <- empty.graph(names(get(noquote(paste("Quantile.",j,sep="")))))
  modelstring(Model.Quantile) <- "[P],[N],[Chla|N:P]"
  assign(paste("fit.Quantile.",j,sep=""), bn.fit(Model.Quantile, get(noquote(paste("Quantile.",j,sep="")))))
  
#   Model evaluation for equal quantile
  assign(paste("logLik.Quantile.",j,sep=""), logLik(get(noquote(paste("fit.Quantile.",j,sep=""))), get(noquote(paste("Quantile.",j,sep="")))))
#   AIC = -2 ( ln ( likelihood )) + 2 K
  assign(paste("AIC.Quantile.",j,sep=""), AIC(get(noquote(paste("fit.Quantile.",j,sep=""))), get(noquote(paste("Quantile.",j,sep=""))), k = 2))
  assign(paste("BIC.Quantile.",j,sep=""), BIC(get(noquote(paste("fit.Quantile.",j,sep=""))), get(noquote(paste("Quantile.",j,sep="")))))
#   AUC for evaluation sample  
  assign(paste("Quantile.Pred.ChlaE.",j,sep=""), predict(get(noquote(paste("fit.Quantile.",j,sep=""))), "Chla", get(noquote(paste("Quantile.E.",j,sep="")))))
#   make confusion matrix for training sample
  assign(paste("Quantile.CM.ChlaP.",j,sep=""), table(get(noquote(paste("Quantile.E.",j,sep="")))[,'Chla'],get(noquote(paste("Quantile.Pred.ChlaE.",j,sep="")))))
#   compute accuracy for training sample
  assign(paste("Quantile.Accuracy.",j,sep=""), sum(diag(table(get(noquote(paste("Quantile.E.",j,sep="")))[,'Chla'],get(noquote(paste("Quantile.Pred.ChlaE.",j,sep="")))))/length(get(noquote(paste("Quantile.Pred.ChlaE.",j,sep=""))))))
 
  assign(paste("Quantile.Pred.ChlaE.",j,sep=""), data.frame(get(noquote(paste("Quantile.Pred.ChlaE.",j,sep="")))))
  assign(paste("AUC.Quantile.",j,sep=""), auc(as.matrix(get(noquote(paste("Quantile.Pred.ChlaE.",j,sep="")))), get(noquote(paste("Quantile.E.",j,sep="")))[ ,'Chla']))
}
#################################################################
# CPTs
#################################################################
pdf("Plot-Quantile-3.pdf", width=8, height=10)
temp$axis.text$cex=2.4
temp$plot.polygon$col <- "darkgray" 
trellis.par.set(temp)
bn.fit.barchart(fit.Quantile.1$P, ylab= "", xlab= list(label="Probabilities", cex=3), main="")
bn.fit.barchart(fit.Quantile.1$N, ylab="", xlab= list(label="Probabilities", cex=3),main="")
temp$axis.text$cex=1.2
trellis.par.set(temp)
bn.fit.barchart(fit.Quantile.1$Chla, ylab=list(label="Chlorophyll a Levels", cex=1.8), xlab= list(label="Probabilities", cex=1.8), main= list(label="Conditional Probabilities", cex=1.8))
invisible(dev.off())
#################################################################
# Quantile.Pred.ChlaE.Trial= predict(fit.Quantile.1, "Chla", Quantile.E.1[Quantile.E.1[,3]=="[3.43,5.99]",])
# Interval.Pred.ChlaE.Trial= predict(fit.Interval.1, "Chla", Interval.E.1[Interval.E.1[,3]=="[3.43,5.11]",])
# MM.Pred.ChlaE.Trial= predict(fit.MM.1, "Chla", MM.E.1[MM.E.1[,3]=="[3.43,5.7]",])
#################################################################
# Sum of Square Errors
#################################################################
SSE.Q= NULL
for (i in 1:1943){
  if(Quantile.Pred.ChlaE.1[i,1]=="[-2.3,1.72]") {SSE.Q[i]=-0.29}
  if(Quantile.Pred.ChlaE.1[i,1]=="(1.72,2.64]") {SSE.Q[i]=2.18}
  if(Quantile.Pred.ChlaE.1[i,1]=="(2.64,6.41]") {SSE.Q[i]=4.525}
}

Square.Q= NULL
for (i in 1:1943){Square.Q[i]=(SSE.Q[i]-Evaluation.1[i+1,"Chla"])^2
}

sum(Square.Q, na.rm=TRUE)
#################################################################
# LATEX tables and figures
# Calculating logLik, AIC, BIC, Accuracy, and AUC
#################################################################
# LATEX OUTPUT by xtable 1.7-3 package
Quantile <- matrix(NA, 0, 5)

for(j in 1:10){
  assign(paste("Quantile",sep="")
         , rbind(get(noquote(paste("Quantile",sep="")))
                 , cbind(get(noquote(paste("logLik.Quantile.",j,sep="")))
                         , get(noquote(paste("AIC.Quantile.",j,sep="")))
                         , get(noquote(paste("BIC.Quantile.",j,sep="")))
                         , get(noquote(paste("Quantile.Accuracy.",j,sep="")))
                         , get(noquote(paste("AUC.Quantile.",j,sep=""))))))
 
}

colnames(Quantile) = c("logLik", "AIC", "BIC","Accuracy", "AUC")

Quantile <- rbind(Quantile, 
                  cbind(mean(Quantile[,"logLik"]), mean(Quantile[,"AIC"]), mean(Quantile[,"BIC"]), mean(Quantile[,"Accuracy"]), mean(Quantile[,"AUC"])))

Quantile.Table <- xtable(Quantile, digits=3)
print(Quantile.Table)
#################################################################
# CPTs to LATEX N
print(xtable(fit.Quantile.1$N$prob), floating=FALSE)
# CPTs to LATEX P
print(xtable(fit.Quantile.1$P$prob), floating=FALSE)
# CPTs to LATEX Chla
print(xtable(fit.Quantile.1$Chla$prob[1:3, 1:3, 1]), floating=FALSE)
print(xtable(fit.Quantile.1$Chla$prob[1:3, 1:3, 2]), floating=FALSE)
print(xtable(fit.Quantile.1$Chla$prob[1:3, 1:3, 3]), floating=FALSE)
#################################################################
pdf("Plot-Quantile.pdf", width=8, height=10)

bn.fit.barchart(fit.Quantile.1$P, ylab="Phosphorus Levels")
bn.fit.barchart(fit.Quantile.1$N, ylab="Nitrogen Levels")
bn.fit.barchart(fit.Quantile.1$Chla, ylab="Chlorophyll a Levels")

# Histogram of data for model building: continuous vs discretized
hist(Model.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
hist(Model.1[ ,"Chla"], breaks=c(-2.3025851, 1.722767, 2.639057, 6.4134590), main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Quantile.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Quantile", 1,line=2,at=c(-3),col="darkgreen")

# Histogram of data for model evaluation: continuous vs discretized
hist(Evaluation.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
hist(Evaluation.1[ ,"Chla"], breaks=c(-2.3025851, 1.722767, 2.639057, 6.4134590), main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Quantile.E.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Quantile", 1,line=2,at=c(-3),col="darkgreen")

# Histogram of continuous data for model evaluation vs discretized model prediction
hist(Evaluation.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
pred.1 <- runif(length(Quantile.Pred.ChlaE.1[Quantile.Pred.ChlaE.1=="[-2.3,1.72]"]), min = -2.3025851, max = 1.722767)
pred.2 <- runif(length(Quantile.Pred.ChlaE.1[Quantile.Pred.ChlaE.1=="(1.72,2.64]"]), min = 1.722767, max = 2.639057)
pred.3 <- runif(length(Quantile.Pred.ChlaE.1[Quantile.Pred.ChlaE.1=="(2.64,6.41]"]), min = 2.639057, max = 6.4134590)
pred <- c(pred.1, pred.2, pred.3)
hist(pred, breaks=c(-2.3025851, 1.722767, 2.639057, 6.4134590), , main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Quantile.E.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Quantile", 1,line=2,at=c(-3),col="darkgreen")
invisible(dev.off())
#################################################################
# Break Points: BN discretized with moment matching method for P
#################################################################
mu.P = mean(Lake.Data[,'P'])
sd.P = sd(Lake.Data[,'P'])
s.P = skewness(Lake.Data[,'P'], na.rm = FALSE)
k.P = kurtosis(Lake.Data[,'P'], na.rm = FALSE)
f.P = moment(Lake.Data[,'P'], order = 5, central=TRUE, na.rm = FALSE)/sd(Lake.Data[,'P'])^5
c0.P = (f.P-2*s.P*k.P+s.P^3)/(k.P-s.P^2-1)
c1.P = (s.P*f.P-k.P^2+k.P-s.P^2)/(k.P-s.P^2-1)
c2.P = (-f.P+s.P*k.P+s.P)/(k.P-s.P^2-1)

# polyroot: Find zeros of a real or complex polynomial.
root.P = polyroot(c(c0.P, c1.P, c2.P, 1))
Im(root.P)
root.P = Re(root.P)

pa.P = (1+root.P[2]*root.P[3])/((root.P[2]-root.P[1])*(root.P[3]-root.P[1]))
pb.P = (1+root.P[1]*root.P[3])/((root.P[1]-root.P[2])*(root.P[3]-root.P[2]))
pc.P = (1+root.P[1]*root.P[2])/((root.P[1]-root.P[3])*(root.P[2]-root.P[3]))
a.P = mu.P+sd.P*min(root.P)
b.P = mu.P+sd.P*median(root.P)
c.P = mu.P+sd.P*max(root.P)

c(a.P,b.P,c.P)
c(pa.P,pb.P,pc.P)
hist(Lake.Data[,'P'])

Break.MM.P.1 = min(Lake.Data[,'P'])
Break.MM.P.2 = a.P+(pa.P/(pa.P+pb.P))*(b.P-a.P)
Break.MM.P.3 = b.P+ (pb.P/(pb.P+pc.P))*(c.P-b.P)
Break.MM.P.4 = max(Lake.Data[,'P'])

Breaks.MM.P <- c(Break.MM.P.1, Break.MM.P.2, Break.MM.P.3, Break.MM.P.4)
#################################################################
# Break Points: BN discretized with moment matching method for N
#################################################################
mu.N = mean(Lake.Data[,'N'])
sd.N = sd(Lake.Data[,'N'])
s.N = skewness(Lake.Data[,'N'], na.rm = FALSE)
k.N = kurtosis(Lake.Data[,'N'], na.rm = FALSE)
f.N = moment(Lake.Data[,'N'], order = 5, central=TRUE, na.rm = FALSE)/sd(Lake.Data[,'N'])^5
c0.N = (f.N-2*s.N*k.N+s.N^3)/(k.N-s.N^2-1)
c1.N = (s.N*f.N-k.N^2+k.N-s.N^2)/(k.N-s.N^2-1)
c2.N = (-f.N+s.N*k.N+s.N)/(k.N-s.N^2-1)

# polyroot: Find zeros of a real or complex polynomial.
root.N = polyroot(c(c0.N, c1.N, c2.N, 1))
Im(root.N)
root.N = Re(root.N)

pa.N = (1+root.N[2]*root.N[3])/((root.N[2]-root.N[1])*(root.N[3]-root.N[1]))
pb.N = (1+root.N[1]*root.N[3])/((root.N[1]-root.N[2])*(root.N[3]-root.N[2]))
pc.N = (1+root.N[1]*root.N[2])/((root.N[1]-root.N[3])*(root.N[2]-root.N[3]))
a.N = mu.N+sd.N*min(root.N)
b.N = mu.N+sd.N*median(root.N)
c.N = mu.N+sd.N*max(root.N)

c(a.N,b.N,c.N)
hist(Lake.Data[,'N'])

Break.MM.N.1 = min(Lake.Data[,'N'])
Break.MM.N.2 = a.N+(pa.N/(pa.N+pb.N))*(b.N-a.N)
Break.MM.N.3 = b.N+ (pb.N/(pb.N+pc.N))*(c.N-b.N)
Break.MM.N.4 = max(Lake.Data[,'N'])

Breaks.MM.N <- c(Break.MM.N.1, Break.MM.N.2, Break.MM.N.3, Break.MM.N.4)
#################################################################
# Break Points: BN discretized with moment matching method for Chla
#################################################################
mu.Chla = mean(Lake.Data[,'Chla'])
sd.Chla = sd(Lake.Data[,'Chla'])
s.Chla = skewness(Lake.Data[,'Chla'], na.rm = FALSE)
k.Chla = kurtosis(Lake.Data[,'Chla'], na.rm = FALSE)
f.Chla = moment(Lake.Data[,'Chla'], order = 5, central=TRUE, na.rm = FALSE)/sd(Lake.Data[,'Chla'])^5
c0.Chla = (f.Chla-2*s.Chla*k.Chla+s.Chla^3)/(k.Chla-s.Chla^2-1)
c1.Chla = (s.Chla*f.Chla-k.Chla^2+k.Chla-s.Chla^2)/(k.Chla-s.Chla^2-1)
c2.Chla = (-f.Chla+s.Chla*k.Chla+s.Chla)/(k.Chla-s.Chla^2-1)

# polyroot: Find zeros of a real or complex polynomial.
root.Chla = polyroot(c(c0.Chla, c1.Chla, c2.Chla, 1))
Im(root.Chla)
root.Chla = Re(root.Chla)

pa.Chla = (1+root.Chla[2]*root.Chla[3])/((root.Chla[2]-root.Chla[1])*(root.Chla[3]-root.Chla[1]))
pb.Chla = (1+root.Chla[1]*root.Chla[3])/((root.Chla[1]-root.Chla[2])*(root.Chla[3]-root.Chla[2]))
pc.Chla = (1+root.Chla[1]*root.Chla[2])/((root.Chla[1]-root.Chla[3])*(root.Chla[2]-root.Chla[3]))
a.Chla = mu.Chla+sd.Chla*min(root.Chla)
b.Chla = mu.Chla+sd.Chla*median(root.Chla)
c.Chla = mu.Chla+sd.Chla*max(root.Chla)

c(a.Chla,b.Chla,c.Chla)

Break.MM.Chla.1 = min(Lake.Data[,'Chla'])
Break.MM.Chla.2 = a.Chla+(pa.Chla/(pa.Chla+pb.Chla))*(b.Chla-a.Chla)
Break.MM.Chla.3 = b.Chla+ (pb.Chla/(pb.Chla+pc.Chla))*(c.Chla-b.Chla)
Break.MM.Chla.4 = max(Lake.Data[,'Chla'])
Breaks.MM.Chla <- c(Break.MM.Chla.1, Break.MM.Chla.2, Break.MM.Chla.3, Break.MM.Chla.4)
#################################################################
# BN discretized with moment matching method
#################################################################
for (j in 1:10){
  assign(paste("MM.",j,sep=""), data.frame(P = cut(get(noquote(paste("Model.",j,sep="")))$P, breaks = Breaks.MM.P, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                                 , Chla = cut(get(noquote(paste("Model.",j,sep="")))$Chla, breaks = Breaks.MM.Chla, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                                 , N = cut(get(noquote(paste("Model.",j,sep="")))$N, breaks = Breaks.MM.N, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)))
  assign(paste("MM.",j,sep=""), get(noquote(paste("MM.",j,sep="")))[complete.cases(get(noquote(paste("MM.",j,sep="")))),])
  
  
  levels(MM.1$N) <- c("Low Nitrogen","Medium Nitrogen","High Nitrogen")
  levels(MM.1$P) <- c("Low Phosphorus","Medium Phosphorus","High Phosphorus")
  levels(MM.1$Chla) <- c("Low","Medium","High")
  assign(paste("MM.E.",j,sep=""), data.frame(P = cut(get(noquote(paste("Evaluation.",j,sep="")))$P, breaks = Breaks.MM.P, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                                   , Chla = cut(get(noquote(paste("Evaluation.",j,sep="")))$Chla, breaks = Breaks.MM.Chla, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)
                                                   , N = cut(get(noquote(paste("Evaluation.",j,sep="")))$N, breaks = Breaks.MM.N, ordered = TRUE, include.lowest=TRUE, include.highest=TRUE)))
  assign(paste("MM.E.",j,sep=""), get(noquote(paste("MM.E.",j,sep="")))[complete.cases(get(noquote(paste("MM.E.",j,sep="")))),])
  levels(MM.E.1$N) <- c("Low Nitrogen","Medium Nitrogen","High Nitrogen")
  levels(MM.E.1$P) <- c("Low Phosphorus","Medium Phosphorus","High Phosphorus")
  levels(MM.E.1$Chla) <- c("Low","Medium","High")
  Model.MM <- empty.graph(names(get(noquote(paste("MM.",j,sep="")))))
  modelstring(Model.MM) <- "[P],[N],[Chla|N:P]"
  assign(paste("fit.MM.",j,sep=""), bn.fit(Model.MM, get(noquote(paste("MM.",j,sep="")))))
  
  #   Model evaluation for moment matching
  assign(paste("logLik.MM.",j,sep=""), logLik(get(noquote(paste("fit.MM.",j,sep=""))), get(noquote(paste("MM.",j,sep="")))))
  #   AIC = -2 ( ln ( likelihood )) + 2 K
  assign(paste("AIC.MM.",j,sep=""), AIC(get(noquote(paste("fit.MM.",j,sep=""))), get(noquote(paste("MM.",j,sep=""))), k = 2))
  assign(paste("BIC.MM.",j,sep=""), BIC(get(noquote(paste("fit.MM.",j,sep=""))), get(noquote(paste("MM.",j,sep="")))))
  #   AUC for evaluation sample  
  assign(paste("MM.Pred.ChlaE.",j,sep=""), predict(get(noquote(paste("fit.MM.",j,sep=""))), "Chla", get(noquote(paste("MM.E.",j,sep="")))))
  #   make confusion matrix for testing sample
  assign(paste("MM.CM.ChlaT.",j,sep=""), table(get(noquote(paste("MM.E.",j,sep="")))[,'Chla'],get(noquote(paste("MM.Pred.ChlaE.",j,sep="")))))
  #   compute accuracy for training sample
  assign(paste("MM.Accuracy.",j,sep=""), sum(diag(table(get(noquote(paste("MM.E.",j,sep="")))[,'Chla'],get(noquote(paste("MM.Pred.ChlaE.",j,sep="")))))/length(get(noquote(paste("MM.Pred.ChlaE.",j,sep=""))))))
  assign(paste("MM.Pred.ChlaE.",j,sep=""), data.frame(get(noquote(paste("MM.Pred.ChlaE.",j,sep="")))))
  assign(paste("AUC.MM.",j,sep=""), auc(as.matrix(get(noquote(paste("MM.Pred.ChlaE.",j,sep="")))), get(noquote(paste("MM.E.",j,sep="")))[ ,'Chla']))
}
#################################################################
pdf("Plot-MM-3.pdf", width=8, height=10)
temp$axis.text$cex=2.4
temp$plot.polygon$col <- "darkgray" 
trellis.par.set(temp)
bn.fit.barchart(fit.MM.1$P, ylab="", xlab= list(label="Probabilities", cex=3), main="")
bn.fit.barchart(fit.MM.1$N, ylab="", xlab= list(label="Probabilities", cex=3),main="")
temp$axis.text$cex=1.2
trellis.par.set(temp)
bn.fit.barchart(fit.MM.1$Chla, ylab=list(label="Chlorophyll a Levels", cex=1.8), xlab= list(label="Probabilities", cex=1.8), main= list(label="Conditional Probabilities", cex=1.8))
invisible(dev.off())
#################################################################
# Sum of Square Errors
#################################################################
SSE.MM= NULL
for (i in 1:1943){
  if(MM.Pred.ChlaE.1[i,1]=="[-2.3,1.91]") {SSE.MM[i]=-0.195} # -2.3+((1.91+2.3)/2)
  if(MM.Pred.ChlaE.1[i,1]=="(1.91,3.39]") {SSE.MM[i]=2.65} # 1.91+((3.39-1.91)/2)
  if(MM.Pred.ChlaE.1[i,1]=="(3.39,6.41]") {SSE.MM[i]=4.9} # 3.39+((6.41-3.39)/2)
}

Square.MM= NULL
for (i in 1:1943){Square.MM[i]=(SSE.MM[i]-Evaluation.1[i+1,"Chla"])^2
}

sum(Square.MM, na.rm=TRUE)
#################################################################
# LATEX tables and figures
#################################################################
pdf("Plot-Comparison.pdf", width=8, height=10)
hist(Model.1[ ,"Chla"], breaks=100, xlim=c(-3,7), xlab="Chlorophyll a", main="")
axis(3, at= Breaks.I.Chla, col="black", col.ticks="black", col.axis="black", line=-2)
mtext("Interval", 3, line=-2, at=-3.5, col="black")
axis(3, at= Breaks.Q.Chla, col="royalblue4", col.ticks="royalblue4", col.axis="royalblue4", line=0)
mtext("Quantile", 3, line=0, at=-3.5, col="royalblue4")
axis(3, at= Breaks.MM.Chla, col="darkgreen", col.ticks="darkgreen", col.axis="darkgreen", line=2)
mtext("Moment Matching", 3, line=2, at=-3.5, col="darkgreen")
invisible(dev.off())
#################################################################
# LATEX OUTPUT by xtable 1.7-3 package
#################################################################
MM <- matrix(NA, 0, 5)

for(j in 1:10){
  assign(paste("MM",sep="")
         , rbind(get(noquote(paste("MM",sep="")))
                 , cbind(get(noquote(paste("logLik.MM.",j,sep="")))
                         , get(noquote(paste("AIC.MM.",j,sep="")))
                         , get(noquote(paste("BIC.MM.",j,sep="")))
                         , get(noquote(paste("MM.Accuracy.",j,sep="")))
                         , get(noquote(paste("AUC.MM.",j,sep=""))))))
}

colnames(MM) = c("logLik", "AIC", "BIC","Accuracy", "AUC")

MM <- rbind(MM, 
            cbind(mean(MM[,"logLik"]), mean(MM[,"AIC"]), mean(MM[,"BIC"]), mean(MM[,"Accuracy"]), mean(MM[,"AUC"])))


MM.Table <- xtable(MM, digits=3)
print(MM.Table)
#################################################################
# CPTs to LATEX N
print(xtable(fit.MM.1$N$prob), floating=FALSE)
# CPTs to LATEX P
print(xtable(fit.MM.1$P$prob), floating=FALSE)
# CPTs to LATEX Chla
print(xtable(fit.MM.1$Chla$prob[1:3, 1:3, 1]), floating=FALSE)
print(xtable(fit.MM.1$Chla$prob[1:3, 1:3, 2]), floating=FALSE)
print(xtable(fit.MM.1$Chla$prob[1:3, 1:3, 3]), floating=FALSE)
#################################################################
pdf("Plot-MM-Compare.pdf", width=8, height=10)

# Histogram of data for model building: continuous vs discretized
hist(Model.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a" , cex.axis=1.5)
hist(Model.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=c(-2.31, 0.6027629, 3.5081109, 6.42), main="",  add=T, border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Interval.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="black", cex=1.5)
mtext("Interval", 1,line=2,at=c(-3),col="black", cex=1.5)

# Histogram of data for model evaluation: continuous vs discretized
hist(Evaluation.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", cex.axis=1.5, xlab="Chlorophyll a")
hist(Evaluation.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=c(-2.31, 0.6027629, 3.5081109, 6.42), main="",  add=T, border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Interval.E.1[, "Chla"]), cex=1.5, 1,line=2,at=c(-1,2,5),col="black")
mtext("Interval", 1,line=2,at=c(-3), cex=1.5, col="black")
# Histogram of data for model building: continuous vs discretized
hist(Model.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", cex.axis=1.5, xlab="Chlorophyll a")
hist(Model.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=c(-2.3025851, 1.722767, 2.639057, 6.4134590), main="",  add=T, border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Quantile.1[, "Chla"]), cex=1.5, 1,line=2,at=c(-1,2,5),col="black")
mtext("Quantile", 1,line=2, cex=1.5, at=c(-3),col="black")

# Histogram of data for model evaluation: continuous vs discretized
hist(Evaluation.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", cex.axis=1.5, xlab="Chlorophyll a")
hist(Evaluation.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=c(-2.3025851, 1.722767, 2.639057, 6.4134590), main="",  add=T, border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(Quantile.E.1[, "Chla"]), cex=1.5, 1,line=2,at=c(-1,2,5),col="black")
mtext("Quantile", 1,line=2, cex=1.5,at=c(-3),col="black")

# Histogram of data for model building: continuous vs discretized
hist(Model.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", cex.axis=1.5, xlab="Chlorophyll a")
hist(Model.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=c(-2.3025851, 1.497835, 3.231148, 6.4134590), main="",  add=T, border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(MM.1[, "Chla"]), cex=1.5, 1,line=2,at=c(-1,2,5),col="black")
mtext("Moment Matching", 1, cex=1.5,line=2,at=c(-3),col="black")

# Histogram of data for model evaluation: continuous vs discretized
hist(Evaluation.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", cex.axis=1.5 , xlab="Chlorophyll a")
hist(Evaluation.1[ ,"Chla"], cex.lab=1.5, cex.axis=1.5, breaks=c(-2.3025851, 1.497835, 3.231148, 6.4134590), main="",  add=T, border="black",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(MM.E.1[, "Chla"]), cex=1.5, 1,line=2,at=c(-1,2,5),col="black")
mtext("Moment Matching", cex=1.5, 1,line=2,at=c(-3),col="black")
invisible(dev.off())
#################################################################
pdf("Plot-MM.pdf", width=8, height=10)
bn.fit.barchart(fit.MM.1$P, ylab="Phosphorus Levels")
bn.fit.barchart(fit.MM.1$N, ylab="Nitrogen Levels")
bn.fit.barchart(fit.MM.1$Chla, ylab="Chlorophyll a Levels")


# Histogram of data for model building: continuous vs discretized
hist(Model.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
hist(Model.1[ ,"Chla"], breaks=c(-2.3025851, 1.497835, 3.231148, 6.4134590), main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(MM.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Moment Matching", 1,line=2,at=c(-3),col="darkgreen")

# Histogram of data for model evaluation: continuous vs discretized
hist(Evaluation.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
hist(Evaluation.1[ ,"Chla"], breaks=c(-2.3025851, 1.497835, 3.231148, 6.4134590), main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(MM.E.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Moment Matching", 1,line=2,at=c(-3),col="darkgreen")

# Histogram of continuous data for model evaluation vs discretized model prediction
hist(Evaluation.1[ ,"Chla"], breaks=100, xlim=c(-3,7), freq = FALSE, main="", xlab="", border="royalblue4", xaxt = 'n')
axis(1, col.axis="royalblue4", , xlab="Chlorophyll a")
pred.1 <- runif(length(MM.Pred.ChlaE.1[MM.Pred.ChlaE.1=="[-2.3,2.56]"]), min = -2.3025851, max =  2.564949)
pred.2 <- runif(length(MM.Pred.ChlaE.1[MM.Pred.ChlaE.1=="(2.56,3.37]"]), min =  2.564949, max = 3.367296 )
pred.3 <- runif(length(MM.Pred.ChlaE.1[MM.Pred.ChlaE.1=="(3.37,6.41]"]), min =  3.367296, max = 6.4134590)
pred <- c(pred.1, pred.2, pred.3)
hist(pred, breaks=c(-2.3025851, 1.497835, 3.231148, 6.4134590), , main="",  add=T, border="darkgreen",  density=t(t(c(20,20,20))), angle=t(t(c(45,45,45))))
mtext(levels(MM.E.1[, "Chla"]), 1,line=2,at=c(-1,2,5),col="darkgreen")
mtext("Moment Matching", 1,line=2,at=c(-3),col="darkgreen")
invisible(dev.off())
#################################################################
# Save History
#################################################################
base <- "~/Dropbox/Research/Discretization/"
setwd(base)
# save your command history 
savehistory(file="11252015") # default is ".Rhistory" 
# recall your command history 
loadhistory(file="11252015") # default is ".Rhistory"