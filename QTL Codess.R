setwd("C:/R")
QTL_group_6 <- read.csv("C:/R/QTL_group_6.csv", header=FALSE)
Data <- as.data.frame(read.table(file = "C:/R/QTL_group_6.csv", header = TRUE))
library(qtl)
Data <-read.cross("csv",".","QTL_group_6.csv")
summary(Data)
nind(Data)
nphe(Data)
nchr(Data)
totmar(Data)
nmar(Data)
plot(Data)
Data
plotMissing(Data)
plotMap(Data)
plotPheno(Data, pheno.col=2)
plotMap(Data, chr=c(1,3), show.marker.names=TRUE)
plotMissing(Data, reorder=TRUE)
Data <-drop.nullmarkers(Data)
totmar(Data)
Data <-est.rf(Data)
plotRF(Data)
plotRF(Data, chr=c(1,3))
summary(Data)
Data <- drop.markers(Data, "Gg_rs15777012")
summary(Data)
totmar(Data)
newmap <-est.map(Data, error.prob=0.01)
plotMap(Data, newmap,show.marker.names = TRUE)
Data <-replace.map(Data, newmap)
Data <-calc.errorlod(Data, error.prob=0.01)
top.errorlod(Data)
plotGeno(Data, chr=3, ind=c(24:34, 71:81))

Comb <- pull.pheno(Data,"comb_g")
sex <- pull.pheno(Data,"sex")
comb_sex <- cbind(Comb,sex)

weight <- pull.pheno(Data, "w212")
weight_sex <-cbind(weight, sex)

head(weight_sex)
##QTL mapping
Data <-calc.genoprob(Data, step=5, error.prob=0.01)

head(weight_sex)

out.hk <-scanone(Data, method="hk", pheno.col = Comb)
summary(out.hk)

Data <-sim.geno(Data, step=5, n.draws=16, error.prob=0.01)
out.imp <-scanone(Data, method="imp", pheno.col = Comb)

summary(out.imp)


summary(out.hk, threshold=3)


summary(out.imp, threshold=3)
plot(out.hk, chr=c(1,3), add=TRUE)
plot(out.hk, chr=c(1,3), col="blue", add=TRUE)

operm.hk <-scanone(Data, method="hk", n.perm=1000)
summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE)
qc <-c(3)
qp <-c(13)
qtl <-makeqtl(Data, chr=qc, pos=qp)

myformula <-y  ̃Q1
out.fq <-fitqtl(Data, qtl=qtl, formula = myformula, pheno.col = Comb)
summary(out.fq)
out.fq <-fitqtl(Data, qtl=qtl, formula = myformula, get.ests=TRUE, pheno.col = Comb)
summary(out.fq)
revqtl <-refineqtl(Data, qtl=qtl, formula = myformula, pheno.col = Comb)
summary(revqtl)
plot(revqtl)
out.fq2 <-fitqtl(Data, qtl=revqtl, formula=myformula, pheno.col = Comb)
summary(out.fq2)
myformula2 <- y~Q1+sex

myformula4 <- y~ Q1+weight+sex
out.fq3 <-fitqtl(Data, qtl=qtl, formula = myformula2, covar = comb_sex, pheno.col = Comb)
summary(out.fq3)
out.fq4 <-fitqtl(Data, qtl=qtl, formula = myformula4, covar = weight_sex, pheno.col = Comb)
summary(out.fq4)

