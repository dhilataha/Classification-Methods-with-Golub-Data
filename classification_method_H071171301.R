data(golub, package = "multtest") 
gol.true <- factor(golub.cl,levels=0:1,labels= c(" ALL","not ALL"))
gol.pred <- factor(golub[1042,]>1.27,levels=c("TRUE","FALSE"), labels=c("ALL","notALL")) 
table(gol.pred,gol.true)

#Ex 2
library(ROCR) 
gol.true <- factor(golub.cl,levels=0:1,labels= c("TRUE","FALSE")) 
pred <- prediction(golub[1042,], gol.true)
perf <- performance(pred, "tpr", "fpr" )
performance(pred,"auc")
plot(perf)


#8.3 Classification Trees
#Ex 1
set.seed(123); n<-10 ; sigma <- 0.5
fac <- factor(c(rep(1,n),rep(2,n),rep(3,n)))
levels(fac) <- c("ALL1","ALL2","AML")
geneA <- c(rnorm(10,0,sigma),rnorm(10,2,sigma),rnorm(10,4,sigma))
dat <- data.frame(fac,geneA)
library(rpart) 
rp <- rpart(fac ~ geneA, method="class",data=dat)
plot(rp, branch=0,margin=0.1); text(rp, digits=3, use.n=TRUE) 

#Ex 2
set.seed(123)
n<-10 ; sigma <- 0.5 
fac <- factor(c(rep(1,n),rep(2,n),rep(3,n)))
levels(fac) <- c("ALL1","ALL2","AML") 
geneA <- c(rnorm(20,0,sigma),rnorm(10,2,sigma))
geneB <- c(rnorm(10,0,sigma),rnorm(20,2,sigma))
geneC <- c(rnorm(30,1,sigma)) 
dat <- data.frame(fac,geneA,geneB,geneC) 
library(rpart) 
rp <- rpart(fac ~ geneA + geneB + geneC, method="class",data=dat)

#Ex 3
library(rpart);data(golub); library(multtest) 
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML")) 
gol.rp <- rpart(gol.fac ~ golub[1042,] , method="class")
predictedclass <- predict(gol.rp, type="class")
table(predictedclass, gol.fac) 
predict(gol.rp,type="class")
summary(gol.rp) 

predict(gol.rp,type="class")
predict(gol.rp, type="prob") 

boxplot(golub[2124,] ~gol.fac)

#Ex 4
library(rpart);data(golub); library(multtest)
row.names(golub)<- paste("gene", 1:3051, sep = "") 
goldata <- data.frame(t(golub[1:3051,])) 
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML")) 
gol.rp <- rpart(gol.fac~., data=goldata, method="class", cp=0.001) 
plot(gol.rp, branch=0,margin=0.1); text(gol.rp, digits=3, use.n=TRUE) 
golub.gnames[896,]


#Ex5
library("hgu95av2.db");library(ALL);data(ALL) 
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")] 
pano <- apply(exprs(ALLB123), 1, function(x) anova(lm(x ~ ALLB123$BT))$Pr[1]) 
names <- featureNames(ALL)[pano<0.000001] 
symb <- mget(names, env = hgu95av2SYMBOL) 
ALLBTnames <- ALLB123[names, ] 
probedat <- as.matrix(exprs(ALLBTnames)) 
row.names(probedat)<-unlist(symb) 

diagnosed <- factor(ALLBTnames$BT) 
tr <- rpart(factor(ALLBTnames$BT) ~ ., data = data.frame(t(probedat)))
plot(tr, branch=0,margin=0.1); text(tr, digits=3, use.n=TRUE)
rpartpred <- predict(tr, type="class") 
table(rpartpred,diagnosed) 

predicted.class <- predict(tr, type="class") 
predicted.probabilities <- predict(tr, type="prob") 
out <- data.frame(predicted.probabilities,predicted.class, diagnosis=factor(ALLBTnames$BT)) 
print(out,digits=2) 

#Ex 6
i <- sample(1:78, 39, replace = FALSE) 
noti <- setdiff(1:78,i) 
df <- data.frame(Y = factor(ALLBTnames$BT), X =t(probedat)) 
rpart.est <- rpart(Y ~ ., data = df, subset=i) 
rpart.pred.t <- predict(rpart.est, df[i,], type="class") 
table(rpart.pred.t,factor(ALLBTnames$BT[i])) 
rpart.pred.t B1 B2 B3 
rpart.pred.v <- predict(rpart.est,df[noti,], type="class") 
table(rpart.pred.v,factor(ALLBTnames$BT[noti])) 


#8.4 Support Vector Machine
#Ex 1
library(e1071) 
df <- data.frame(Y = factor(ALLBTnames$BT), X =t(probedat)) 
Y <- factor(ALLBTnames$BT);X <- t(probedat) 
svmest <- svm(X, Y, data=df, type = "C-classification", kernel = "linear") 
svmpred <- predict(svmest, X, probability=TRUE) 
table(svmpred, factor(ALLBTnames$BT)) 

summary(svmest) 
dim(svmest$SV)
dim(svmest$coefs)

#Ex 2
Yt <- factor(ALLBTnames$BT)[i]; Yv <- factor(ALLBTnames$BT)[noti]
X <- t(probedat); Xt <- X[i,]; Xv <- X[noti,] 
svmest <- svm(Xt, Yt, type = "C-classification", kernel = "linear") 
svmpredt <- predict(svmest, Xt, probability=TRUE) 
table(svmpredt, Yt) 

svmpredv <- predict(svmest, Xv, probability=TRUE) 
table(svmpredv, Yv) 


#8.5 Neural Network
#Ex 1
Y <- factor(ALLBTnames$BT);X <- t(probedat)
library(nnet) 
df <- data.frame(Y = Y, X = X[, sample(ncol(X), 20)]) 
nnest <- nnet(Y ~ .,data = df, size = 5, maxit = 500, decay = 0.01, + MaxNWts = 5000)
pred <- predict(nnest, type = "class") 
table(pred, Y) # prints confusion matrix

#Ex 2
nnest.t <- nnet(Y ~ ., data = df,subset=i, size = 5,decay = 0.01, + maxit=500) 
prednnt <- predict(nnest.t, df[i,], type = "class") 
table(prednnt,Ytrain=Y[i]) 
prednnv <- predict(nnest.t, df[noti,], type = "class") 
table(prednnv, Yval= Y[noti]) 


#8.6 Generalized Liniear Models
#Ex 1
library(faraway) 
logitmod <- glm((-golub.cl + 1) ~ golub[1042,], family=binomial(link = "logit")) 
pchisq(deviance(logitmod),df.residual(logitmod),lower=FALSE) 
plot((-golub.cl + 1) ~ golub[1042,], xlim=c(-2,5), ylim = c(0,1), xlab="CCND3 expression values ", ylab="Probability of ALL") 
x <- seq(-2,5,.1) 
lines(x,ilogit(-4.844124 + 4.439953*x)) 
pchisq(deviance(logitmod),df.residual(logitmod),lower=FALSE)
summary(logitmod)

pred <- predict(logitmod,type="response") > 0.5
pred.fac <- factor(pred,levels=c(TRUE,FALSE),labels=c("ALL","not ALL")) 
table(pred.fac,gol.fac) 

#Ex 2
library(nnet);library("hgu95av2.db");library(ALL);data(ALL)
probe.names <- c("1389_at","35991_at","40440_at") 
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")] 
probedat <- exprs(ALLB123)[probe.names,] 
row.names(probedat) <- unlist(mget(probe.names, env = hgu95av2SYMBOL)) 
fac <- factor(ALLB123$BT,levels=c("B1","B2","B3")) 
dat <- data.frame(fac,t(probedat)) 
mnmod <- multinom(fac ~ ., family=binomial(link = "logit"),data=dat)
summary(mnmod)

predmn <- predict(mnmod,type="class") 
table(predmn,fac) 






