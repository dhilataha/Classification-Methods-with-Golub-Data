#Chapter 8
#1
library(multtest);data(golub); gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML")) 
maxgol <- apply(golub[,gol.fac=="ALL"], 1, function(x) max(x)) 
mingol <- apply(golub[,gol.fac=="AML"], 1, function(x) min(x))
sum(maxgol < mingol) 
which.min(maxgol - mingol)
golub.gnames[2124,] 
boxplot(golub[2124,] ~gol.fac)

#plot 1
gol.rp <- rpart(gol.fac ~ golub[2124,], method="class", cp=0.001) 
plot(gol.rp, branch=0,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)

#plot 2
grep("Gdf5",golub.gnames[,2]) 
gol.rp <- rpart(gol.fac ~ golub[2058,], method="class", cp=0.001)
plot(gol.rp, branch=0,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)

#plot 3
gol.rp <- rpart(gol.fac ~., data.frame(t(golub)), method="class", cp=0.001)
plot(gol.rp, branch=0,margin=0.1); text(gol.rp, digits=3, use.n=TRUE)




#2a
library(multtest);library(ROCR);data(golub) 
golub.clchanged <- -golub.cl +1 
pred <- prediction(golub[1042,], golub.clchanged) 
perf <- performance(pred, "sens", "spec")
plot(perf) 

#3a
library(rpart) 
predictors <- matrix(rnorm(100*4,0,1),100,4)
colnames(predictors) <- letters[1:4] 
groups <- gl(2,50)
simdata <- data.frame(groups,predictors)
rp<-rpart(groups ~ a + b + c + d,method="class",data=simdata) 
predicted <- predict(rp,type="class")
table(predicted,groups)
plot(rp, branch=0,margin=0.1); text(rp, digits=3, use.n=TRUE)

library(e1071) 
svmest <- svm(predictors, groups, data=df, type = "C-classification", kernel = "linear") 
svmpred <- predict(svmest, predictors, probability=TRUE)
table(svmpred, groups) 

library(nnet) 
nnest <- nnet(groups ~ ., data = simdata, size = 5,maxit = 500, decay = 0.01, MaxNWts = 5000) 
pred <- predict(nnest, type = "class") 

table(pred, groups) # prints confusion matrix

#4
library(ALL); library(hgu95av2.db); library(rpart); data(ALL) 
ALLrem <- ALL[,which(pData(ALL)$remission %in% c("CR","REF"))] 
remfac <-factor(pData(ALLrem)$remission) 
pano <- apply(exprs(ALLrem),1,function(x)
  t.test(x ~ remfac)$p.value)
names <- featureNames(ALLrem)[pano<.001] 
ALLremsel<- ALLrem[names,] 
data <- data.frame(t(exprs(ALLremsel)))
all.rp <- rpart(remfac ~., data, method="class", cp=0.001) 
plot(all.rp, branch=0,margin=0.1); text(all.rp, digits=3, use.n=TRUE)
rpart.pred <- predict(all.rp, type="class")
table(rpart.pred,remfac) 
7/(93+1+6+14)
mget(c("1840_g_at","36769_at","1472_g_at","854_at"), env = hgu95av2GENENAME)


#5
library(ROCR); data(golub, package = "multtest") 
gol.true <- factor(golub.cl,levels=0:1,labels= c("TRUE","FALSE")) 
auc.values <- apply(golub,1, 
    function(x) performance(prediction(x, gol.true),"auc")@y.values[[1]]) 
o <- order(auc.values,decreasing=TRUE) 
golub.gnames[o[1:25],2]

ade4TkGUI()
