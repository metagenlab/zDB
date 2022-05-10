library('e1071')
load('BPBAac.Rdata')
testdata <- read.csv('sample.data',head=TRUE)
testall <- subset(testdata,select=c(-Name,-Name0))
Nametestall <- subset(testdata,select=c(Name))
test3 <- testall
NameTest <- Nametestall
pred <- predict(model,test3,decision.values = TRUE)
values<-attr(pred,"decision.values")
Lst <- list(NameTest,values)
write.table(Lst,file="values.csv")
system("perl Classify.pl 0.5")


