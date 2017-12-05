library("e1071")
library(readr)


## piece of code of interest, extracted from the slides
#data(iris)
#model <-naiveBayes(Species ~., data=iris)
#table(predict(model, iris), iris[,5])
#pairs(iris[1:4], main = "Edgar Anderson's Iris Data", pch = 21, bg = c("red", "green3", "blue")[unclass(iris$Species)])

## Load the datasets for training and testing
train_df <- read_csv("bioinfo1.train.expr.csv")
test_df <- read_csv("bioinfo1.test.expr.csv")

## load the labels, and define the x columns as factors, do not really what this does, but is neede
## just googled the error and found it on stackoverflow
labels_df <- read_csv("bioinfo1.train.labels.csv")
labels_df$x <- as.factor(labels_df$x)


## Define the naiveBayes model based on the formula and the data
## the formula is just the class associated to each element of the train dataset
model <-naiveBayes(labels_df$x ~., data=train_df)

## Now we can predict over the same dataset, in order to measure the accuracy of our model
prediction <- predict(model,train_df)
summary(prediction)
## Ouptuts the confussion matrix for the prediction
table(prediction, labels_df$x)
cat("Accuracy of naive bayes", sum(prediction!=labels_df$x, na.rm=TRUE) / length(prediction) * 100, "%")


## Let's predict over the test dataset
real_prediction <- predict(model, test_df)
show(real_prediction)
write.csv(real_prediction, file = "output.csv", row.names=FALSE)



library(randomForest)
##let's try with random forests
rf_model = randomForest(labels_df$x ~., data=train_df)
rf_prediction <- predict(rf_model,train_df)
table(rf_prediction, labels_df$x)
#summary(rf_model)
# real prediction for rf model, and let's assume it is also 100% of accuracy
rf_real_prediction <-predict(rf_model, test_df)

# now we can show how the naive bayes prediction for the test data performs
table(real_prediction, rf_real_prediction)

cat("Accuracy of naive bayes", sum(prediction!=rf_prediction, na.rm=TRUE) / length(prediction) * 100, "%")
## Lets try SVM
#svm_model <- svm(labels_df$x ~., data=train_df)
#svm_prediction <- predict(svm_model,train_df)
#table(svm_prediction, labels_df$x)
#summary(svm_model)
