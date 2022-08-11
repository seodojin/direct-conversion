# caret 
# brief introduction 

install.packages("caret", dependencies = c("Depends", "Suggests"))
install.packages("caret")
library(caret)

library(mlbench)
data(Sonar)
set.seed(107)

inTrain = createDataPartition(y = Sonar$Class, ## target variable
                              p = 0.75, ## 75% -> training data 
                              list = FALSE) # format of the results 
str(inTrain) # a set of integers for the rows that belong to the training dataset 

training = Sonar[inTrain, ]
testing = Sonar[-inTrain, ]

ctrl = trainControl(method = "repeatedCV", # resampling method
                    number = 10,  # number of fold
                    repeats = 3,  # number of repeats 
                    summaryFunction = twoClassSummary, # performance measures
                    classProbs = T) # in case of two class problems, to calcular probabilities

plsfit = train(Class ~ ., 
               data = training, 
               method = "pls", 
               tuneLength = 15, 
               # a set of tuning parameter values, all integers between 1 and 15
               # tuneGrid for specific values are desired 
               trControl = ctrl,
               metric = "ROC", # performance measure to select the optimal model
               preProc = c("center", "scale")
                 )

plot(plsfit)

plsClass = predict(plsfit, newdata = testing)
plsProbs = predict(plsfit, newdata = testing, type = "prob")

confusionMatrix(plsClass, testing$Class)

# regularized discriminant analysis 

ldafit = train(Class ~ ., 
               data = training, 
               method = "lda", 
               # a set of tuning parameter values, all integers between 1 and 15
               # tuneGrid for specific values are desired 
               trControl = ctrl,
               metric = "ROC", # performance measure to select the optimal model
               preProc = c("center", "scale")
               )

ldaClass = predict(ldafit, newdata = testing)
ldaProbs = predict(ldafit, newdata = testing, type = "prob")

confusionMatrix(ldaClass, testing$Class)

# Compared models 

resamps = resamples(list(pls = plsfit, lda = ldafit))
summary(resamps)

diffs = diff(resamps) # paired t-test for each resample 
summary(diffs)

bwplot(resamps)

