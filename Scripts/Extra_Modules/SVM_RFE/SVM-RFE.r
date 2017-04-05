library(e1071)

################################################
# Feature Ranking with SVM-RFE
################################################
svmrfeFeatureRanking = function(x,y){
    n = ncol(x)
    
    survivingFeaturesIndexes = seq(1:n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    
    while(length(survivingFeaturesIndexes)>0){
        #train the support vector machine
        svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" )
        
        #compute the weight vector
        w = t(svmModel$coefs)%*%svmModel$SV
        
        #compute ranking criteria
        rankingCriteria = w * w
        
        #rank the features
        ranking = sort(rankingCriteria, index.return = TRUE)$ix
        
        #update feature ranked list
        featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
        rankedFeatureIndex = rankedFeatureIndex - 1
        
        #eliminate the feature with smallest ranking criterion
        (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
        
    }
    
    return (featureRankedList)
}
