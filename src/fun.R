library(glmnet)
library(dplyr)
load("./example/SBS_signature.RData")
load("./example/actual_lm_regression_result.RData")

pred_TME = function(data, cancer_type) {

    # data preprocessing

    # convert
    x = xxx(data)

    # predict
    y = predict(model, newdata = x)

    # out
    out = data.frame(y)

    return(out)

}