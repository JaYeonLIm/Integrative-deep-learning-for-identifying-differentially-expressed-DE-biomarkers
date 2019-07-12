# Integrative-deep-learning-for-identifying-differentially-expressed-DE-biomarkers
## Description
This package is aim to find biomarkers via integrative deep learning. We propose adding an integrated layer to the deep learning structure, which would enable the effective analysis of genetic data as well as the discovery of significant biomarkers of disease. `INTEGRATED_GANN(DList, y, number_of_function, out_fct,
                       err_fct,lrate, number_lambdas,
                       lambda_max, weights_on_input_layer, epsilon_lambda,
                       maxiter = 500, ms, epsilon_iterations = 1e-5)` is main fucntion where we train datasets and test the performance of trained model.

## Arguments

`DList` : a list. It consists of multiple gene expression datasets. Each dataset contains same dimension. We use the function `PREDICTOR(DList)` to make the list to predictor of models. 

`y`: a matrix. It contains clinical data which consisits of 1 or 0 as the model aims at classification.

`number_of_function`:	the number of activation function in each integrated X 

`out_fct`	: options - "ReLU", "linear", "logistic"

`err_fct` :	Regression - "sse" , Classification - "ce"

`lrate` : learning rate	

`number_lambdas`:	the number of lamdas

`epsilon_lambda` :	the minimum value of lambda

`maxiter`:	maximum iterations

`ms` :	"cv" - cross validation

`epsilon_iterations` : threshold in training

## Examples 

First, we need to run functions which are components of main function `INTEGRATED_GANN`. 
```{r}
library(coda)
library(MASS)
library(foreach)
library(iterators)
library(parallel)
library(doMC)
library(e1071)
library(MCMCpack)
library(penalized)
library(glmnet)

PREDICTOR = function(DList){
  tmp = c()
  for(i in 1:length(DList)){
    tmp = rbind(tmp,DList[[i]])
  }
  x = t(tmp)
  return(x)
}


INTEGRATE_INDEX=function(DList){
  G = nrow(DList[[1]])
  K = length(DList)
  integration_index = list()
  for (i in 1:G) {
    integration_index[[i]] =c(i)
    for (k in 2:K){
      integration_index[[i]]= append(integration_index[[i]],i+(k-1)*G)
    }
  }
  return(integration_index)
}


SETUPd = function (x, y, number_of_function, out_fct, err_fct,
                   integration_index, lrate, weights_on_input_layer,
                   maxiter, epsilon_iterations = 1e-5)
{
  
  sm = list()
  sm$bias = 0
  sm$v = list()  
  sm$w = list()  
  sm$beta = list() 
  sm$n = length(y) 
  sm$K = length(integration_index)  
  
  sm$integration_index = integration_index 
  sm$number_of_function = number_of_function  
  sm$out_fct = out_fct  
  sm$err_fct = err_fct  
  sm$lrate = lrate     
  sm$maxiter = maxiter  
  sm$weights_on_input_layer = weights_on_input_layer  
  sm$sse = rep(0, maxiter) 
  sm$iteration = 0
  sm$index_active = c(1 : sm$K)
  sm$epsilon_iterations = epsilon_iterations 
  sm$lambda_list = list()
  sm$beta_C = list()
  sm$youden = vector()
  sm$youdenin = list()
  return(sm)
}

INTEGRATE_PREDICTORS = function(x, integration_index)  
{
  N = nrow(x)
  integratedX = list()
  for (i in 1 : length(integration_index))
    integratedX[[i]] = matrix(x[, integration_index[[i]]], N, length(integration_index[[i]]))
  return(integratedX)
}

DATA_SPLIT = function(response, predictors, ms, train_rate = 0.7)      
{
  number_integratedX = length(predictors)   
  if (ms == "cv")
  {
    
    data_size = length(response)
    train_index = sort(sample(data_size, train_rate * data_size)) 
    validation_index = (1 : data_size)[-train_index]              
    train = list(response = response[train_index])               
    validation = list(response   = response[validation_index])    
    for (i in 1 : number_integratedX)
    {
      if (NCOL(predictors[[i]]) == 1)                   
      {
        train$integratedX[[i]] = matrix(predictors[[i]][train_index, ])   
        validation$integratedX[[i]] = matrix(predictors[[i]][validation_index, ]) 
      }
      else
      {
        train$integratedX[[i]] = predictors[[i]][train_index, ]      
        validation$integratedX[[i]] = predictors[[i]][validation_index, ]  
      }
    }
    return(list(train = train, validation = validation))
  }
  else
  {
    train = list(response = response, integratedX = predictors)
    validation = train
    return(list(train = train, validation = validation))
  }
}

LAMBDAS_GIVEN_LAMBDA_MAX = function(lambda_max, number_lambdas,
                                    epsilon_lambda)     
{
  lambda_max = lambda_max    
  lambdas_all = rep(0, number_lambdas) 
  lambda_min = epsilon_lambda
  ratio_max_min = 1.0 / epsilon_lambda
  lambdas_all = seq(lambda_min,lambda_max,(lambda_max-lambda_min)/(number_lambdas-1))
  return(lambdas_all)
}

FUNCTION_SELECTION = function(sm)
{
  RELU = list()
  RESULT_FUNCTIONS_NONRANDOM = list()
  
  K = sm$K
  number = sm$number_of_function
  
  RELU[[1]] = function(x){ log(1.0 + exp(x)) }
  
  ALL_FUNCTIONS = list(RELU)
  
  FUNCTION_EACH_PREDICTOR = list()
  FUNCTION_EACH_PREDICTOR[[1]] = ALL_FUNCTIONS[[1]]
  RESULT_FUNCTIONS_NONRANDOM = FUNCTION_EACH_PREDICTOR
  sm$selected_function = RESULT_FUNCTIONS_NONRANDOM
  return(sm)
}

SYSTEMATIC_INITIAL_WEIGHT_GANN = function(sm)
{
  
  K = sm$K   
  y = sm$y   
  sm$n=length(y)
  number = sm$number_of_function
  integration_index = sm$integration_index  
  selected_function = sm$selected_function  
  C = matrix(0, sm$n, K)
  
  consolidated_values = rep(0, sm$n) 
  variance = 1
  for(i in 1 : K)
  {
    dimension = length(integration_index[[i]])   
    sm$v[[i]] = matrix(1, nrow = dimension, ncol = number[i])  
    sm$w[[i]] = matrix(1, number[i], 1)                       
    
    consolidated_values = consolidated_values + selected_function[[1]][[1]](
      apply(sm$integratedX[[i]], 1, sum))
    C[, i] = consolidated_values  
    sm$beta_C[[i]]=0
  }
  design_matrix = as.matrix(data.frame(y, C)) 
  fit <- glmnet(x=design_matrix, y,alpha = 0,lambda = 0)
  ridge_coef = coef(fit)
  sm$beta = matrix(ridge_coef[3 : length(ridge_coef)], length(ridge_coef) - 2, 1, dimnames = NULL) 
  sm$bias = ridge_coef[1]
  return(sm)
}

FEEDFOWARD=function(sm, lambda, index_active)
{
  x = sm$integratedX  
  y = sm$y
  v = sm$v
  w = sm$w
  normalized_w = sm$w
  beta = sm$beta  
  bias = sm$bias 
  n = length(sm$y)
  K = sm$K   
  out_fct = sm$out_fct
  err_fct = sm$err_fct
  number = sm$number_of_function  
  selected_function = sm$selected_function  
  consolidated_value = matrix(0, n, K)  
  if (out_fct == "linear")
  {
    output_function = function(x)(x)
  }
  else if (out_fct == "ReLU")
  {
    output_function = function(x)
    {
      log(1.0 + exp(x))
    }
  }
  else
  {
    if (out_fct == "logistic")
    {
      output_function = function(x)
      {
        1 / (1 + exp(-x))
      }
    }
  }
  
  if (err_fct == "ce")   
  {
    error_function = function(x, Y)
    {
      
      -(Y * log(1/(1+exp(-x))) + (1 - Y) * log(1 - (1/(1+exp(-x)))))  
    }
  }
  else
  {
    if (err_fct == "sse")
    {
      error_function = function(x, y)
      {
        (x - y)^2/2
      }
    }
  }
  
  if (sm$weights_on_input_layer == TRUE)  
  {
    for (j in 1 : length(index_active))
    {
      i = index_active[j]
      activation = matrix(0, n, number[i])  
      sm$a[[i]] = x[[i]] %*% v[[i]]   
      for (m in 1 : number[i])
      {
        activation[, m] = selected_function[[1]][[1]](sm$a[[i]][, m])  
      }
      sm$O[[i]] = activation  
      normalized_w[[i]] = w[[i]] / sqrt(sum(w[[i]]^2))   
      consolidated_value[, i] = sm$O[[i]] %*% normalized_w[[i]] 
    }
    yhat = consolidated_value %*% beta + bias  
  }
  #---------------------------------------------------------------------------------
  if (sm$weights_on_input_layer == FALSE)   
  {
    for(j in 1 : length(index_active))
    {
      i = index_active[j]
      activation = matrix(0,n,number[i])
      predictors_summation = matrix(apply(x[[i]], 1, sum))  
      sm$a[[i]] = matrix(rep(predictors_summation, times = number[i]), nrow = n,
                         ncol = number[i])      
      for (m in 1 : number[i])
        activation[, m] = selected_function[[1]][[1]](predictors_summation)   
      sm$O[[i]] = activation   
      normalized_w[[i]] = w[[i]] / sqrt(sum(w[[i]]^2))
      consolidated_value[, i] = sm$O[[i]] %*% normalized_w[[i]] 
    }
    yhat = consolidated_value %*% beta + bias
  }
  #---------------------------------------------------------------------------------
  sm$yhat = output_function(yhat)         
  sm$loss_value = mean(error_function(sm$yhat, y)) + lambda * sum(abs(beta))
  sm$output_function = output_function
  sm$loss_function = error_function
  sm$C = consolidated_value
  sm$normalized_w = normalized_w
  return(sm)
}


differentiate <-function (orig.fct, hessian = FALSE)
{
  body.fct <- deparse(body(orig.fct))
  if (body.fct[1] == "{")
    body.fct <- body.fct[2]
  text <- paste("y~", body.fct, sep = "")
  text2 <- paste(deparse(orig.fct)[1], "{}")
  temp <- deriv(eval((parse(text = text))), "x", func = eval(parse(text = text2)),
                hessian = hessian)
  temp <- deparse(temp)
  derivative <- NULL
  #---------------------------------------------------------------------------------
  if (!hessian)
    for (i in 1:length(temp))
    {
      if (!any(grep("value", temp[i])))
        derivative <- c(derivative, temp[i])
    } else
      for (i in 1:length(temp)) {
        if (!any(grep("value", temp[i]), grep("grad", temp[i]), grep(", c", temp[i])))
          derivative <- c(derivative, temp[i])
      }
  number <- NULL
  for (i in 1:length(derivative)) {
    if (any(grep("<-", derivative[i])))
      number <- i
  }
  #---------------------------------------------------------------------------------
  if (is.null(number))
  {
    return(function(x)
    {
      matrix(0, nrow(x), ncol(x))
    })
  }
  else
  {
    derivative[number] <- unlist(strsplit(derivative[number], "<-"))[2]
    derivative <- eval(parse(text = derivative))
  }
  return(derivative)
}

#------------------------------------------------------------------------------------

BACKPROP_VW = function(sm)
{
  x = sm$integratedX
  C = sm$C
  O = sm$O
  a = sm$a
  K = sm$K
  n = length(sm$y)
  number = sm$number_of_function
  selected_function = sm$selected_function
  v = sm$v
  w = sm$w
  index_active = sm$index_active
  beta = sm$beta
  bias = sm$bias
  lrate = sm$lrate
  output_function = sm$output_function
  loss_function = sm$loss_function
  loss_function_deriv = differentiate(loss_function)
  output_function_deriv = differentiate(output_function)
  C_deriv = matrix(1, n, K)
  if (sm$out_fct == "linear")
    pderiv_yhat = loss_function_deriv(sm$yhat, sm$y) 
  else
  {
    C_deriv = output_function_deriv(C)
    pderiv_yhat = loss_function_deriv(sm$yhat, sm$y) * output_function_deriv(sm$yhat)  
  }
  
  #---------------------------------------------------------------------------------
  if (sm$weights_on_input_layer == TRUE)
  {
    for (i in 1 : length(index_active))
    {
      j = index_active[i]
      a_deriv = matrix(0, n, number[j])
      error_hidden = matrix(0, n, number[j])
      error_consolidated = matrix(0, n, 1)
      error_consolidated = beta[j] * pderiv_yhat * C_deriv[, j]
      wt_w = sum(w[[j]]^2)
      w_norm = sqrt(wt_w)
      for (l in 1 : number[j])
        a_deriv[, l] = differentiate(selected_function[[1]][[1]])(a[[j]][, l])
      error_hidden = (error_consolidated %*% t(w[[j]] / w_norm)) * a_deriv
      pderiv_v = t(x[[j]]) %*% error_hidden
      v[[j]] = v[[j]] - lrate * pderiv_v / n
      deriv_of_w_wrt_c = (wt_w^-0.5) * O[[j]] -
        (wt_w^-1.5) * (O[[j]] %*% (w[[j]]) %*% t(w[[j]]))
      pderiv_w = t(deriv_of_w_wrt_c) %*% error_consolidated
      w[[j]] = w[[j]] - lrate * pderiv_w / n
    }
  }
  #---------------------------------------------------------------------------------
  if (sm$weights_on_input_layer == FALSE)
  {
    for (i in 1 : length(index_active))
    {
      j = index_active[i]
      error_consolidated = matrix(0, n, 1)
      error_consolidated = beta[j] * pderiv_yhat * C_deriv[, j]
      wt_w = sum(w[[j]]^2)
      deriv_of_w_wrt_c = (wt_w^-0.5) * O[[j]] -
        (wt_w^-1.5) * (O[[j]] %*% (w[[j]]) %*% t(w[[j]]))
      pderiv_w = t(deriv_of_w_wrt_c) %*% error_consolidated
      w[[j]] = w[[j]] - lrate * pderiv_w / n
    }
  }
  sm$v = v
  sm$w = w
  sm$bias = bias
  return(sm)
}

#------------------------------------------------------------------------------------
BACKPROP_BETA = function(sm, lambda)
{
  x = sm$integratedX
  C = sm$C
  O = sm$O
  a = sm$a
  K = sm$K
  NR_iterations = 500   
  halving_iterations = 50
  n = length(sm$y)
  new_yhat = sm$yhat
  old_R = sm$loss_value
  number = sm$number_of_function
  selected_function = sm$selected_function
  tilde_beta = beta = sm$beta
  tilde_bias = bias = sm$bias
  output_function = sm$output_function
  loss_function = sm$loss_function
  loss_function_deriv = differentiate(loss_function)   
  output_function_deriv = differentiate(output_function)
  #---------------------------------------------------------------------------------
  scale = 1e-4
  bias = tilde_bias - lrate*mean(new_yhat - sm$y)
  for (j in 1 : K)
  {
    for(NR_iter in 1 : NR_iterations)
    { 
      d = mean((new_yhat - sm$y) * C[, j]) + lambda * sign(tilde_beta[j])
      tmp_power = (tilde_beta[j]) #/ scale
      
      if(exp(tmp_power)==Inf){
        cdf = 1
      } else {
        cdf = exp(tmp_power) / (1+exp(tmp_power))  
      }
      
      w = (1/scale) * cdf * (1 - cdf) 
      d2 = mean(C[, j]^2) + lambda * w 
      diff = d / d2
      for (iters in 1 : halving_iterations)
      {
        beta[j] = tilde_beta[j] - diff
        new_yhat = C %*% beta + bias
        new_R = mean(loss_function(new_yhat, sm$y)) + lambda * sum(abs(beta))
      }
      if (abs(beta[j] - tilde_beta[j]) < 1e-6){
        break
      }
      tilde_beta[j] = beta[j]
    }
    if (abs(beta[j]) < 1e-5)
      beta[j] = 0
  }
  sm$beta = beta
  sm$bias = bias
  sm$index_active = which(beta != 0)
  return(sm)
}

#------------------------------------------------------------------------------------

TRAIN_GANN = function(sm, lambda)
{
  lrate_halve = sm$lrate     
  threshold = sm$epsilon_iterations    
  sm$warning = "Algorithm converges"  
  sm$lambda = lambda   
  initial_index_active = sm$index_active  
  cat("Active index = ", sm$index_active, "\n")
  if (length(initial_index_active) == 0)  
  {
    cat("All betas are shinked to 0", "\n")
    return(sm)
  }
  D = length(initial_index_active)
  
  for(i in 1 : sm$maxiter)
  {
    sm = FEEDFOWARD(sm, lambda, sm$index_active)  
    sm$sse[i] = sm$loss_value      
    if (i == 1)
      sse_difference = sm$sse[i]
    else
      sse_difference = sm$sse[i-1] - sm$sse[i]
    if(abs(sse_difference) < threshold) break   
    sm = BACKPROP_BETA(sm, lambda) 
    sm = FEEDFOWARD(sm, lambda, sm$index_active)
    sm = BACKPROP_VW(sm)
    sm$iteration = sm$iteration + 1
  }
  return(sm)
}

LIVE_PARAMETERS = function(training_result, weights_on_input_layer)
{
  live_betas = 1 * (training_result$bias != 0) + sum(abs(training_result$beta) > 1e-6)
  live_vw = 0
  K = length(training_result$w)
  if (weights_on_input_layer == TRUE)
  {
    for(k in 1 : K)
      live_vw = live_vw + sum(training_result$w[[k]] != 0) +
        sum(training_result$v[[k]] != 0)
  }
  else
  {
    for(k in 1 : K)
      live_vw = live_vw + sum(training_result$w[[k]] != 0)
  }
  cat("Live beta = ", live_betas, "Live vw = ", live_vw, "\n")
  return(live_betas + live_vw)
}


MODEL_SELECTION = function(train, validation, ms, fits)  
{
  test_fitted_values = list()
  number_fits = length(fits)
  save_fits = matrix(0, nrow = length(validation$response), ncol = (number_fits / 2))
  if (ms == "cv")
  {
    rss = rep(0, number_fits)
    for (k in 1 : number_fits)
    {
      test_fit = GANN_TEST(validation$integratedX, validation$response, fits[[k]])
      rss[k] = test_fit$loss_value
      test_fitted_values[[k]] = test_fit$yhat
      if (round(k/2) == 0.5 * k)
        save_fits[, (k/2)] = test_fit$yhat
    }
    return(list(opt_index = which.min(rss), rss = rss[which.min(rss)],
                test_fit = test_fitted_values[[which.min(rss)]], save_fits = save_fits))
  }
  if (ms == "ic")  #information criteria
  {
    sample_size = fits[[1]]$n
    penalty = log(sample_size)
    ic = rep(0, number_fits)
    for (k in 1 : number_fits)
    {
      rss_without_penalty = fits[[k]]$loss_value -
        fits[[k]]$lambda * sum(abs(fits[[k]]$beta))
      if (fits[[k]]$dimension == 1)
      {
        rss_without_penalty = 0.5 * mean((fits[[1]]$y - mean(fits[[1]]$y))^2)
        fits[[k]]$yhat = rep(mean(fits[[1]]$y), sample_size)
      }
      ic[k] = sample_size * log(rss_without_penalty) + penalty * fits[[k]]$dimension
      if (round(k/2) == 0.5 * k)
        save_fits[, (k/2)] = fits[[k]]$yhat
    }
    return(list(opt_index = which.min(ic), rss = ic[which.min(ic)],
                test_fit = test_fitted_values, save_fits = save_fits))
  }
  if (ms == "gcv")  
  {
    sample_size = fits[[1]]$sample_size
    gcv = rep(0, number_fits)
    for (k in 1 : number_fits)
    {
      dimension = fits[[k]]$dimension + 3 * (fits[[k]]$dimension - 1) / 2
      denumerator = sample_size * (1 - dimension / sample_size)^2
      gcv[k] = fits[[k]]$loss_value / denumerator
    }
    return(list(opt_index = which.min(gcv), rss = gcv[which.min(gcv)],
                test_fit = test_fitted_values))
  }
  
}

#------------------------------------------------------------------------------------
# Test the performance
GANN_TEST = function(newx, newy, sm)
{
  result = list()
  sm_test = sm
  sm_test$integratedX = newx
  sm_test$y = newy
  sm_test$n = length(newy)
  sm_test = FEEDFOWARD(sm_test, lambda = 0, sm$index_active)
  if(sm$err_fct == "ce")
  {
    result$yhat = 1 * (sm_test$yhat >= 0.5) 
    result$loss_value = sm_test$loss_value
  }
  if(sm$err_fct == "sse")
  {
    result$yhat = sm_test$yhat
    result$loss_value = sm_test$loss_value*2
  }
  return(result)
}

```
Then we run the main function.

```{r}
INTEGRATED_GANN = function(DList, y, number_of_function, out_fct,
                           err_fct, lrate, number_lambdas,
                           lambda_max, weights_on_input_layer, epsilon_lambda,
                           maxiter = 500, ms="cv", epsilon_iterations = 1e-5)
{
  integration_index = INTEGRATE_INDEX(DList)
  x = PREDICTOR(DList)
  ##Make setup function 
  sm = SETUPd(x, y, number_of_function, out_fct, err_fct,
              integration_index, lrate, weights_on_input_layer,
              maxiter, epsilon_iterations)
  ## Integration_index = list of indices that indicate which predictors are integrated
  integratedX = INTEGRATE_PREDICTORS(x, integration_index)  
  ##Split the data into train / validation
  data_split = DATA_SPLIT(y, integratedX, ms="cv", train_rate = 0.7)  
  sm$integratedX =data_split$train$integratedX
  sm$y=data_split$train$response
  sm$n = length(data_split$train$response)  
  ## Generate list of lambda
  lambda_list = LAMBDAS_GIVEN_LAMBDA_MAX(lambda_max, number_lambdas, epsilon_lambda)
  sm$lambda_list = lambda_list
  sm = FUNCTION_SELECTION(sm)
  sm = SYSTEMATIC_INITIAL_WEIGHT_GANN(sm)
  C_matrices = list()   
  ## Make empty list to save results for each lambda
  training_result = list()
  sm$beta_C=sm$beta_C
  for (lambdas in 1 : number_lambdas)
  {
    cat("===============================================", "\n")
    cat(lambdas, "th lambda", lambda_list[lambdas], "\n")
    training_result[[lambdas]] = TRAIN_GANN(sm, lambda_list[lambdas])
    training_result[[lambdas]]$dimension = LIVE_PARAMETERS(training_result[[lambdas]],
                                                           weights_on_input_layer)
    cat("Beta =", training_result[[lambdas]]$beta, "\n")
    cat("Number of live parameters = ", training_result[[lambdas]]$dimension, "\n")
    cat("Corresponding Training RSS =", training_result[[lambdas]]$loss_value, "\n")
    
    for(i in 1:sm$K)
    {
      sm$beta_C[[i]][lambdas] = training_result[[lambdas]]$beta[i]
    }
    
    ## Pass estimated coefficients to list sm, so that at next lambda, those becomes
    ## initial weights
    sm$bias = training_result[[lambdas]]$bias
    sm$beta = training_result[[lambdas]]$beta
    sm$v = training_result[[lambdas]]$v
    sm$w = training_result[[lambdas]]$w
    sm$index_active = training_result[[lambdas]]$index_active
    ## If we want to check the fit of each consolidated values, these five lines should be put
    beta_vector = as.vector(training_result[[lambdas]]$beta)
    if (training_result[[lambdas]]$dimension == 1)
      C_matrices[[lambdas]] = matrix(0, sm$n, sm$K)
    else if(length(beta_vector) == 1)
      C_matrices[[lambdas]] = training_result[[lambdas]]$C %*% matrix(beta_vector)
    else
      C_matrices[[lambdas]] = training_result[[lambdas]]$C %*% diag(beta_vector)
  }
  model_search = MODEL_SELECTION(data_split$train, data_split$validation, ms,
                                 training_result)
  optimal_fit = model_search$opt_index
  opt_model = training_result[[optimal_fit]]
  opt_model$testdata_fitted = model_search$test_fit
  opt_model$save_fits = model_search$save_fits
  opt_model$C_matrices = C_matrices
  opt_model$beta_C=sm$beta_C
  opt_model$data_y = data_split$validation$response
  cat("===============================================", "\n")
  cat("Model selection criterion =", ms, "\n")
  cat("Optimized at", optimal_fit, "th iteration", "\n")
  cat("Optimal lambda = ", opt_model$lambda, "\n")
  cat("Corresponding Test rss =", model_search$rss, "\n")
  cat("===============================================", "\n")
  return(opt_model)
}
```
We can find biomarkers (`selected.gene`) by following step.

```{r}
DList<-get(load("BComics.Rdata"))
y<-read.csv("y.csv")
y<-as.matrix(y)
out_fct="logistic" # out_fct: ReLU, linear, logistic
err_fct="ce"
lrate=0.001
weights_on_input_layer="TRUE" 
maxiter=900
lambda_max=1
number_lambdas=20
epsilon_lambda=0.001
number_of_function = rep(5,828)
ms="cv"
result=INTEGRATED_GANN(DList, y, number_of_function, out_fct,
                       err_fct,lrate, number_lambdas,
                       lambda_max, weights_on_input_layer, epsilon_lambda,
                       maxiter = 500, ms, epsilon_iterations = 1e-5)
save(result,file="bcresult.Rdata")
selected.gene<-gene_name[c(result$index_active)]
selected.gene
```
