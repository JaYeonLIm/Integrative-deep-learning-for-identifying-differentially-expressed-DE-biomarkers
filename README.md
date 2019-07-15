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

## Functions

`INTEGRATE_PREDICTORS(x, integration_index)` : makes 'integratedX'

`LAMBDAS_GIVEN_LAMBDA_MAX(lambda_max, number_lambdas,epsilon_lambda)` : makes lambda lists. Optimal lamda is selected in the function `MODEL_SELECTION`.

`TRAIN_GANN(sm, lambda)`: trains the models with training datasets. It contains the functions `FEEDFOWARD`,`BACKPROP_BETA` and `BACKPROP_VW`.

`MODEL_SELECTION(train, validation, ms, fits)` : First,`GANN_TEST(newx, newy, sm)`functions test the accuracy of trained models. Then, find the optimal model (opitmal lambda).


## Examples 

First, we need to copy the pacakge "DeepOmics2" in current dirrectory. 
Then, load the package. 
```{r}
library(DeepOmics2)
```


We can find biomarkers (`selected.gene`) by following step.

```{r}
# Reference Arguments 
DList<-get(load("BComics.Rdata"))   # load multiple gene expression datasets
y<-read.csv("y.csv")                # clinical data which consisits of 1 or 0
y<-as.matrix(y)            
out_fct="logistic"                  # out_fct: ReLU, linear, logistic
err_fct="ce"                        # "sse", "ce"
lrate=0.001                         # learning rate
weights_on_input_layer="TRUE" 
maxiter=900                         # maximum iterations
lambda_max=1
number_lambdas=20                   # the number of lamdas
epsilon_lambda=0.001                # threshold in training
number_of_function = rep(5,828)     # the number of activation function in each integrated X 
ms="cv"                             # 'cv' : cross validation

# Reference Functions
result=INTEGRATED_GANN(DList, y, number_of_function, out_fct,
                       err_fct,lrate, number_lambdas,
                       lambda_max, weights_on_input_layer, epsilon_lambda,
                       maxiter = 500, ms, epsilon_iterations = 1e-5)
save(result,file="bcresult.Rdata")
selected.gene<-gene_name[c(result$index_active)]
selected.gene
```
