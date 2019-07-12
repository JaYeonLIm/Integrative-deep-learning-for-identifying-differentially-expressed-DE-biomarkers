# Integrative-deep-learning-for-identifying-differentially-expressed-DE-biomarkers
## Description
This package is aim to find biomarkers via integrative deep learning. We propose adding an integrated layer to the deep learning structure, which would enable the effective analysis of genetic data as well as the discovery of significant biomarkers of disease. `INTEGRATED_GANN(DList, y, number_of_function, out_fct,
                       err_fct,lrate, number_lambdas,
                       lambda_max, weights_on_input_layer, epsilon_lambda,
                       maxiter = 500, ms, epsilon_iterations = 1e-5)` is main fucntion where we train datasets and test the performance of trained model.

## Arguments

`DList` : a list. It consists of multiple gene expression datasets. Each dataset contains same dimension. We use the function `PREDICTOR(DList)` to make the list to predictor of models. 

`y`: a matrix. It contains clinical data which consisits of 1 or 0 as the model aims at classification.




```{r cars}
summary(cars)
```

## Examples 

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
