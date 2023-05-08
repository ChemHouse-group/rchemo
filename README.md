## rchemo - Dimension reduction, Regression and Discrimination for Chemometrics  

**rchemo** is a [package](https://github.com/ChemHouse-group/rchemo/blob/main/inst/rchemo_functions_github.md) for **data exploration and prediction** with focus on **high dimensional data** and **chemometrics**. 

The package was initially designed about **partial least squares regression and discrimination models** and variants, in particular locally weighted PLS models (**LWPLS**) (e.g. https://doi.org/10.1002/cem.3209).
Then, it has been expanded to many other methods for 
analyzing high dimensional data. 

The name **rchemo** comes from the fact that the package is orientated to chemometrics, but most of the provided methods are fully **generic to other domains**. 

Functions such as **transform**, **predict**, **coef** and **summary** are available. 
**Tuning the predictive models** is facilitated by generic functions **gridscore** (validation dataset) and 
**gridcv** (cross-validation). Faster versions are also available for models based on latent variables (LVs) 
(**gridscorelv** and **gridcvlv**) and ridge regularization (**gridscorelb** and **gridcvlb**).

All the functions have a **help page** with a documented example. 

**NOTE**: This repository replaces the previous [rchemo repository](https://github.com/mlesnoff/rchemo) that now is archived. 

## <span style="color:green"> **News** </span> 

Click [**HERE**](https://github.com/ChemHouse-group/rchemo/blob/main/inst/NEWS.md) to see **what changed** in the previous versions. 

or write in the R console
```{r}
news(package = "rchemo")
```

## <span style="color:green"> **Installation** </span> 

Using [**Rstudio**](https://posit.co/download/rstudio-desktop/) is recommended for installation and usage.

**rchemo** can be installed from **CRAN**, or from Github using the following steps. 

#### <span style="color:green"> 1.  Install package **'remotes'** from **CRAN** </span>

Use the **Rstudio** menu 

or write in the R console
```{r}
install.packages("remotes")
```

#### <span style="color:green"> 2. Install package **'rchemo'** </span> 

**a) Most recent version**

Write in the R console
```{r}
remotes::install_github("ChemHouse-group/rchemo", dependencies = TRUE)
```
In case of the following question during installation process:
```{r}
These packages have more recent versions available.
Which would you like to update?"
```
it is recommended to skip updates (usually choice **3** = None)

**b) Any given tagged version**

e.g. with tag "v0.1-0", write in the R console
```{r}
remotes::install_github("ChemHouse-group/rchemo@v0.1-0", dependencies = TRUE)
```

## <span style="color:green"> **Usage** </span> 

Write in the R console
```{r}
library(rchemo)
```

## <span style="color:green"> **How to cite** </span> 

R package rchemo: Dimension Reduction, Regression and Discrimination for Chemometrics. https://github.com/ChemHouse-group/rchemo.




