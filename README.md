# hemiphoto2LAI (a beta version):  An R package for estimating leaf area index & leaf angle distribution from hemiphotos
----
## Description
Digital hemispherical photography (i.e., hemiphoto) is a convenient, rapid tool to estimate two key canopy structure parameters: leaf area index (LAI) and leaf angle distribution (LAD). Underlying such use of hemiphotos is a gap probability model widely known as either the Poisson model or Beer's Law for describing light-vegetation interactions. Based on this gap formula, numerous algorithms have been developed to convert hemiphotos into LAI and LAD. The hemiphoto2LAI package implements a total of 135 LAI estimation models, including the majority of classical algorithms proposed during the past several decades and more importantly, a newly proposed binary nonlinear regression algorithm (BNR). The implemented algorithms are adopted to a total of 19 common leaf angle distribution models. Details on the BNR algorithm and the package can be found in Zhao et al. (2019).
 
 ----
## Installation

There are three ways to install the **hemiphoto2LAI** R package.

### 1. CRAN

**hemiphoto2LAI** will be submitted to CRAN soon. Once it is available there, you can install it with:

```R
install.packages("hemiphoto2LAI")
```

### 2. GitHub

You can alternatively install the development version of **hemiphoto2LAI** from [GitHub](https://github.com/zhaokg/hemiphoto2LAI) as follows:

```R
if (!require(devtools)) install.packages('devtools')
devtools::install_github("zhaokg/hemiphoto2LAI")
```

Note that the avove will install "hemiphoto2LAI" from source. Becaues hemiphoto2LAI was written in the mixed use of C/C++ and Fortran. You need to make sure your machine is able to have a C and a Fotran compiler appropriately set up. For example, see [Package Development Prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) for the tools needed for your operating system. In particular, on Windows platforms, the most convenient option to go with the Rtools toolkit.

### 3. Pre-compiled binary

For those users or machines that cann't build the package from source, we will provide pre-compiled binary packages for you to install

#### Windows x86 or x64
```R
install.packages("https://github.com/zhaokg/hemiphoto2LAI/blob/master/precompiled_binary/hemiphoto2LAI_0.1.zip" ,repos=NULL)
```

#### Debian Linux x86 x64
```R
install.packages("https://github.com/zhaokg/hemiphoto2LAI/blob/master/precompiled_binary/hemiphoto2LAI_0.1_R_x86_64-pc-linux-gnu.tar.gz" ,repos=NULL)
```

#### Ubuntu Linux x86 x64
```R
install.packages("https://github.com/zhaokg/hemiphoto2LAI/blob/master/precompiled_binary/hemiphoto2LAI_0.1_R_x86_64-pc-linux-gnu_ubuntu.tar.gz" ,repos=NULL)
```
#### Mac
```R
install.packages("the link is coming" ,repos=NULL)
```
## Usage
