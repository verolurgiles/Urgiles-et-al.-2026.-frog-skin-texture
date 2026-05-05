<h3><p align="center">Supplementary Code for:</p></h3>
<h2><p align="center">Warts and All: Evolutionary and ecological patterns in frog skin texture</p></h2>

<I>Veronica L. Urgiles<sup>1,2</sup>,  Dylan K. Wainwright <sup>3</sup>, Anna E. Savage <sup>1</sup>, Eric A. Riddell <sup>4 </sup>, Molly C. Womack<sup>5</sup></I>

1. Department of Biological Sciences, University of Central Florida, 4110 Libra Drive, Orlando, FL, USA, 32816

2. Instituto Nacional de Biodiversidad del Ecuador, Pasaje Rumipamba 341 y Av. de los Shyris, Quito, Pichincha, ECU

3. Department of Biological Sciences, Purdue University, 915 Mitch Daniels Blvd, West Lafayette, Indiana, USA

4. Department of Biology, University of North Carolina – Chapel Hill, Chapel Hill, NC 27599, USA

5. Ecology and Evolutionary Biology, Cornell University, 215 Tower Rd, Ithaca, NY, USA, 14853

------

**This repository contains the R code used in the study to promote transparency and reproducibility. Raw data and full model outputs are provided in the supplemental materials and in Dryad**. 

[TOC]

### **1. Raw Data Summary**

#### 	1.1. Species means dataset summary: 

| Variable           | **Type**  | **Description**                                              | **Min** | **1st Quartile** | **Median** | **Mean** | **3rd Quartile** | **Max** |
| ------------------ | --------- | ------------------------------------------------------------ | ------- | ---------------- | ---------- | -------- | ---------------- | ------- |
| **phylo_names**    | Character | Species names match to taxonomic ID's of the phylogeny of Portik et al., (2023) | —       | —                | —          | —        | —                | —       |
| **body_region**    | Factor    | Body region category                                         | —       | —                | —          | —        | —                | —       |
| **Ecology**        | Factor    | Microhabitat type                                            | —       | —                | —          | —        | —                | —       |
| **mean_roughness** | Numerical | Mean skin texture roughness                                  | 0.26    | 1.80             | 3.01       | 3.58     | 4.79             | 19.56   |
| **mean_peak**      | Numerical | Mean maximum peak height                                     | 1.98    | 9.26             | 13.56      | 16.38    | 20.9             | 62.42   |
| **sd_roughness**   | Numerical | Standard deviation of roughness                              | 0.0069  | 0.57             | 0.98       | 1.29     | 1.63             | 12.2    |
| **sd_peak**        | Numerical | Standard deviation of maximum peak height                    | 0.15    | 2.8              | 4.3        | 5.8      | 7.1              | 34.7    |

#### 

### **2. R packages** 

All analyses were performed in R version 4.4.2.

| Package   | Reference                   |
| --------- | --------------------------- |
| Tidyverse | Wickham et al. (2019)       |
| Rgdal     | Bivand et al. (2017)        |
| ape       | Paradis and  Schliep (2019) |
| phytools  | Revell (2024)               |
| mvMORPH   | Clavel et al., (2015)       |
| Raster    | Hijmans et al. (2020)       |
| RRPP      | Collyer and Adams (2025)    |
| lmerTest  | Kuznetsova et al. (2017)    |



### **3. Data preparation** 

```R

														###Data_preparation######

dat <- read.csv("Supplemental_Table1.csv")

#clean and transform data to log scale 
dat <- dat %>%
  mutate(
    log_Sp = ifelse(!is.na(mean_peak) & mean_peak > 0, log(mean_peak), NA_real_),
    log_Sq = ifelse(!is.na(mean_roughness) & mean_roughness > 0, log(mean_roughness), NA_real_)
  ) %>%
  mutate(
    across(c(mean_SVL),
           ~ as.numeric(scale(.)), .names = "{.col}_z")
  )
#remove species with missing microhabitat information
dat_clean <- dat %>%
  mutate(
    body_region = toupper(str_squish(body_region)),
    ecology     = str_to_lower(str_squish(as.character(ecology))),
    ecology     = na_if(ecology, ""),                 
    ecology     = na_if(ecology, "na")                
  )
#separate data by metric
comp <- dat_clean %>%
  group_by(phylo_names) %>%
  summarise(
    n_Sp    = n_distinct(body_region[!is.na(mean_peak)]),
    n_Sq    = n_distinct(body_region[!is.na(mean_roughness)]),
    n_both  = n_distinct(body_region[!is.na(mean_peak) & !is.na(mean_roughness)]),
    .groups = "drop"
  )
#include species with complete data for all body regions

species_sp4    <- comp %>% filter(n_Sp   == 4) %>% pull(phylo_names)
species_sq4    <- comp %>% filter(n_Sq   == 4) %>% pull(phylo_names)
species_both4  <- comp %>% filter(n_Sp   == 4, n_Sq == 4) %>% pull(phylo_names)

# Metric-specific datasets 
dat_sp_complete <- dat %>%
  filter(phylo_names %in% species_sp4, !is.na(mean_peak))

dat_sq_complete <- dat %>%
  filter(phylo_names %in% species_sq4, !is.na(mean_roughness))

#phylogenetic tree to prunn from Portik et al. (2023)

#portik_tree <- read.tree("Supplementary_File_S3.tre")

# names for matching 
tips_raw  <- portik_tree$tip.label
tips_norm <- tips_raw %>%
  str_squish() %>% str_replace_all(" +", "_") %>% str_replace_all("_+", "_") #may have white spaces

sp_sp_norm <- species_sp4 %>%
  str_squish() %>% str_replace_all(" +", "_") %>% str_replace_all("_+", "_")#may have white spaces

sp_sq_norm <- species_sq4 %>%
  str_squish() %>% str_replace_all(" +", "_") %>% str_replace_all("_+", "_")#may have white spaces

# Keep lists and prune tree
keep_sp_norm <- intersect(sp_sp_norm, tips_norm)
keep_sp_raw  <- tips_raw[match(keep_sp_norm, tips_norm)]
tree_sp_pruned <- keep.tip(portik_tree, keep_sp_raw)

keep_sq_norm <- intersect(sp_sq_norm, tips_norm)
keep_sq_raw  <- tips_raw[match(keep_sq_norm, tips_norm)]
tree_sq_pruned <- keep.tip(portik_tree, keep_sq_raw)

#keep relevant columns and convert data to wide format
sq_log_long <- dat_sq_complete %>%
  dplyr::select(phylo_names, body_region, ecology, log_Sq) %>%
  filter(!is.na(log_Sq))

sp_log_long <- dat_sq_complete %>%
  dplyr::select(phylo_names, body_region, log_Sp) %>%
  filter(!is.na(log_Sp))

# Wide format: roughness log_Sq (DA/DP/VA/VP)
X_logSq <- sq_log_long %>%
  group_by(phylo_names, body_region) %>%
  summarise(log_Sq = mean(log_Sq, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = body_region, values_from = log_Sq) %>%
  as.data.frame()
rownames(X_logSq) <- X_logSq$phylo_names; X_logSq$phylo_names <- NULL
X_logSq <- X_logSq[, intersect(c("DA","DP","VA","VP"), colnames(X_logSq)), drop = FALSE]

# Wide format: maximum peak height (DA/DP/VA/VP)
X_logSp <- sp_log_long %>%
  group_by(phylo_names, body_region) %>%
  summarise(log_Sp = mean(log_Sp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = body_region, values_from = log_Sp) %>%
  as.data.frame()
rownames(X_logSp) <- X_logSp$phylo_names; X_logSp$phylo_names <- NULL
X_logSp <- X_logSp[, intersect(c("DA","DP","VA","VP"), colnames(X_logSp)), drop = FALSE]
          
 
```



### **3. ** **Skin texture across body regions**

```R
##Note: Code examples below use the roughness dataset for illustration.
##Analyses for maximum peak heightfollow the same workflow and differ only in the input dataset.
#This specific analysis was performed using our individual level dataset
ind<- read.csv("Supplemental_Table2.csv")
ind$species <- ind$phylo_names
specimen.data <- ind>% 
            dplyr::select(museum_code, species, region, parameter, mean) %>% 
            dplyr::filter(parameter=="Sq") %>% 
            na.omit()

#reorder body region factor

specimen.data$region <- factor(specimen.data$region,
                     levels = c("VP", "VA", "DP", "DA"))

#linear mixed model
fm.region <- lmer(log(mean) ~ region  + (1 | species) +(1 | museum_code), data=specimen.data)
anova(fm.region)

#phylogenetic signal
lambda_values <- sapply(1:ncol(X_logSq), function(i) {
  phylosig(tree_u, X_logSq[, i], method = "lambda")$lambda
})


```

### **4. Rates of skin texture evolution** **across body regions**

```R

#####Evolutionary rates roughness########
####Function of Adams (2013)######
#Because body-region traits are inherently correlated the trait covariance matrix had very small eigenvalues. We therefore added a minimal ridge term to the covariance matrix to prevent issues arising from those very small eigenvalues.
##Note: Code examples below use the roughness metric
##Analyses for maximum peak height (X_logSp) follow the same workflow and differ only in the input dataset.
         
CompareRates.multTrait <- function(phy, x, TraitCov = TRUE, ms.err = NULL, ms.cov = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Need Matrix")
  
  x <- as.matrix(x)
  N <- nrow(x)
  p <- ncol(x)
  
  ##  VCV and alignment 
  C <- ape::vcv.phylo(phy)
  C <- C[rownames(x), rownames(x)]
  C <- as.matrix(Matrix::forceSymmetric(C))
  
  evC <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  if (min(evC) <= 0) {
    C <- C + diag(abs(min(evC)) + 1e-10, nrow(C))
  }
  
  ##  Optional measurement-error
  if (is.matrix(ms.err)) ms.err <- as.matrix(ms.err[rownames(x), , drop = FALSE])
  if (is.matrix(ms.cov)) ms.cov <- as.matrix(ms.cov[rownames(x), , drop = FALSE])
  
  ## Cholesky builder used in constrained (TraitCov=TRUE) case 
  build.chol <- function(b) {
    c.mat <- matrix(0, nrow = p, ncol = p)
    c.mat[lower.tri(c.mat)] <- b[-1]
    c.mat[p, p] <- exp(b[1])
    c.mat[1, 1] <- sqrt(sum((c.mat[p, ])^2))
    if (p > 2) {
      for (i in 2:(p - 1)) {
        val <- c.mat[1, 1]^2 - sum((c.mat[i, ])^2)
        c.mat[i, i] <- ifelse(val > 0, sqrt(val), 0)
      }
    }
    c.mat
  }
  
  ##  inverse + logdet via Cholesky
  safe_chol_inv_logdet <- function(M, jitter0 = 1e-10, jitter_max = 1e-4) {
    M <- as.matrix(Matrix::forceSymmetric(M))
   
    for (k in 0:8) {
      jj <- min(jitter_max, jitter0 * 10^k)
      Mj <- M + diag(jj, nrow(M))
      ch <- try(chol(Mj), silent = TRUE)
      if (!inherits(ch, "try-error")) {
        invM   <- chol2inv(ch)
        logdet <- 2 * sum(log(diag(ch)))
        return(list(inv = invM, logdet = logdet, jitter = jj))
      }
    }
    
   
    ev <- eigen(M, symmetric = TRUE)
    vals <- ev$values
    vecs <- ev$vectors
    
    tol <- max(dim(M)) * max(abs(vals)) * .Machine$double.eps
    vals_inv <- ifelse(abs(vals) > tol, 1/vals, 0)
    
    invM <- vecs %*% (vals_inv * t(vecs))
    logdet <- sum(log(pmax(abs(vals), .Machine$double.eps)))
    
    list(inv = invM, logdet = logdet, jitter = NA_real_)
  }
  
  ##  GLS means and observed rate matrix R.obs 
  Sinv  <- solve(C)
  a.obs <- colSums(Sinv) %*% x / sum(Sinv)
  
  D <- matrix(0, N * p, p)
  for (i in 1:(N * p)) {
    for (j in 1:p) {
      if ((j - 1) * N < i && i <= j * N) D[i, j] <- 1.0
    }
  }
  
  y   <- as.matrix(as.vector(x))
  one <- matrix(1, N, 1)
  
  R.obs <- t(x - one %*% a.obs) %*% Sinv %*% (x - one %*% a.obs) / N
  R.obs <- as.matrix(Matrix::forceSymmetric(R.obs))
  
  # if TraitCov == FALSE, enforce diagonal R.obs
  if (TraitCov == FALSE) R.obs <- diag(diag(R.obs), p)
  
  # make R.obs numerically PD
  evR <- eigen(R.obs, symmetric = TRUE, only.values = TRUE)$values
  if (min(evR) <= 0) {
    R.obs <- R.obs + diag(abs(min(evR)) + 1e-10, p)
  }
  
  ## 6. Build within-species ME block m.e if provided 
  m.e <- NULL
  
  if (is.matrix(ms.err) && !is.null(ms.err)) {
    m.e <- diag(as.vector(ms.err))
  }
  
  if (is.matrix(ms.err) && is.matrix(ms.cov) && !is.null(ms.cov)) {
    within.spp <- cbind(ms.err, ms.cov)
    rc.label <- NULL
    for (i in 1:p) rc.label <- rbind(rc.label, c(i, i))
    for (i in 1:p) {
      for (j in 2:p) {
        if (i != j && i < j) rc.label <- rbind(rc.label, c(i, j))
      }
    }
    m.e <- NULL
    for (i in 1:p) {
      tmp <- NULL
      for (j in 1:p) {
        for (k in 1:nrow(rc.label)) {
          if (setequal(c(i, j), rc.label[k, ]) == TRUE) {
            tmp <- cbind(tmp, diag(within.spp[, k]))
          }
        }
      }
      m.e <- rbind(m.e, tmp)
    }
  }
  
  ##  Observed (unconstrained) log-likelihood 
  LLik.obs <- {
    K <- kronecker(R.obs, C)
    if (!is.null(m.e)) K <- K + m.e
    IL <- safe_chol_inv_logdet(K)
    quad <- crossprod(y - D %*% t(a.obs), IL$inv %*% (y - D %*% t(a.obs)))
    as.numeric(-0.5 * (quad + N * p * log(2 * pi) + IL$logdet))
  }
  
  ## Reasonable start for diagonal value 
  sigma.mn <- mean(diag(R.obs))
  
  ##  Likelihood under no trait covariation (TraitCov = FALSE) 
  lik.covF <- function(sigma) {
    R <- R.obs
    diag(R) <- sigma
    K <- kronecker(R, C)
    if (!is.null(m.e)) K <- K + m.e
    IL <- safe_chol_inv_logdet(K)
    quad <- crossprod(y - D %*% t(a.obs), IL$inv %*% (y - D %*% t(a.obs)))
    LLik <- as.numeric(-0.5 * (quad + N * p * log(2 * pi) + IL$logdet))
    if (!is.finite(LLik)) LLik <- -1e+10
    return(-LLik)  # minimizer
  }
  
  ## Likelihood for no trait covariation (TraitCov = FALSE)
  lik.covT <- function(sigma) {
    low.chol <- build.chol(sigma)
    R <- low.chol %*% t(low.chol)
    R <- as.matrix(Matrix::forceSymmetric(R))
    K <- kronecker(R, C)
    if (!is.null(m.e)) K <- K + m.e
    IL <- safe_chol_inv_logdet(K)
    quad <- crossprod(y - D %*% t(a.obs), IL$inv %*% (y - D %*% t(a.obs)))
    LLik <- as.numeric(-0.5 * (quad + N * p * log(2 * pi) + IL$logdet))
    if (!is.finite(LLik)) LLik <- -1e+10
    return(-LLik)  # minimizer
  }
  
  ## Optimize 
  if (TraitCov == FALSE) {
    model1 <- optim(sigma.mn, fn = lik.covF, method = "L-BFGS-B",
                    hessian = TRUE, lower = 0)
    R.constr <- diag(model1$par, p)
  } else {
    R.offd <- rep(0, (p * (p - 1) / 2))
    model1 <- optim(par = c(sigma.mn, R.offd), fn = lik.covT,
                    method = "L-BFGS-B")
    chol.mat <- build.chol(model1$par)
    R.constr <- chol.mat %*% t(chol.mat)
  }
  
  optimmessage <- if (model1$convergence == 0) {
    "Optimization has converged."
  } else {
    "Optim may not have converged. Consider changing start value or limits."
  }
  
  ##  LRT (df = p - 1) 
  LRT      <- -2 * ((-model1$value - LLik.obs))
  LRT.prob <- pchisq(LRT, (p - 1), lower.tail = FALSE)
  
  AIC.obs    <- -2 * LLik.obs        + 2 * p + 2 * p
  AIC.common <- -2 * (-model1$value) + 2 + 2 * p
 
  list(
    Robs             = R.obs,
    Rconstrained     = R.constr,
    Lobs             = LLik.obs,
    Lconstrained     = (-model1$value),
    LRTest           = LRT,
    Prob             = LRT.prob,
    AICc.obs         = AIC.obs,
    AICc.constrained = AIC.common,
    optimmessage     = optimmessage
  )
}


stopifnot(all(c("DA","DP","VA","VP") %in% colnames(X_logSq)))
X_sq <- as.matrix(X_logSq[, c("DA","DP","VA","VP")])
X_sq <- X_sq[tree_u$tip.label, , drop = FALSE]   # reorder to tree
stopifnot(identical(rownames(X_sq), tree_u$tip.label))
stopifnot(!anyNA(X_sq))

res_sq_full <- CompareRates.multTrait(tree_u, X_sq, TraitCov = TRUE)

#####Evolutionary rates maximum peak height######## 
X_logSp_aligned

stopifnot(all(c("DA","DP","VA","VP") %in% colnames(X_logSp_aligned)))
X_sp <- as.matrix(X_logSp_aligned[, c("DA","DP","VA","VP")])
X_sp <- X_sp[tree_sp_u$tip.label, , drop = FALSE]   # reorder to tree
stopifnot(identical(rownames(X_sp), tree_sp_u$tip.label))
stopifnot(!anyNA(X_sq))

res_Sq_cov <- CompareRates.multTrait(tree_u, X_sp, TraitCov = TRUE)
```



### **5. Models of evolution** 

```R
####MVMOPRH EVOLUTIONARY MODELS########
##Note: Code examples below use the roughness metric
##Analyses for maximum peak height (X_logSp) follow the same workflow and differ only in the input dataset.

X_logSq <- X_logSq[tree_sq_pruned$tip.label, , drop = FALSE]
stopifnot(identical(rownames(X_logSq), tree_sq_pruned$tip.label))
tree_u <- try(force.ultrametric(tree_sq_pruned, method = "nnls"), silent = TRUE)
if (inherits(tree_u, "try-error") || !is.ultrametric(tree_u)) {
  tree_u <- force.ultrametric(tree1, method = "extend")
}
stopifnot(is.ultrametric(tree_u))

###regimes mapping and matrix ####
reg_per_sp <- sq_log_long %>%
  select(phylo_names, ecology) %>%
  mutate(
    phylo_names = str_squish(phylo_names),
    regime = ecology %>%
      str_squish() %>% str_to_lower() %>% str_replace_all("[- ]", "_")
  ) %>%
  filter(!is.na(regime), regime != "") %>%
  distinct(phylo_names, regime) %>%
  group_by(phylo_names) %>%
  summarise(regime = first(regime), .groups = "drop")
 ###align to tree tips 
reg_aligned <- tibble(phylo_names = tree_u$tip.label) %>%
  left_join(reg_per_sp, by = "phylo_names")
  
#### vector factor for mvMORPH ###
regimes <- setNames(factor(reg_aligned$regime), reg_aligned$phylo_names)
regimes <- droplevels(regimes)
lev <- levels(regimes); k <- length(lev)
  
set.seed(47)
K <- 50
pi_root <- setNames(rep(1/k, k), lev) 
maps <- lapply(seq_len(K), function(i)
  phytools::make.simmap(tree, regimes, model = "ER", nsim = 1,
                        pi = pi_root, message = FALSE))
class(maps) <- "multiSimmap"  

 ###running multivariate mvmorph models###
fit_bm1_sq <- mvBM(tree_u, X_logSq, model = "BM1", scale.height = TRUE)
fit_ou1_sq <- mvOU(tree_u, X_logSq, model = "OU1", scale.height = TRUE)

lev <- levels(regimes); k <- length(lev)
pi_root <- setNames(rep(1/k, k), lev)
               
########running 100 simmap#######            
maps <- phytools::make.simmap(tree_u, regimes, model="ER", nsim=100,
                              pi=pi_root, message=FALSE)

maps_list <- if (inherits(maps, "multiSimmap")) maps else list(maps)

maps_ok <- Filter(function(m){
  me <- m$mapped.edge
  all(lev %in% colnames(me)) && all(colSums(me)[lev] > 0)
}, maps_list)
if (!length(maps_ok)) stop("error")

####making OUM for each map #####
fits_oum <- lapply(seq_along(maps_ok), function(i){
  cat(sprintf("Fitting OUM on map %d/%d...\n", i, length(maps_ok))) 
  try(mvOU(maps_ok[[i]], X_logSq, model = "OUM", scale.height = TRUE,
           control = list(maxit = 1000)), silent = TRUE)
})
oum_ok <- Filter(function(f) !is.null(f) && !inherits(f, "try-error"), fits_oum)
if (!length(oum_ok)) stop("error")
                 
n <- nrow(X_logSq)

AICc_from_fit <- function(f) {
  ll <- stats::logLik(f)                         
  k  <- as.numeric(attr(ll, "df"))              
  llv <- as.numeric(ll)
  aic  <- -2*llv + 2*k
  aicc <- aic + (2*k*(k+1)) / max(n - k - 1, 1e-8)
  c(AIC = aic, AICc = aicc, k = k, logLik = llv)
}

bm1 <- AICc_from_fit(fit_bm1_sq)
ou1 <- AICc_from_fit(fit_ou1_sq)
oum <- t(sapply(oum_ok, AICc_from_fit))
                 
                 
oum_df <- as.data.frame(oum, stringsAsFactors = FALSE)
rownames(oum_df) <- NULL
                 
ll_bm1 <- as.numeric(bm1)                 
ll_ou1 <- as.numeric(ou1)                
ll_oum <- as.numeric(oum)
                 
 ## counts
p <- ncol(X_logSq)                        
r <- nlevels(regimes)                     
n <- nrow(X_logSq)            

k_bm1 <- p + p*(p+1)/2
k_ou1 <- p + p*(p+1)
k_oum <- r*p + p*(p+1)                    

aic  <- function(ll, k) -2*ll + 2*k
aicc <- function(ll, k, n) aic(ll,k) + (2*k*(k+1)) / max(n - k - 1, 1e-8)

AICc_bm1 <- aicc(ll_bm1, k_bm1, n)
AICc_ou1 <- aicc(ll_ou1, k_ou1, n)
AICc_oum <- aicc(ll_oum, k_oum, n)   

comp <- data.frame(
  model = c("BM1","OU1", paste0("OUM_map", seq_along(ll_oum))),
  logLik = c(ll_bm1, ll_ou1, ll_oum),
  k      = c(k_bm1,  k_ou1,  rep(k_oum, length(ll_oum))),
  AICc   = c(AICc_bm1, AICc_ou1, AICc_oum),
  row.names = NULL
)


                 
  
                             
```

### **6. Correlations between skin texture variation and climate and body size**

#### 6.1 Summary of individual-level analysis of arboreal, leaflitter and burrower frogs

| Variable        | **Type**  | **Description**                                              | **Min** | **1st Quartile** | **Median** | **Mean** | **3rd Quartile** | **Max** |
| --------------- | --------- | ------------------------------------------------------------ | ------- | ---------------- | ---------- | -------- | ---------------- | ------- |
| **phylo_names** | Character | Species names match to taxonomic ID's of the phylogeny of Portik et al., (2023) | —       | —                | —          | —        | —                | —       |
| **Museum_code** | Factor    | Body region category                                         | —       | —                | —          | —        | —                | —       |
| **Region**      | Factor    |                                                              | —       | —                | —          | —        | —                | —       |
| **Ecology**     | Factor    | Microhabitat type                                            | —       | —                | —          | —        | —                | —       |
| **Parameter**   | Numerical | Sq= roughness and Sp= maximum peak height                    | —       | —                | —          | —        | —                | —       |
| **log_value**   | Numerical | Log transformed values of roughness and maximum peak height () | 0.327   | 2.94             | 6.6        | 10       | 13.27            | 65      |
| **SVL**_**z**   | Numerical | Body size (mm) (scaled)                                      | 13.3    | 26.8             | 37.6       | 43.3     | 56.18            | 178.8   |
| **Aridity_z**   | Numerical | Aridity index (scaled)                                       | 0.005   | 0.66             | 1.16       | 1.29     | 1.82             | 4.73    |
| **evapo_z**     | Numerical | Potential evapotranspiration (scaled)                        | 0.007   | 0.13             | 0.14       | 0.15     | 0.17             | 0.27    |

#### 6.2 Individual-level phylogenetic generalized least-squares model (E-pgls)

```R
##Note: Code examples below use the roughness metric
##Analyses for maximum peak height follow the same workflow and differ only in the input dataset.

ind <- read.csv("Supplemental_Table_1.csv")


ind.portik <- ind.portik %>%
  dplyr::select(
    phylo_names, parameter, mean, SVL,
    aridity, evapo, region, museum_code, ecology
  ) %>%
  dplyr::filter(parameter %in% c("Sq")) %>%
  dplyr::filter(
    !is.na(mean),
    !is.na(SVL),
    !is.na(evapo),
    !is.na(aridity),
    !is.na(phylo_names), 
    !is.na(ecology), 
  )

unique(is.na(ind.portik))


life <- c("terrestrial_leaflitter", "terrestrial_burrower", "arboreal")

######life######

ind.life <- ind.portik %>%
  filter(ecology %in% life)

##### scale variables ###
ind.d.f <- ind.life %>%
  mutate(log_value = log(mean),
         SVL_z = scale(SVL),
         aridity_z = scale(aridity),
         evapo_z = scale(evapo)
  )

######crop phylogeny####

my.species <- unique(ind.d.f$phylo_names)

as.data.frame(my.species)

tree.crop<-drop.tip(portik,portik$tip.label[which(portik$tip.label%in%my.species==FALSE)])
tree.crop$tip.label


ind.portik <- ind.d.f[order(match(ind.d.f$phylo_names, tree.crop)),]


fit1 <- lm.rrpp.ws(log_value ~ (region*SVL_z) + (region*aridity_z) +
                     (region*evapo_z), data = ind.portik,
                   subjects = "phylo_names",
                   Cov = vcv(tree.crop),
                   gamma = "sample",
                   print.progress = FALSE)

pairwise_test <- pairwise(fit1, fit.null = NULL, groups = ind.portik$region, covariate = NULL)

summary(pairwise_test)
```

### **7. References**

[Adams DC. Comparing evolutionary rates for different phenotypic traits on a phylogeny using likelihood. Syst Biol. 2013 Mar;62(2):181–92.](http://paperpile.com/b/3M0dzJ/N9UQ)

[Bivand R, Keitt T, Rowlingson B, Pebesma E, Sumner M, Hijmans R. Rgdal: Bindings for the geospatial data abstraction library 2017. URL https://CRAN R-project org/package= rgdal R package version [Internet\]. 2018 [cited 2024 Dec 6];1](http://paperpile.com/b/3M0dzJ/rOQv)

[Clavel J, Escarguel G, Merceron G. Mvmorph: An r package for fitting multivariate evolutionary models to morphometric data. Methods Ecol Evol. 2015 Nov 1;6(11):1311–9](http://paperpile.com/b/3M0dzJ/HJh3)

[Collyer MC, Adams DC. Linear Model Evaluation with Randomized Residuals in a Permutation Procedure “RRPP”. CRAN; 2025. Available from: ](http://paperpile.com/b/3M0dzJ/WTK2)https://cran.r-project.org/web/packages/RRPP/RRPP.pdf

[Hijmans RJ, van Etten J. Raster: geographic data analysis and modeling. R package. 2020;](http://paperpile.com/b/3M0dzJ/EQTB)

[Kuznetsova A, Brockhoff PB, Christensen RHB. LmerTest package: Tests in linear mixed effects models. J Stat Softw. 2017 Dec 6;82(13):1–26.](http://paperpile.com/b/3M0dzJ/FkAf)

Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.” *Bioinformatics*, **35**, 526-528. [doi:10.1093/bioinformatics/bty633](https://doi.org/10.1093/bioinformatics/bty633).

[Portik DM, Streicher JW, Wiens JJ. Frog phylogeny: A time-calibrated, species-level tree based on hundreds of loci and 5,242 species. Mol Phylogenet Evol. 2023 Nov;188(107907):107907.](http://paperpile.com/b/3M0dzJ/IE4M)

[Revell LJ. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ. 2024 Jan 5;12(e16505):e16505.](http://paperpile.com/b/3M0dzJ/TSfi)

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” *Journal of Open Source Software*, **4**(43), 1686. [doi:10.21105/joss.01686](https://doi.org/10.21105/joss.01686).
