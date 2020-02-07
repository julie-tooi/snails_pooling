# In with code we use librarydown, please download them for running this code
library(plyr)
library(dplyr)
library(limma)
library(ggplot2)
library(tidyr)
library(dbframe) # ## remotes::install_github("grayclhn/dbframe-R-library") (first you need download remotes) or another variant for download libraries frim Github https://rdrr.io/github/grayclhn/dbframe-R-library/f/README.md 

#also you should set a number in set.seed() for correct randomisation
set.seed()

# Function for creation metadata for your samples
# n_g - vector with amount of biological samples in each group
# n_tech - vector with amount of technical repeats of each sample
gen_groups <- function(n_g, n_tech) {
  # n_biol - overall amount of biological samples
  n_biol <- sum(n_g)
  if (length(n_tech) != n_biol) stop("'length(n_tech)' must be == 'n_biol'")
  # Factor specifying groups of samples (tritment)
  G <- unlist(mapply(rep, x = paste0("G", seq_len(length(n_g))), each = n_g))
  G <- unlist(mapply(rep, x = G, each = n_tech))
  # Biological repeats
  B <- unlist(mapply(rep, x = paste0("B", seq_len(n_biol)), each = n_tech))
  R <- rep(1, length(G))
  res <- data.frame(G, B, R)
  return(res)
}
#Example. 120 samples, 60 in each group
meta_dat <- gen_groups(n_g = c(60, 60), n_tech = rep(1, 120))

# Function which add to metadata, data about samples poolingД
# dat - metadata
# n_pools - vector with amount of pools - should be a natural number and number of samples in each group must be divisible by it
add_pools <- function(dat, n_pools){
  tol <- .Machine$double.eps^0.5
  if(!all(n_pools > tol & abs(n_pools - round(n_pools)) < tol)) stop("n_pools must be a vector of even natural numbers")
  if(any(n_pools %% 2 != 0)) stop("n_pools must be a vector of even numbers")
  if(any(nrow(dat) %% n_pools != 0)) stop("nrow(dat) must be divisible by n_pools")
  n_rep <- nrow(dat) / n_pools
  pools <- mapply(function(a, b) rep(paste0("P", sprintf("%03.0f", 1:a)), each = b), a = n_pools, b = n_rep)
  dimnames(pools)[[2]] <- paste0("Pools_", n_pools)
  data.frame(dat, pools)
}
#WARNING! - In futher function metadata with pools ALWAYS named meta_pools, because in some function we should make it edfault value of arguments

#Example 120 Biological samples, and we pools them - 6 pools (20 samples in each), 8 pools(15 samples in each), 10 pools (12 samples in each), 12 pools (10 samples in each), 20 pools (6 samples in each), 30 pools (4 samples in each)
meta_pools <- add_pools(meta_dat, number_of_pools)

#Function which create data about each biological sample
#beta_0 - level of expression in group 1 (base level)
#beta_1 - how expression in group 2 is differ from group 1
#sig_B - level of biological dispersion
#common_meta_dat - meta_data (we use metadata with information about pools (made by add_pools function))
#return vector with value for each biological sample
data_creation <- function(beta_0, beta_1, sig_B, common_meta_dat){
  b_i <- rnorm(nlevels(common_meta_dat$B), mean = 0, sd = sig_B)
  Y_start <- beta_0 + (as.numeric(common_meta_dat$G)-1)*beta_1+ b_i
  return(Y_start)
}

#Example: 3 - level in group, in average level in group 2 more than in 1st,  level of biological dispersion - 1, meta_data are in meta_pools
#test_exp1 <- data_creation(beta_0 = 3, beta_1 = 2,  sig_B = 1, common_meta_dat = meta_pools)

#Function. Repeat data_creation n_strings times
#Return matrix, where each column is data about sample, and each string is value of data we simulated
data_collecting <- function(n_strings, beta_0, beta_1, sig_B, common_meta_dat){
  Y_results <- RepParallel(n_strings, data_creation(beta_0, beta_1, sig_B, common_meta_dat))
  return(Y_results)
}

#Example: 100 - number of repeats, 3 - level in group, in average level in group 2 more than in 1st,  level of biological dispersion - 1, meta_data are in meta_pools
#test_exp2 <- data_collecting(100, 3, 2, 1, meta_pools) 

#Function which make pools - mean value for biological samples in pool. Taking into account technical repetition
#Y_result - vectors of values for samples we pool
#index - vector same lenght then Y_result, how we pooled samples (for each sample number of pool)
#sig_e - technical error (variability determed technical errors)
#return vector with values for each pool
pools <- function(Y_result, index, sig_e){
  Y_pull <- rep(as.vector(tapply(Y_result, INDEX = index, FUN = mean, simplify = TRUE)), each = 2)
  Y_pull <- Y_pull + rnorm(length(Y_pull), mean = 0, sd = sig_e)
  return (as.array(Y_pull))
}

#Example. test_exp1 - vector of values, meta_pools$Pools_10 - index, sig_e = 1)
#test_exp3 <- pools(test_exp1, meta_pools$Pools_10, 1)

#Function. Make all variants of pools from metadata with pools
#Y_res - matrix (from data_collecting)
#Return List containing list with matrix of values for each variant of pooling
pulling <- function(Y_res, common_meta_dat, sig_e){
  res <- list()
  for (i in c(4:ncol(common_meta_dat))) {
    pool_res <- data.frame(do.call(rbind, list(apply(Y_res, 2, pools, index = common_meta_dat[i], sig_e = sig_e))))
    res<- append(res, list(pool_res))
  }
  return (res)
}

#Example
#test_exp4 <- pulling(test_exp2, meta_pools, 1)

#Function for taking into account technical repetition and technical varuability for data without pooling
#Y_result - vectors of values for samples
#sig_e - technical error (variability determed technical errors)
#return vactor of values with technical repetition and technical varuability (twice longer, because we use 2 technical repetition)
error_add <- function(Y_result, sig_e){
  Y_res <- rep(Y_result, each = 2)
  Y_res <- Y_res + rnorm(length(Y_res), mean = 0, sd = sig_e)
  return(Y_res)
}
#Example
#test_exp5 <- error_add(test_exp1, 1)

#Function. Combune work of previous functions: data_collection, pulling and errors add
#n_strings - number of spots (number of replication)
#beta_0 - level of expression in group 1 (base level)
#beta_1 - how expression in group 2 is differ from group 1
#sig_B - level of biological dispersion
#common_meta_dat - meta_data (we use metadata with information about pools (made by add_pools function))
#sig_e - technical error (variability determed technical errors)
#return: List containing list with matrix of values for each variant of pooling (including no pooling)
data_with_pools <- function(n_strings, beta_0, beta_1, sig_B, common_meta_dat, sig_e){
  Y_result <- data_collecting(n_strings, beta_0, beta_1, sig_B, common_meta_dat)
  Y_no_pool <- data.frame(apply(Y_result, 2, error_add, sig_e = sig_e))
  Y_matrix_pooled <- pulling(Y_result, common_meta_dat, sig_e)
  spot_data <- append(Y_matrix_pooled, list(Y_no_pool))
  return(spot_data)
}
#Example
#exp_test1 <- data_with_pools(n_strings = 100, beta_0 = 0, beta_1 = 1.5, sig_B = 1, common_meta_dat = meta_pools, sig_e = 0.5)


#Function. Moderated t-test (limma package)
#exp_sn - expression data list
#df_for_count - metadata (without pooling but with technical repetition)
#return power - percentage of cases where the difference was statistically significant
random_sample_limma <- function(exp_sn, df_for_count){
  biolrep <- df_for_count$B
  corfit <- duplicateCorrelation(exp_sn, ndups = 1, block = biolrep)
  X <- model.matrix(~ G, data = df_for_count)
  fit <- lmFit(exp_sn, design = X, block = biolrep, cor = corfit$consensus, method = "robust", maxit = 10000)
  efit <- eBayes(fit)
  a <- topTable(efit, coef = 2, number = nrow(exp_sn), adjust = "BH")
  value <- (sum(a$adj.P.Val<0.05))/nrow(exp_sn)
  return(value)
}

#Function - create metadata for each variant of pooling
#meta_pools - metadata with pools (made by add_pools function)
#Return list containng lists with metadata
meta_data_pools_generation <- function(meta_pools){
  no_pool <- meta_pools %>% dplyr::select(G, B, R)
  pool_start <- data.frame(rbind(apply(no_pool, 2, rep, each = 2)))
  no_pool <- pool_start
  no_pool$pool <- "No_pool"
  pool_end <- list()
  for (i in 4:(ncol(meta_pools))) {
    pool <- colnames(meta_pools)[i]
    res_pool <- unique(meta_pools %>% dplyr::select(G, R, pool))
    res_pool <- data.frame(rbind(apply(res_pool, 2, rep, each = 2)))
    colnames(res_pool)[3] <- "B"
    res_pool$pool <- pool
    pool_end <- append(pool_end, list(res_pool))
  }
  pool_end <- append(pool_end, list(no_pool))
  return (pool_end)
}

#Example
#pools_all <- meta_data_pools_generation(meta_pools)

#Function (Intermediate), allow to use list for exp_sn in random_sample_limma function
#list_dat - list containg data of the expression
#meta - list of metadata
list_to_limma <- function(list_dat, meta){
  df<- do.call(rbind, list_dat)
  power <- random_sample_limma(t(df), as.data.frame(meta))
  return(power)
}

#Function, combine previous function for lists of data
#Return vector with values of power for each variant of pooling
power_limma <- function(list_all_data, list_all_meta){
  res <- array()
  for (i in (1:length(list_all_data))){
    res <- append(res, list_to_limma(list_all_data[i], list_all_meta[i]))
  }
  return(res)
}

#Example
#exp_test2 <- power_limma(exp_test1, pools_all) #first NA - is OK

#Function/ Combuned most previous function (after meta_pools creration)
##Return vector with values of power for each variant of pooling
from_data_to_power <- function(n_strings, beta_0, beta_1, sig_B, common_meta_dat, sig_e){
  spot_data <- data_with_pools(n_strings, beta_0, beta_1, sig_B, common_meta_dat, sig_e)
  meta_data <- meta_data_pools_generation(common_meta_dat)
  result <- power_limma(spot_data, meta_data)
  return(result)
}

#Example
#exp_test3 <- from_data_to_power(100, 0, 2, 1, meta_pools, 1) #first NA - is OK

#Function which calculated mean and quantiles (2,5% and 97,5%)
#pow - vector of value
#return dataframe
mean_and_q <- function(pow){
  pow <- as.vector(pow)
  pow_mean <- mean(pow)
  perc_025 <- quantile(pow, probs = 0.025) 
  perc_975 <- quantile(pow, probs = 0.975)
  res_df <- data.frame(pow_mean = pow_mean, perc_025= perc_025, perc_975=perc_975)
  return(res_df)
}

#Function replicate from_data_to_power function times times and calculate mean and quantiles of power for each variant of pooling
#Warning - meta_pools is default value for argument common_meta_pools
#Return list
var_power <- function(times, n_strings, beta_0, beta_1, sig_B, common_meta_dat = meta_pools, sig_e){
  power_matrix <- RepParallel(times, from_data_to_power(n_strings, beta_0, beta_1, sig_B, common_meta_dat, sig_e))
  power_matrix <- power_matrix[2:nrow(power_matrix),]
  pow_variety<- apply(as.matrix(power_matrix), 1, mean_and_q)
  pow_variety_frame<- data.frame(do.call(rbind, do.call(cbind,pow_variety)))
  return(pow_variety_frame)
}

#Example
##Warning function could take some time!
#test3 <- var_power(20, 100, 0, 2, 1, meta_pools, 0.5)

#To check how power change with different start parameters. You can change them for ever you want
##Make dataframe contaning all combination of parameters
beta_0 <- 0 #base level of expression (expression in group 1)
beta_1 <- c(1,2,4,6) #difference between expression in groups
sig_B <- c(1, 1.5, 2.5, 5) #Biological dispersion
sig_e <- c(0.1, 0.5, 1, 1.5, 2) #technical variability
k <- 100 #number of spots creating (number of rows in each expression set)
times <- 100 #amount of times we repeat for each scenario variant
data_scen <- expand.grid(times = times, k = k, beta_0 = beta_0, beta_1 = beta_1, sig_B = sig_B, sig_e = sig_e) 
#

#Function - var_power applied for each row in data_scen dataframe and combined results into new dataframe
#data_scen - dataframe with scenario parameters for function
#common_meta_dat = meta_pools Default parameter for matedata with pools
#return dataframe
resulting_function <- function(data_scen, common_meta_dat = meta_pools){
  res <- mapply(var_power, data_scen$times, data_scen$k, data_scen$beta_0, data_scen$beta_1, data_scen$sig_B, sig_e = data_scen$sig_e)
  d_res <- as.data.frame(do.call(rbind,res))
  result <- list()
  for (i in (seq(from = 1, to = ncol(d_res), by = 3))){
    result <- append(result, list(data.frame(mean = d_res[,i], q_025 = d_res[,(i+1)], q925 = d_res[,(i+2)], pools = i))) 
  }
  result_df <- do.call(rbind, result)
  scen <- as.data.frame(apply(data_scen, 2, rep, ncol(d_res)/3))
  result_df <- cbind(scen, result_df)
  return(result_df)
}

#Warning function work quiet long (depend of technical parameters of computer and times and k in data_scen and also amount of rows at dat_scen)
result_test <- resulting_function(data_scen, common_meta_dat = meta_pools)
#in pools column at dataframe contain coded information about pools but not amount of sample in pool, which we want
samples_in_pool <- c(number_of_samples/number_of_pools, 1)
pool_level <- rep(samples_in_pool, each = nrow(data_scen))
result_test$pools_decode <- pool_level

#presentation of result at graph
sub_title <- "Parameters of generation: k = 100, n = 100 \nsig_e - technical variation, sig_B - biological variation \nbeta_1 - difference between expression in groups"
ggplot(result_test, aes(x = pools_decode, y = mean, color = factor(beta_1))) + 
  geom_errorbar(aes(ymax=q_025, ymin=q925, width=I(0.5), size=I(0.2)))+
  geom_point(aes(size = I(2)))+
  labs(title="Result of generation - Statistical power of moderated t-test",
       subtitle = sub_title) +
  facet_grid(sig_B~sig_e, labeller = label_both) + 
  scale_x_continuous(name='The number of samples in the pool')+
  scale_y_continuous(name='Statistical power of moderated t-test', labels = scales::percent, limits = c(0, 1))+
  scale_color_discrete(name = "Difference \nbetween \ngroups \n(beta 1)")+
  theme_bw()
