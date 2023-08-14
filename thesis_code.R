library(dplyr)
library(grf)
library(fastDummies)

data <- readxl::read_excel('artea.xlsx', sheet = 'AB_test')

# One-hot encoding categorical channel of acquisition variable
data <- data |>
  dummy_cols(select_columns = c("channel_acq"), remove_selected_columns = TRUE) 

summary(data)

# Defining outcome and treatment variables and the covariates
y <- data$trans_after 
T <- data$test_coupon 
X <- data[, c('channel_acq_1', 'channel_acq_2', 'channel_acq_3', 'channel_acq_4', 'channel_acq_5', 'num_past_purch', 'spent_last_purchase', 'weeks_since_visit', 'browsing_minutes', 'shopping_cart')]

df <- data.frame(y = y, T = T)

# Exploratory analysis ----------------------------------------------------

# Calculating the summary statistics of each treatment condition
summary_stats <- df |> 
  group_by(T) |> 
  summarise( mean_y = mean(y), median_y = median(y), sd_y = sd(y), min_y = min(y), max_y = max(y), Q1_y = quantile(y, 0.25), Q3_y = quantile(y, 0.75), n = n() )

# Computing the ATE and p-value
ate <- summary_stats$mean_y[summary_stats$T == 1] - summary_stats$mean_y[summary_stats$T == 0]
t_test_result <- t.test(y ~ T, data = df)
p_value <- t_test_result$p.value

# Adding the ATE and p-value as new variables in summary_stats
summary_stats$ATE <- ifelse(summary_stats$T == 1, ate, NA)
summary_stats$p_value <- ifelse(summary_stats$T == 1, p_value, NA)

print(summary_stats)

# Calculating the average order value during the experimentation period
total_trans <- sum(data$trans_after)
total_rev <- sum(data$revenue_after)
aov_exp <- total_rev/total_trans

# Calculating the total discount costs
discount <- 0.2
discount_cost <- discount*aov_exp

set.seed(716) 

# Splitting the data into three equal random samples
n <- nrow(data)
sample_size <- n / 3
index <- sample(1:n, size = n, replace = FALSE)
exp_index <- index[1:sample_size] # experimentation sample index
al_index <- index[(sample_size + 1):(2 * sample_size)] # active learning sample index
hold_index <- index[(2 * sample_size + 1):n] # holdout sample index

# Creating variables for covariates, outcome and treatment of experimentation sample
X_exp <- as.matrix(X[exp_index, ]) 
y_exp <- y[exp_index]
T_exp <- T[exp_index] # use pre-existing assignment (assuming random assignment)

df_exp <- data.frame(y = y_exp, T = T_exp)

# Creating a data frame for the active learning sample
df_al <- data.frame(y = y[al_index], T = T[al_index])

# Creating a data frame for the holdout sample
df_hold <- data.frame(y = y[hold_index], T = T[hold_index])

# Calculating the summary statistics of each treatment condition for each sample
summary_stats_exp <- df_exp |> group_by(T) |> summarise( mean_y = mean(y), median_y = median(y), sd_y = sd(y), min_y = min(y), max_y = max(y), Q1_y = quantile(y, 0.25), Q3_y = quantile(y, 0.75), n = n() )
summary_stats_al <- df_al |> group_by(T) |> summarise( mean_y = mean(y), median_y = median(y), sd_y = sd(y), min_y = min(y), max_y = max(y), Q1_y = quantile(y, 0.25), Q3_y = quantile(y, 0.75), n = n() ) 
summary_stats_hold <- df_hold |> group_by(T) |> summarise( mean_y = mean(y), median_y = median(y), sd_y = sd(y), min_y = min(y), max_y = max(y), Q1_y = quantile(y, 0.25), Q3_y = quantile(y, 0.75), n = n() )

summary_stats_exp$sample <- "exp"
summary_stats_al$sample <- "al"
summary_stats_hold$sample <- "hold"

# Combining the summary statistics of each sample into one table
summary_stats_all <- rbind(summary_stats_exp, summary_stats_al, summary_stats_hold)
print(summary_stats_all)

# Checking whether the assignment of the randomly selected units is equally split
freq_table_exp <- table(T_exp) # frequency table of treatment and control labels 
prop_table_exp <- prop.table(freq_table_exp) # proportion table of treatment and control labels

# Model training ----------------------------------------------------------

# Pre-fitting experimentation model for Y and W 
forest.Y <- regression_forest(X_exp, y_exp, tune.parameters = "all") 
forest.W <- regression_forest(X_exp, T_exp, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions 
W.hat <- predict(forest.W)$predictions

# Fitting the causal forest on the experimentation sample 
causal_forest <- grf::causal_forest(X = X_exp, Y = y_exp, W = T_exp, Y.hat = Y.hat, W.hat = W.hat, num.trees = 5000, tune.parameters = "all")

# Creating variables for covariates, outcome and treatment of active learning sample
X_al <- as.matrix(X[al_index, ]) 
y_al <- y[al_index]  
T_al <- T[al_index] 

# Predicting the CATEs and standard errors for the units in the active learning sample
var_pred <- predict(causal_forest, X_al, estimate.variance = TRUE) # predicted CATEs and variances 
se_pred <- sqrt(var_pred$variance.estimates) # predicted standard errors

# Selecting the units based on a standard error cutoff of the estimated CATEs
quantile_80 <- quantile(se_pred, 0.8) # 80th percentile of standard errors 
selected_index <- al_index[se_pred >= quantile_80] # selected units index

# Creating variables for covariates, outcome and treatment of subsample of active learning sample
X_sel <- as.matrix(X[selected_index, ]) 
y_sel <- y[selected_index]  
T_sel <- T[selected_index] 

# Checking whether the assignment of the AL selected units is equally split
freq_table_al <- table(T_sel)  # frequency table of treatment and control labels
prop_table_al <- prop.table(freq_table_al) # proportion table of treatment and control labels (aiming for 0.4-0.6)

# Adding the selected units and labels to the initial experimentation sample
X_combined <- rbind(X_exp, X_sel)  
y_combined <- c(y_exp, y_sel) 
T_combined <- c(T_exp, T_sel) 

# Pre-fitting experimentation model for Y and W 
forest.Y <- regression_forest(X_combined, y_combined, tune.parameters = "all") 
forest.W <- regression_forest(X_combined, T_combined, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions 
W.hat <- predict(forest.W)$predictions

# Fitting the causal forest on the experimentation sample and actively selected subsample  
causal_forest_combined <- grf::causal_forest(X = X_combined, Y = y_combined, W = T_combined, Y.hat = Y.hat, W.hat = W.hat, num.trees = 5000, tune.parameters = "all")

# Creating variables for covariates, outcome and treatment of holdout sample
X_hold <- as.matrix(X[hold_index, ]) 
y_hold <- y[hold_index]  
T_hold <- T[hold_index] 

# Predicting CATEs and standard errors for the units in the holdout sample
var_pred_hold <- predict(causal_forest_combined, X_hold, estimate.variance = TRUE) # predicted CATEs and variances

# Counting the number of positive and negative predicted CATEs
positive_count <- sum(var_pred_hold$predictions > 0) 
negative_count <- sum(var_pred_hold$predictions < 0)

# Defining the benchmark models --------------------------------------------------------------

# First benchmark - predicting the CATEs and variances by fitting a causal forest only on the randomly selected experimentation sample
var_pred_exp <- predict(causal_forest, X_hold, estimate.variance = TRUE) 

# Second benchmark - predicting the CATEs and variances by fitting a causal forest on the experimentation sample and a randomly selected subsample of the AL sample with the same number of units as actively selected subsample 
n_random <- length(selected_index)  # size of random sample same as selected units sample
random_index <- sample(al_index, size = n_random, replace = FALSE)  # randomly select units from AL sample 

# Creating variables for covariates, outcome and treatment randomly selected subample of AL sample
X_random <- as.matrix(X[random_index, ])
y_random <- y[random_index]  
T_random <- T[random_index]  

# Adding the randomly selected units to the initial experimentation sample
X_combined_random <- rbind(X_exp, X[random_index, ]) 
y_combined_random <- c(y_exp, y[random_index]) 
T_combined_random <- c(T_exp, T[random_index]) 

# Pre-fitting experimentation model for Y and W 
forest.Y <- regression_forest(X_combined_random, y_combined_random, tune.parameters = "all") 
forest.W <- regression_forest(X_combined_random, T_combined_random, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions 
W.hat <- predict(forest.W)$predictions

# Fitting the causal forest on the experimentation sample and randomly selected subsample 
causal_forest_random <- grf::causal_forest(X = X_combined_random, Y = y_combined_random, W = T_combined_random, Y.hat = Y.hat, W.hat = W.hat, num.trees = 5000, tune.parameters = "all")

# Predicting the CATEs and variances of the units in the holdout sample
var_pred_random <- predict(causal_forest_random, X_hold, estimate.variance = TRUE)  

# Specifying the evaluation metric (Qini) ---------------------------------------------

# Defining the function to calculate the Qini coefficient and curve for a given model
Qini_model <- function(model, X_hold, y_hold, T_hold, n = 10) {
  
  # Predicting the CATEs for the holdout sample
  cate_pred <- predict(model, X_hold)$predictions

  # Creating a data frame with the predicted CATEs and the true outcome and treatment
  df <- data.frame(cate_pred = cate_pred, y = y_hold, T = T_hold)
  
  # Creating groups based on the predicted CATEs
  df$predrank <- rank(-df$cate_pred) # rank the predicted CATEs in descending order
  brk <- unique(quantile(df$predrank, probs = seq(0, 1, 1 / n))) # create n equal intervals based on the rank
  df$group <- cut(df$predrank, breaks = brk, labels = NULL, include.lowest = TRUE) # assign each unit to a group
  
  # Calculating outcomes per group
  Rc <- tapply(df[df$T == 0, ]$y, df[df$T == 0, ]$group, sum) # sum of outcomes for control units in each group
  Rt <- tapply(df[df$T == 1, ]$y, df[df$T == 1, ]$group, sum) # sum of outcomes for treated units in each group
  RcMean <- tapply(df[df$T == 0, ]$y, df[df$T == 0, ]$group, mean) # mean of outcomes for control units in each group
  RtMean <- tapply(df[df$T == 1, ]$y, df[df$T == 1, ]$group, mean) # mean of outcomes for treated units in each group
  Nc <- tapply(df[df$T == 0, ]$y, df[df$T == 0, ]$group, length) # number of control units in each group
  Nt <- tapply(df[df$T == 1, ]$y, df[df$T == 1, ]$group, length) # number of treated units in each group
  
  # Combining the outcomes per group into a data frame
  PM <- merge(cbind(Rc ,RcMean ,Nc), cbind(Rt ,RtMean ,Nt), by ="row.names", all = TRUE)
  PM$Row.names<- as.numeric(PM$Row.names)
  PM[, c(2 ,4 ,5 ,7)][is.na(PM[, c(2 ,4 ,5 ,7)])] <-0 #replace missing values with zero (implying zero counts)
  
  PM<- PM[order(PM$Row.names),] #order the data frame by group number
  
  res<- cbind(
    group= PM$Row.names,
    Nt= PM$Nt,
    Nc= PM$Nc,
    Rt= PM$Rt,
    Rc= PM$Rc,
    RtMean= PM$RtMean,
    RcMean= PM$RcMean
  )
  
  # Cumulative incremental gain
  PM$cig<- (res[, "Rt"] - res[, "Rc"] * sum(res[, "Nt"]) / sum(res[, "Nc"])) / sum(res[, "Nt"]) # incremental gain per group
  PM$cig_cum<- cumsum(PM$cig) # cumulative incremental gain across groups
  
  # Random cumulative incremental gain
  Ioa<- sum(res[, "Rt"]) / sum(res[, "Nt"]) - sum(res[, "Rc"]) / sum(res[, "Nc"]) # overall incremental gain
  PM$Ioa<- Ioa
  groupNum<- nrow(PM)
  PM$Irand1<- Ioa / groupNum # random incremental gain per group
  PM$Irand_cum<- cumsum(PM$Irand1) # random cumulative incremental gain across groups
  
  # Area under the uplift curve
  x<- seq(1 / groupNum ,1 ,1 / groupNum)
  y<- PM$cig_cum
  AUC<- sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)])) /2
  
  # Area under the random curve
  y.randcumincrgain<- PM$Irand_cum
  AUCrand<- sum((x[-1] - x[-length(x)]) * (y.randcumincrgain[-1] + y.randcumincrgain[-length(y)])) /2
  
  # Difference of the areas: Qini-coefficient
  Qini<- AUC - AUCrand
  
  # Defining the output of the function
  res<- cbind(
    group= PM$Row.names,
    Nt= PM$Nt,
    Nc= PM$Nc,
    Rt= PM$Rt,
    Rc= PM$Rc,
    RtMean= PM$RtMean,
    RcMean= PM$RcMean,
    cumincrgain= PM$cig_cum,
    randcumincrgain= PM$Irand_cum,
    Qini,
    groupNum
  )
  
  res<- round(res, 6)
  
  miny<- min(c(res[, "randcumincrgain"], res[, "cumincrgain"]))
  maxy<- max(c(res[, "randcumincrgain"], res[, "cumincrgain"]))
  plot(
    res[, "cumincrgain"] ~ seq(100 / res[1 , "groupNum"],100 ,100 / res[1 , "groupNum"]),
    type ="b",
    col ="blue" ,
    lty=2 ,
    xlab ="Proportion of population targeted (%)",
    ylab ="Cumulative incremental gains",
    ylim = c(miny ,maxy)
  )
  lines(
    res[, "randcumincrgain"] ~ seq(100 / res[1 , "groupNum"],100 ,100 / res[1 , "groupNum"]),
    type ="l",
    col ="red" ,
    lty=1 
  )
  legend(
    "topright",
    c("Model", "Random Allocation"),
    col=c("blue", "red"),
    lty=c(2 ,1)
  )
  return(as.data.frame(res))
}

# Testing the models --------------------------------

# AL model 
qini_AL <- Qini_model(model = causal_forest_combined,X_hold= X_hold,y_hold= y_hold,T_hold= T_hold,n=10)

# First benchmark - randomly selected experimentation sample
qini_exp<- Qini_model(model= causal_forest,X_hold= X_hold,y_hold= y_hold,T_hold= T_hold,n=10)

# Second benchmark - randomly selected experimentation sample + randomly selected subsample 
qini_random<- Qini_model(model= causal_forest_random,X_hold= X_hold,y_hold= y_hold,T_hold= T_hold,n=10)

# Performance evaluation of the models --------------------------------

# Comparing the Qini coefficients 
qini_table<- data.frame(Model=c("AL", "First benchmark", "Second benchmark"), Qini_coefficient=c(qini_AL[nrow(qini_AL),]$Qini,qini_exp[nrow(qini_exp),]$Qini,qini_random[nrow(qini_random),]$Qini))
print(qini_table)

# Plotting the Qini curves 
plot(qini_AL$cumincrgain ~ seq(100 / qini_AL[1,"groupNum"],100 ,100 / qini_AL[1,"groupNum"]),type="b",col="blue" ,lty=2,xlab="Proportion of population targeted (%)",ylab="Cumulative incremental gains")
lines(qini_exp$cumincrgain ~ seq(100 / qini_exp[1,"groupNum"],100 ,100 / qini_exp[1,"groupNum"]),type="l",col="red" ,lty=1)
lines(qini_random$cumincrgain ~ seq(100 / qini_random[1,"groupNum"],100 ,100 / qini_random[1,"groupNum"]),type="l",col="green" ,lty=3)
legend("bottomright",legend=c("AL", "First benchmark", "Second benchmark"),col=c("blue", "red","green"),lty=c(2 ,1 ,3))

# Model tuning  --------------------------------------------------------

set.seed(716)

# Defining a function to perform cross-validation for a given number of trees and standard error cutoff
cv_model <- function(X, y, T, num_trees, prob, exp_index, al_index, eval_X, eval_y, eval_T) {
  
  # Creating variables for covariates, outcome and treatment of experimentation sample
  X_exp <- as.matrix(X[exp_index, ]) 
  y_exp <- y[exp_index] 
  T_exp <- T[exp_index]
  
  # Creating variables for covariates, outcome and treatment of active learning sample
  X_al <- as.matrix(X[al_index, ]) 
  y_al <- y[al_index] 
  T_al <- T[al_index]
  
  # Pre-fitting experimentation model for Y and W
  forest.Y <- regression_forest(X_exp, y_exp, tune.parameters = "all", num.trees = num_trees) 
  forest.W <- regression_forest(X_exp, T_exp, tune.parameters = "all", num.trees = num_trees) 
  Y.hat <- predict(forest.Y)$predictions 
  W.hat <- predict(forest.W)$predictions
  
  # Fitting the causal forest on the experimentation sample  
  causal_forest <- grf::causal_forest(X= X_exp,Y= y_exp,W= T_exp,Y.hat= Y.hat,W.hat= W.hat,tune.parameters="all",num.trees= num_trees)
  
  # Predicting the CATEs and standard errors for the units in the active learning sample
  var_pred_al<- predict(causal_forest,X_al ,estimate.variance=TRUE) 
  se_pred_al<- sqrt(var_pred_al$variance.estimates) 
  
  # Selecting the units based on a standard error cutoff of the estimated CATEs
  quantile_al <- quantile(se_pred_al, prob)  
  selected_index_al<- al_index[se_pred_al >= quantile_al] 
  
  # Creating variables for covariates, outcome and treatment of subsample of active learning sample
  X_sel_al<- as.matrix(X[selected_index_al, ]) 
  y_sel_al<- y[selected_index_al]
  T_sel_al<- T[selected_index_al]
  
  # Adding the selected units and labels to the initial experimentation sample
  X_combined_al<- rbind(X_exp,X_sel_al)
  y_combined_al<- c(y_exp,y_sel_al) 
  T_combined_al<- c(T_exp,T_sel_al)
  
  # Pre-fitting combined model for Y and W
  forest.Y_combined <- regression_forest(X_combined_al, y_combined_al, tune.parameters="all", num.trees= num_trees)
  forest.W_combined <- regression_forest(X_combined_al, T_combined_al, tune.parameters="all", num.trees= num_trees)
  Y.hat_combined <- predict(forest.Y_combined)$predictions
  W.hat_combined <- predict(forest.W_combined)$predictions
  
  # Fitting the causal forest on the experimentationn sample and actively selected subsample    
  causal_forest_combined <- grf::causal_forest(X= X_combined_al,Y= y_combined_al,W= T_combined_al,Y.hat= Y.hat_combined,W.hat= W.hat_combined,tune.parameters="all",num.trees= num_trees)
  
  # Compute the Qini coefficient for the model using the predictions and the true outcomes of the evaluation fold
  qini_model <- Qini_model(model = causal_forest_combined, X_hold = eval_X, y_hold = eval_y, T_hold = eval_T, n = 10) 
  Qini <- qini_model[nrow(qini_model),]$Qini
  
  # Return the Qini coefficient for this fold
  return(Qini)
}

# Defining a list of possible number of trees to try
tree_list <- seq(1000,2000,500)

# Defining a list of possible standard error cutoffs to try
quantile_list <- seq(0.7, 0.9, 0.1)

# Creating a data frame to store the Qini coefficients for each number of trees and standard error cutoff combination
Qini_table <- expand.grid(num_trees = tree_list, prob = quantile_list, Qini = NA)

# Defining number of folds for cross-validation (k=5) 
k <- 5 
n <- nrow(X) 
index <- sample(1:n, size = n, replace = FALSE) 
fold_size <- floor(n / k)

# Iterating over Qini table
for (i in 1:nrow(Qini_table)) {
  
  # Extracting the number of trees and standard error cutoff for the current iteration
  num_trees <- Qini_table$num_trees[i] 
  prob <- Qini_table$prob[i]
  
  # Performing cross-validation for the current number of trees and standard error cutoff and calculating the average Qini coefficient
  Qini_cv <- rep(NA, k) # vector to store the Qini coefficients for each fold
  
  # Iterating over folds
  for (j in 1:k) {
    
    # Splitting the data into training and validation sets for the current fold
    if (j < k) {
      holdout_index <- index[((j - 1) * fold_size + 1):(j * fold_size)]
    } else {
      # For the last fold, take all remaining data points
      holdout_index <- index[((j - 1) * fold_size + 1):length(index)]
    }
    train_index <- index[-holdout_index]
    
    # Creating variables for covariates, outcome and treatment of training set
    X_train <- as.matrix(X[train_index, ])
    y_train <- y[train_index]
    T_train <- T[train_index]

    # Creating variables for covariates, outcome and treatment of holdout set
    X_holdout <- as.matrix(X[holdout_index, ])
    y_holdout <- y[holdout_index]
    T_holdout <- T[holdout_index]
    
    # Splitting training data in experimentation and AL sample
    half_len <- floor(length(train_index) / 2)
    abs_exp_index <- train_index[1:half_len]
    abs_al_index <- train_index[(half_len + 1):length(train_index)]
    rel_exp_index <- 1:half_len
    rel_al_index <- (half_len + 1):length(train_index)
    
    # Performing cross-validation for the current number of trees and standard error cutoff and fold
    Qini_cv[j] <- cv_model(X = X_train, y = y_train, T = T_train, num_trees = num_trees, prob = prob, exp_index = rel_exp_index, al_index = rel_al_index, eval_X = X_holdout, eval_y = y_holdout, eval_T = T_holdout)
  }
  
  # Calculating the average Qini coefficient across folds
  Qini_mean <- mean(Qini_cv)
  
  # Storing the average Qini coefficient in the Qini table
  Qini_table$Qini[i] <- Qini_mean
  
}

# Printing the Qini table for each number of tree and standard error cutoff combination
print(Qini_table)

# Finding the best number of trees and standard error cutoff based on the highest Qini coefficient
best_row <- which.max(Qini_table$Qini) 
best_num_trees <- Qini_table$num_trees[best_row] 
best_prob <- Qini_table$prob[best_row] 

# Printing the best Qini coefficient, number of trees, standard error cutoff
print(Qini_table$Qini[best_row])
print(best_num_trees)
print(best_prob)

# Using aggregate to compute mean Qini per number of trees
average_Qini <- aggregate(Qini ~ num_trees, data = Qini_table, FUN = mean)
print(average_Qini)

# Using aggregate to compute mean Qini per percentage standard error cutoff
average_Qini_per_prob <- aggregate(Qini ~ prob, data = Qini_table, FUN = mean)
print(average_Qini_per_prob)

# Testing the models with the optimised hyperparameters  ------------------------------------------------

# Splitting the data into experimentation, active learning and holdout samples using predefined indices
n <- nrow(X) 
sample_size <- n / 3 

# Creating vectors to store Qini coefficients
qini_AL_tuned_vec <- c() 
qini_exp_tuned_vec <- c() 
qini_random_tuned_vec <- c() 

# Setting the seeds for the model iterations, using the pre-cross-validation seed for the first iteration
seeds <- c(716, 212, 257, 879, 892, 136, 664, 712, 263, 821)
# seeds <- c(716, sample(1:1000, 9))

# Running the model 10 times with optimised parameters by iterating over the seeds
for (i in 1:10) {
    set.seed(seeds[i])
    index <- 1:n  # resetting the index to its original state before shuffling
    index <- sample(index, size = n, replace = FALSE) 
    exp_index <- index[1:sample_size]  
    al_index <- index[(sample_size + 1):(2 * sample_size)]  
    hold_index <- index[(2 * sample_size + 1):n] 
  
  # Creating variables for covariates, outcome and treatment of experimentation sample
  X_exp <- as.matrix(X[exp_index, ]) 
  y_exp <- y[exp_index]  
  T_exp <- T[exp_index] 
  
  # Pre-fitting experimentation model for Y and W for first benchmark
  forest.Y<- regression_forest(X_exp,y_exp,tune.parameters="all") 
  forest.W<- regression_forest(X_exp,T_exp,tune.parameters="all") 
  Y.hat<- predict(forest.Y)$predictions 
  W.hat<- predict(forest.W)$predictions
  
  # Fitting the causal forest on the experimentation sample  
  causal_forest_exp_tuned <- grf::causal_forest(X= X_exp,Y= y_exp,W= T_exp,Y.hat= Y.hat,W.hat= W.hat,tune.parameters="all",num.trees= best_num_trees)
  
  # Creating variables for covariates, outcome and treatment of active learning sample
  X_al <- as.matrix(X[al_index, ])
  y_al <- y[al_index]
  T_al <- T[al_index]
  
  # Calculating the Qini coefficient for the first benchmark
  qini_exp_tuned <- Qini_model(model= causal_forest_exp_tuned,X_hold= X_hold,y_hold= y_hold,T_hold= T_hold,n=10) 
  
  # Predicting the CATEs and standard errors for the units in the AL sample
  var_pred_al_tuned <- predict(causal_forest_exp_tuned,X_al ,estimate.variance=TRUE) 
  se_pred_al_tuned <- sqrt(var_pred_al_tuned$variance.estimates) 

  # Selecting the units based on the optimal standard error cutoff of the estimated CATEs
  quantile_al_tuned <- quantile(se_pred_al_tuned, best_prob) # quantile of standard errors best_prob
  selected_index_al_tuned <- al_index[se_pred_al_tuned >= quantile_al_tuned] # selected units index
  
  # Creating variables for covariates, outcome and treatment of subsample of active learning sample
  X_sel_al <- as.matrix(X[selected_index_al_tuned, ]) 
  y_sel_al <- y[selected_index_al_tuned]
  T_sel_al <- T[selected_index_al_tuned]
  
  # Adding the selected units and labels to the initial experimentation sample
  X_combined <- rbind(X_exp,X_sel_al)
  y_combined <- c(y_exp,y_sel_al) 
  T_combined <- c(T_exp,T_sel_al)
  
  # Pre-fitting combined model for Y and W
  forest.Y<- regression_forest(X_combined,y_combined,tune.parameters="all") 
  forest.W<- regression_forest(X_combined,T_combined,tune.parameters="all") 
  Y.hat<- predict(forest.Y)$predictions 
  W.hat<- predict(forest.W)$predictions
  
  # Fitting the causal forest on the experimentation sample and actively selected subsample  
  causal_forest_AL_tuned <- grf::causal_forest(X= X_combined,Y= y_combined,W= T_combined,Y.hat= Y.hat,W.hat= W.hat,tune.parameters="all",num.trees= best_num_trees)
  
  # Calculating the Qini coefficient for the AL model
  qini_AL_tuned <- Qini_model(model= causal_forest_AL_tuned,X_hold= X_hold,y_hold= y_hold,T_hold= T_hold,n=10)
  
  # Selecting random subsample of the AL sample
  n_random <- length(selected_index_al_tuned)  # size of random sample same as selected units sample
  random_index <- sample(al_index, size = n_random, replace = FALSE)  # randomly select units from AL sample 
  
  # Creating variables for covariates, outcome and treatment of randomly selected subsample of AL sample
  X_random <- as.matrix(X[random_index, ])
  y_random <- y[random_index]  
  T_random <- T[random_index]  
  
  # Adding the randomly selected units to the initial experimentation sample
  X_combined_random <- rbind(X_exp, X_random) 
  y_combined_random <- c(y_exp, y_random) 
  T_combined_random <- c(T_exp, T_random) 
  
  # Pre-fitting experimentation model for Y and W (for benchmark 2)
  forest.Y<- regression_forest(X_combined_random,y_combined_random,tune.parameters="all") 
  forest.W<- regression_forest(X_combined_random,T_combined_random,tune.parameters="all") 
  Y.hat<- predict(forest.Y)$predictions 
  W.hat<- predict(forest.W)$predictions
  
  # Fitting the causal forest on the experimentation sample and randomly selected subsample  
  causal_forest_combined_tuned <- grf::causal_forest(X= X_combined_random,Y= y_combined_random,W= T_combined_random,Y.hat= Y.hat,W.hat= W.hat,tune.parameters="all",num.trees= best_num_trees)
  
  # Calculating the Qini coefficient
  qini_random_tuned <- Qini_model(model= causal_forest_combined_tuned,X_hold= X_hold,y_hold= y_hold,T_hold= T_hold,n=10) 
  
  # Storing Qini coefficients in vectors
  qini_AL_tuned_vec <- c(qini_AL_tuned_vec, qini_AL_tuned[nrow(qini_AL_tuned),]$Qini)
  qini_exp_tuned_vec <- c(qini_exp_tuned_vec, qini_exp_tuned[nrow(qini_exp_tuned),]$Qini)
  qini_random_tuned_vec <- c(qini_random_tuned_vec, qini_random_tuned[nrow(qini_random_tuned),]$Qini)
  
}

# Performance evaluation of the tuned models --------------------------------------------------

# Creating a data frame with Qini coefficients per iteration for each model 
iteration <- 1:10
qini_df_tuned <- data.frame(iteration,seeds, qini_AL_tuned_vec, qini_exp_tuned_vec, qini_random_tuned_vec)
names(qini_df_tuned) <- c( "Iteration","Seed", "AL", "First benchmark", "Second benchmark")
print(qini_df_tuned)

# Calculating the mean of the Qini coefficients
mean_qini_AL <- mean(qini_df_tuned$AL)
mean_qini_exp <- mean(qini_df_tuned$'First benchmark')
mean_qini_random <- mean(qini_df_tuned$'Second benchmark')

# Comparing the Qini coefficients using statistical test
t_test_AL_exp <- t.test(qini_AL, qini_exp) # first benchmark 
t_test_AL_random <- t.test(qini_AL, qini_random) # second benchmark 

# Printing the results
print(round(mean_qini_AL, 7))
print(round(mean_qini_exp, 7))
print(round(mean_qini_random, 7))
print(round(t_test_AL_exp$p.value, 7))
print(round(t_test_AL_random$p.value, 7))



