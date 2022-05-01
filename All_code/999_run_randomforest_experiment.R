adjust_R2 = function(R2, n, k){
	#Function to convert R2 to adjusted R2
	#R2: R-squared accuracy of model
	#n: Number of samples in training set
	#k: Number of predictors in model
	
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	
	return(adjusted_R2)
}

RF_CV = function(input, output, CV, train_fraction, n_tree){
      
    print(CV)
    #Takes in and input data.table/dataframe and runs a random forest experiment
    
    #Convert input and output to data.table
    input = as.data.table(input)
    output = as.numeric(output)
    
    #Split data into training and testing 
    set.seed(CV)
    train_rows = sample(nrow(input), round(nrow(input)*train_fraction))
    test_rows = setdiff(seq(nrow(input)), train_rows)

    #Train randomForest
    rf = randomForest(x = input[train_rows], y = output[train_rows], keep.forest = T, 
					  xtest = input[test_rows], ytest = output[test_rows], ntree = n_tree, 
					  do.trace = F, importance = T)

    #Get test predictions
    test_R2 = (cor(rf$test$predicted, output[test_rows]))^2
    
	#Get adjusted R2
	adj_R2 = adjust_R2(test_R2, nrow(input), ncol(input))
	
    return(c(adj_R2, importance(rf, type = 1)))
}



run_RF_and_get_pred_importances = function(input, output, n_cv, n_tree=250, cores=8, train_split=0.8){
    
    #Runs a random forest regression experiment using Monte Carlo cross-validation n_cv times. Returns data table for R2 and predictor importances.
    
    #input: dataframe or table of input/predictors. Rows must correspond to observations and columns must correspond to predictors.
    #Output: Numeric vector of response/true model output, with length equal to nrows of input.
    #n_cv: Number of Monte-Carlo cross-validations to be performed.
    #cores= Number of CPU cores to use for parallelization.
    #train_split= Fraction to split data into training and testing sets.
    
    #Returns data table where each row is a cross-validation run, the first column is the Test_R2 for that run, 
    #and the following columns are predictor importances for that run.
    
    
    #Run random forest for each CV and make into dataframe where each row is a single CV run (shows progress bar)
    RF_results = as.data.table(do.call("rbind.data.frame", 
									   mclapply(seq(n_cv), 
												function(i) RF_CV(input, output, i, train_split, n_tree), 
												mc.cores = cores)))

    colnames(RF_results) = c("Test_adj_R2", colnames(input))
    
    return(RF_results)                                                             
                                                                 
                                                                     

}

