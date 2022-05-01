library(reticulate)
shap = import("shap")
sklearn_ensemble=import("sklearn.ensemble")
plt = import('matplotlib.pyplot')

run_RF_and_get_SHAP=function(input, output, n_tree=250L){
    
    #Takes in and input data.table/dataframe and runs a random forest experiment
    
    #Convert input and output to data.table
    input=data.matrix(input)
    output=as.numeric(output)
    

    #Set up model and hyperparameters
    model = sklearn_ensemble$RandomForestRegressor(n_estimators=1000L, 
                               bootstrap = T,
                               max_features = 'sqrt',
                              oob_score=T)
    #Train model
    model$fit(input, output)
    
    #Set up SHAP for model
    explainer=shap$TreeExplainer(model)
    
    #Get SHAP values for each observation
    shap_values=as.data.table(explainer$shap_values(input))
    colnames(shap_values)=colnames(input)                                                      
    
    return(shap_values)                                                                 

}
