# CDNFS
This is the R realization of the feature selection method proposed in *CDNFS: A New Filter Method for 
Nonlinear Feature Selection in High Dimensional Data*.

In the paper, we raised two simulation studies to demonstrate  the effectiveness, efficiency, and interpretability of the proposed CDNFS in handling both regression and classification problems. In addition, experiments on four real datasets were adopted to show the predictive capability and interpretability of  CDNFS combined with some popular learners.  

The source code is provided in `code` folder, and the generated simulations are in `data` folder. The versions of R 
packages are listed in `README.md`.

The `code` folder include two parts, one is for `experiments`, and another is for `evaluation`.

* ***Experiments***: Files for running the experiments on simulations and real datasets

    * **main_simulation_data_generation.R**: generation of the simulation examples

    * **main_simulation_regression_FS_performance.R**: Feature selection performance on the regression simulation,  create *regression_simulation\_variables\_selected\_(learners)\_(IterationTimes).csv* under the path `../FS_performance/regression_simulation/`
     * **main_simulation_classification_FS_performance.R**: Feature selection performance on the classification simulation,  create *classification_simulation\_variables\_selected\_(learners)\_(IterationTimes).csv* under the path `../FS_performance/classification_simulation/`
    * **main_simulation_regression_predictive_performance.R**: The out-of-sample MSE of 5-fold cross validation on the regression simulation, create *result_regression_simulation\_(FeatureSelectionMethod)\_(IterationTimes)\_(FoldNumber).csv* under the path `../regression_simulation/repeat_(IterationTimes)/`
    * **main_simulation_classification_predictive_performance.R**: The out-of-sample accuracy of 5-fold cross validation on the classification simulation, create *result_classification_simulation\_(FeatureSelectionMethod)\_(IterationTimes)\_(FoldNumber).csv* 
    under the path `../classification_simulation/repeat_(IterationTimes)/`
 
    * **main_realdata_regression_predictive_performance.R**: Prediction performance on the regression real dataset, create *result\_realdata\_regression\_(FeatureSelectionMethod)\_(DatasetName)\_(FoldNumber).csv* under the path `../(DatasetName)/repeat_(InterationTimes)/`
    * **main_realdata_classification_predictive_performance.R**: Prediction performance on the classification real datasets, create *result\_realdata\_classification\_(FeatureSelectionMethod)\_(DatasetName)\_(FoldNumber).csv* under the path `../(DatasetName)/(InterationTimes)/`
    * **main_realdata_visualization.R**: Visualization of the feature selection result on real datasets, create tables of the uninformative, relevant, and redundant features detected in each interation under the path `../visualization/(DatasetName)/`
    * **function_CDNFS.R**: All the functions for CDNFS and baselines, and the learners
   
    * **neuralnet_1.44.2.tar.gz**: Version of neuralnet that was used to run all experiments

> The file **function_CDNFS.R** are sourced in all the main functions above.

* ***Evaluation***: Files for generation of tables and plots
    * **main_simulation_evaluation.R**: The feature selection results (S_min, P_{S_i}, P_a and time) of filter methods on simultations, create *evaluation\_(min/p)\_(regression/classification)\_simulation.csv*
    
    * **simulation_result_analysis_CV.R**: The summary of the prediction results on simulations, create *(MSE/ACC)\_outsample\_(learner)\_(regression/classification)_simulation.csv*
    * **realdata_result_analysis_repeated_CV.R**: The summary of the prediction results on real datasets, create *(MSE/ACC)\_outsample\_(learner)\_(DatasetName).csv*
    
    * **read_four_types_of_features.R**: Feature selection result of CDNFS on real datasets, create *(DatasetName)_CDNFS_selected_variables.csv*
    * **ORL_visual_analysis.R**: Visualization of the feature selection result of CDNFS on the ORL dataset, create 400 images 

**Notice**: The programs in `Experiments` can be executed directly, and the programs in `Evaluation` are executed based on the results generated with programs in `Experiments`.

---
**Environment requirement**:

    R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
    Platform: x86_64-w64-mingw32/x64 (64-bit)

**Versions of R packages**:

    cowplot 1.1.0
    e1071 1.7-4
    FOCI 0.1.2
    ggplot2 3.3.2
    infotheo 1.2.0
    lightgbm 3.1.1.99
    maptree 1.4-7
    neuralnet 1.44.2
    NeuralNetTools 1.5.2
    nlnet 1.4
    plyr 1.8.6
    reshape2 1.4.4
    sampling 2.9
    stringr 1.4.0
    rpart 4.1-15
    
    



