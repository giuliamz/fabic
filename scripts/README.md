This folder contains:

For experiment presentation:
* Fabic_el : script to present each run of the experiment (with eyelink)
* Fabic_aud : script to present the unisensory auditory run for the exclusion criteria in Experiment2

For model-free analyses:
* audloc_regress_fabic : script to analyse the results of the auditory localiser (RMSE)
* WAV_rep_fabic : script to compute the average audiovisual weight index
* WAV_rep_subj_spli : script to compute the average audiovisual weight index using participant-specific and condition-specific responses as scale factor (See supplementary materials)
* WAV_perm_test_fabic : script to perform permutation testing on the audiovisual weight index

For model-based analyses:
* exp1_run_model_pool, exp2_run_model_pool : script to run model estimation for BCI pooled, FF pooled, FC pooled and SF pooled
  - dependencies: BCI/bci_fit_model.m, BCI/fitModelFus.m, BCI/fixCrit_fitmodel.m and BCI/stoFus_fitmodel.m
* exp1_run_model_sep, exp2_run_model_sep : script to run model estimation for BCI separated, FF separated, FC separated and SF separated
  - dependencies: BCI/bci_fit_model.m, BCI/fitModelFus.m, BCI/fixCrit_fitmodel.m and BCI/stoFus_fitmodel.m
* exp1_model_comparison, exp2_model_comparison: script to run model comparison in the 2 (BCI vs FF model) x 2 (pooled vs separated) model space
* exp1_model_comparison_8mod, exp2_model_comparison_8mod : script to run model comparison in the 4 (BCI vs FF vs FC vs SF model) x 2 (pooled vs separated) model space
* exp1_coefficient_determination, exp2_coefficient_determination, exp1_coefficient_determination_8mod, exp2_coefficient_determination_8mod : script to compute the coefficient of determination (R^2)
  - dependencies: BCI/fitModelNull
* exp1_parameters_perm_testing, exp2_parameters_perm_testing : permutation testing of winning model parameters
* correlation_WAV_parameters : Correlation analysis between WAV and parameters winning model
* WAV_predicted : Compute WAV based on predictions of BCI and FC models (separated)
* correlation_params_bci_fixcrit : Correlation analysis pcommon and criterion parameters

Model-free vs Model-based:
* proportion_responses : Compute proportions of participants' responses and model predictions
* correlation_WAV_params : Correlation analysis WAV and BCI/FC model parameters 

Scripts are listed in order of execution.

Analysis of response times was performed on JASP using participants median response times per condition ('rt_experiment1.csv' and 'rt_experiment2.csv' in the data folder).
