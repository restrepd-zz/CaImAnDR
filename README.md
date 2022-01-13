# CaImAnDR

Code to process GCaMP imaging data with CaImAn for experiments described in Ma et al.: doi: https://doi.org/10.1101/2019.12.14.876201

The code uses CaImAn to analyze GCaMP dF/F imaging data.

## Getting Started


### Prerequisites

The code uses the following MATLAB toolboxes:

Image
Parallel Computing
Signal
Statistics


You also need to have CaImAn installed because drgCaImAn_script.m uses some CaImAn functions

https://github.com/flatironinstitute/CaImAn-MATLAB

And some functions from GitHub/restrepd/drgMaster

All the data was deposited in GigaDB.

## Running the code to generate the figures in Ma et al. 2020

The image files were processed with CaImAn by running drgCalmAn_script.m

The output from experiments acquired with dropcspm_hf.m was then processed using drgCaImAn_batch_dropc.m that creates per trial analysis and saves an intermediary file with a suffix pre_per.mat used for the rest of the analysis. drgCaImAn_batch_dropc.m needs as input the INTAN .rhd files, the olfactometer .mat files and the CaImAn-processed .mat files. The location of these files is saved in a drgCaImAn_dropc_choices file such as drgCaImAn_dropc_choices_GRIN1_CA1.m.

drgCalmAn_script.m Finds components (ROIs) using CaImAn code. It ouputs a file called file_name_CalmAn.mat. It is important to enter the CaImAn parameters using options = CNMFSetParms

drgCaImAn_batch_dropc.m This code takes as input the output of drgCalmAn_script (file_name_CalmAn.mat), an output file from the INTAN RHD2000 board (.rhd) that has digital output sent by the olfactometer and lick recordings and a .mat file from the olfactometer with metadata on each trial. The code generates and intermediary file called file_name_CalmAn_batch_pre_per.mat.

readdropcspm_hf.m computes percent correct behavior as a function of trial number using a window of 20 trials. It takes as input several olfactometer output files (please name files 1-textID_spm.mat, 2-textID_spm.mat, etc).

Most of the code is not well commented (sorry), but we provide a well-documented file (drgCaImAnLDAtimecourse.m)that has been used outside our group (yes!) by David Protter and Zoe Donaldson to implement our LDA algorithm.

To run this reasonably well-commented example run LDA_demo.m that loads data and runs drgCaImAnLDAtimecourse.m.

Below are per figure instructions on how we generated each panel. For this code you often have to enter a choices file that has all the information (where the files are found, their names, the time windows to be used, etc). We have also provided the choice files. You will have to modify them to have the code run in your computer.

### Fig. 1c

For Fig.1c (and Supplementary Fig. 4e) use drgCaImAn_curate_components.m using:
5-cerebellum-mmg04-spm_180910_180910_161628.rhd
20180910_mmPVG04_Cerebellum-5-Registered_CalmAn.mat
5-mmPVG04-cerebellum-spm20180910T161640spm.mat

### Fig. 1d

For Fig. 1d use drgCaImAn_rainbow_mod.m with the same input files.

### Figs. 2a,b,c,d

For Figs. 2a,b,c,d run drgCaImAnBatchPerSessionPerTrial.m using the choices file drgCaImAnChoices_20180515_mmPVG02_Cerebellum.m

### Fig. 2e

For Fig. 2e use drgCaImAnBatchPerSessionPerROI.m to generate a violin plot using the choices file drgCaImAnChoices_20180515_mmPVG02_Cerebellum.m

### Fig. 2f,g

For Figs. 2f,g first process data with drgCaImAnBatchPerSession.m

Using the following choices files:

drgCaImAnChoicesDiego20180910_mmPVG04_Cerebellum.m
that yields 20180910_mmPVG04_Cerebellum_bpsout.mat

drgCaImAnChoices_20181017_19mmPVG05_Cerebellum.m
yields 20181017_19_mmPVG05_Cerebellum_out_bpsout.mat

drgCaImAnChoices_20180515_18_mmPVG02_Cerebellum.m
yields 20180515_18_mmPVG02_Cerebellum_out_bpsout.mat

drgCaImAnChoicesDiego20180702_05_mmG7f09_Cerebellum.m
yields 20180702_05_mmG7f09-Cerebellum_out_bpsout.mat

Then run drgCaImAn_average_dFF_multi_session.m

### Fig. 3b

For Fig. 3b. Example of principal component analysis generated with drgCaImAnBatchPerSessionReversalPerTrialPCA.m using the choices file drgCaImAnChoicesDiego20180910_mmPVG04_CerebellumPCA.m

### Fig. 3c

For Fig. 3c run drgCaImAnBatchPerSessionReversalPerTrialLDA.m
using drgCaImAnChoicesDiego20180910_mmPVG04_CerebellumLDA.m

### Fig. 3d

For Fig. 3d first run drgCaImAnBatchPerSessionReversalPerTrialLDA.m with the following choice files:

drgCaImAnChoicesDiego20180910_mmPVG04_CerebellumLDA.m
yields 20180910_mmPVG04_Cerebellum_lda.mat

drgCaImAnChoices_20180515_18_mmPVG02_CerebellumLDA.m
yields 20180515_18_mmPVG02_Cerebellum_out_lda_new_sum_lda.mat.

drgCaImAnChoicesDiego20180702_05_mmG7f09_CerebellumLDA.m generated
yields 20180702_05_mmG7f09-Cerebellum_lda_sum_lda

drgCaImAnChoices_20181017_19mmPVG05_CerebellumLDA.m generated
20181017_19mmPVG05_Cerebellum_out_lda_sum_lda.mat

drgCaImAnBatchPerSessionReversalPerTrialLDA.m performs a linear decoding analysis

Then run summary_lda_recalc.m to generate Fig. 3d

### Fig. 3e

Run drgCaImAnBatchPerSessionPerTrialLDASubsample.m to perform LDA analysis with subsets of ROIs using the choice file drgCaImAnChoicesDiego20180910_mmPVG04_CerebellumLDA.m.

### Fig. 3f

Run drgCaImAnBatchPerSessionPerTrialLDASubsample with the following choices files:

drgCaImAnChoicesDiego20180910_mmPVG04_CerebellumLDA.m
yielding 20180910_mmPVG04_Cerebellum_new_out_ldasub.mat

drgCaImAnChoices_20180515_18_mmPVG02_CerebellumLDA.m yielding 20180515_18_mmPVG02_Cerebellum_out_lda_new_ldasub.mat

drgCaImAnChoicesDiego20180702_05_mmG7f09_CerebellumLDA.m
yielding 20180702_05_mmG7f09-Cerebellum_lda_sum_ldasub.mat

drgCaImAnChoices_20181017_19mmPVG05_CerebellumLDA.m
yielding 20181017_19mmPVG05_Cerebellum_out_lda_sum_ldasub

Then run summary_lda_subset.m

### Figs. 4a,b,c

Use drgCaImAnBatchPerSessionReversalPerTrial.m
with choice file drgCaImAnChoicesDiego20180917and19_mmPVG04_Cerebellum.m

### Fig. 4d

Run All_mice_reversal_recalc.m that uses data processed by drgCaImAnBatchPerSessionReversalPerTrial.m

### Fig. 4e

Use drgCaImAnBatchPerSessionReversalPerTrialLDA.m with choice file drgCaImAnChoicesDiego20180917and19_mmPVG04_CerebellumLDA.m

### Fig. 4f

First use drgCaImAnBatchPerSessionReversalPerTrialLDA.m with the following choice files:

drgCaImAnChoicesDiego20180917and19_mmPVG04_CerebellumLDA.m
yields 20180917and19_mmPVG04_Cerebellum_lda.mat

drgCaImAnChoicesDiego20180608_11_15_mmG7f09_CerebellumLDA.m yields
20180608_mmG7f09_Cerebellum_lda.mat

drgCaImAnChoicesDiego20180419and23_mmG06_cerebellumLDA.m
yields 20180419and23_mmG06_cerebellum_lda.mat

Then run All_mice_lda_reversal_recalc.m

### Figs. 5a-c

Run drgCaImAnBatchPerSessionEventsPerTriallickvsdFF.m with choices file
drgCaImAnChoicesDiego20180917_mmPVG04_Cerebellum_licks.m

### Fig. 5d

For Hit, CR and FA run drgCaImAnBatchPerSessionEventsPerTriallickvsdFF.m with choices file drgCaImAnChoicesDiego20180917_mmPVG04_Cerebellum_licks.m. The trials shown are fund in the time periods 60-90, 200-230, 270-300 seconds

For Miss run drgCaImAnBatchPerSessionEventsPerTriallickvsdFF.m with choices file drgCaImAnChoices_20181017_mmPVG05_Cerebellum.m 1320-1350 seconds.

### Fig. 5e,f

Run drgCaImAnBatchPerSessionEventsPerTriallickvsdFF.m with choices file
drgCaImAnChoicesDiego20180917_mmPVG04_Cerebellum_licks.m

### Fig. 5g,h

First run drgCaImAnBatchPerSessionEventsPerTriallickvsdFF.m with the following choices files�

drgCaImAnChoicesDiego20180917_mmPVG04_Cerebellum_licks.m
drgCaImAnChoicesDiego20180910_mmPVG04_Cerebellum.m
drgCaImAnChoices_20180515_18_mmPVG02_Cerebellum.m
drgCaImAnChoicesDiego20180702_05_mmG7f09_CerebellumLDA.m
drgCaImAnChoices_20181017_19mmPVG05_Cerebellum.m
drgCaImAnChoicesDiego20180419_mmG06_cerebellum_LDA_events.m

Then run summary_derivativesdFFlick_recalc.m

### Figs 6c,d

First run drgCaImAnBatchPerSessionEventsPerTriallickvsdFFOptFlow.m with the following choices files:

drgCaImAnChoicesDiego20180910_mmPVG04_CerebellumLDAOpFlow.m
drgCaImAnChoices_20180518_mmPVG02_CerebellumLDAOpFlow.m
drgCaImAnChoicesDiego20180702_05_mmG7f09_CerebellumLDAOptFlow.m
drgCaImAnChoices_20181017_mmPVG05_CerebellumLDAOpFlow.m
drgCaImAnChoicesDiego20180917_mmPVG04_Cerebellum_for_PCAOpFlow.m
drgCaImAnChoicesDiego20180419_mmG06_cerebellum_for_PCA_OptFlow.m

Then, using the output from drgCaImAnBatchPerSessionEventsPerTriallickvsdFFOptFlow.m run drgCaImAn_glm.m

20180910_mmPVG04_Cerebellumopflow_slopestest.mat to yield
20180910_mmPVG04_Cerebellumopflow_slopestest_glmtest.mat

2018018_mmPVG02_Cerebellum_OpFlowout_slopestest.mat to yield
2018018_mmPVG02_Cerebellum_OpFlowout_slopestest_glmtest.mat

20180702_05_mmG7f09-Cerebellumbatch_slopestest.mat to yield
20180702_05_mmG7f09-Cerebellumbatch_slopestest_glmtest.mat

20181017_mmPVG05_Cerebellum_out_slopestest.mat to yield
20181017_mmPVG05_Cerebellum_out_slopestest_glmtest

20180917_mmPVG04_Cerebellum_PCA_slopestest.mat to yield
20180917_mmPVG04_Cerebellum_PCA_slopestest_glmtest.mat

20180419_mmG06_cerebellumPCAevents_slopestest.mat to yield
20180419_mmG06_cerebellumPCAevents_slopestest_glmtest.mat

summary_glm_recalc.m uses the output from drgCaImAn_glm.m to calculate Figs. 6c,d

### Fig. 7aii,bii

Percent_correct_dreadds.m calculates the statistics for Figs. 7aii and 7bii

### Fig. 7c,d

Run drgCalmAn_batch_dropc_no_microscope run with
drgCaImAn_dropc_choices_PVhM4Di_081519.m
Generates *_batch_licks.mat


Then run drgCaImAnBatchMultiSessionLicksv3.m which uses the output file from the last run.

If you would like to get information on how to generate the Supplementary Figures we will be glad to add those to readme.md

## For processing Fabio's CA1 data

drgCaImAn_batch_dropc_fsdz loads dF/F from csv files generated through ImageJ. It needs as input a drgCaImAn_dropc_choicesfile such as drgCaImAn_dropc_choices_Grin3_fsds1_2_20210302.m
If the imaging file is .csv it assumes that it is ROIs processed by Fabio in ImageJ, if it is .mat it assumes it is the output from EXTRACT.

drgCaImAn_pre_per_to_pydec reads the pre_per file and saves data for Kording's lab python neural_decoding

drgCaImAnBatchPerSessionReversalPerTrialLDA_fsdz reads the pre_per file and does an LDA analysis for the go-no go spm

drgCaImAn_pre_per_to_LDA_fsdz does a decoding analysis for Splus vs. Sminus

drgCaImAn_batch_pre_per_to_LDA_fsdz performs a batch decoding analysis for Splus vs. Sminus. It needs as an input a choices file such as drgCaImAn_LDAfsdz_choices_05242021.m
Discriminant is calculated per group.

drgCaImAn_batch_analysis_pre_per_to_LDA_fsdz dispays the results of drgCaImAn_batch_pre_per_to_LDA_fsdz

drgCaImAn_batch_pre_per_to_KLD_fsdz with input from a file such as drgCaImAn_KLDfsdz_choicesMing_12042021 shows pseucolor plots of thersholded dF/F, timecourse of the overall thresholded dF/F and KL divergence.

## To processing Ming's go-no go CA1 data

## Authors

* Diego Restrepo (diego.restrepo@cuanschutz.edu) with help from the authors of the manuscript.

## License

This project is licensed under the GNU General Public License v3.0
 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thanks to Jesse Gilmer for help with the code for dimensionality
