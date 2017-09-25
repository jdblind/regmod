* ./optim_routine_exsarcroi.R implements the main optimization body for the proposed regularization model.

* Root script ./synthetic_data.R generates synthetic 3D data and carries out a performance analysis for the optimsation (regularized model fitting) process for varying regularization values and SNR levels.

* Root script ./quant_cohort.R demonstrates how qantitation was carried out for the sarcoma cohort. Note that it has been edited to hide patient information, which we cannot share at this point, and will not run because it loads and uses patient data. The code is shared to demonstrate replicability of results.

* Folders scripts/, funs/ and utils/ contain functions and R code modules used for quantitation and post-analysis.

* Folder Cpack/ contains useful C-implementations. 

* Scripts ./scripts/_plot.R, ./scripts/_plot_figures.R and ./scripts/_plot_case_studies.R  provide examples of figures like those used in the paper.

* Note that script ./scripts/script_analysis_single_cases.R will not run either because it loads patient data, which we cannot share at this point. However the code is shared to demonstrate replicability of results.

* NOTE:
Tarballs and zipballs 
Cpack/clapack-3.2.1-CMAKE.tgz
Cpack/f2c.tar
Cpack/levmar-2.6.tgz
Cpack/libf2c.zip
may need to be unzipped/untarred in same directory.