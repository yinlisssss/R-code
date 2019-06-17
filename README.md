# R-code
R code for reproducing the results for manuscript submitted to BMC cancer.

The code covers key steps in our manuscript. 

Step1: Data download

Step2: Gene selection

Step3: Lasso regression

Step4: Roc and survival analysis

Step5: Stratification analysis

Step6: Nomogram

To recover GEO test and validation set data, run this in shell or Macos terminal.

```
cat external_test_part* > external_test.RData
cat external_validation_part* > external_validation.RData
```
