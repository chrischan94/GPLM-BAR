# GPLM-BAR

In this GitHub respository, we provide the sample code to our novel approch for simulataneous variable selection and estimation under the context of Generalized Partly Linear Models (GPLM). In order to replicate our results, you would need to install the 'BrokenAdaptiveRidge' package. To do this, run the following:

```{R}
install.packages("devtools")
install.packages("Cyclops")
install.packages("BrokenAdaptiveRidge")
```
The Cyclops package
We provide the sample code for the three scenarios (Scenario 1, Scenario 2 and Scenario 4) in our manuscript. For scenario 1, where the strong signals in beta are considered, please refer to LPLMNEWn600p300.R. For scenario 2, where a mixture of strong and weak signals in beta are considered, please refer to LPLMWEAKn600p300.R.
