# masters-thesis
R-code for computing confidence intervals for the crash intensity and crash probability using the profile- and modified profile likelihood.

# Functions
Functions for computing confidence intervals for the (near) crash-intensity and the difference in (near) crash probability can be found in this repository. The two functions for near-crash intensity can be found in pl_lambda.R and mpl_lambda.R, and the two functions for crash probability are pl_diffpi.R and mpl_diffpi.R. If the function name includes pl it implies that profile likelihood was used, and if it contains mpl it implies the modified profile likelihood was used.

# Formatting input 
The functions are all take the following elements as arguments:

  - A dataframe with named columns "Date", and at least one of the following "TTCmin", "T2min", "PET" measurements. For the diffpi functions, you have to enter two dataframes as input.
  - A numeric threshold corresponding to the threshold of the GPD. The input name is either thresh, uNE or uDK. 
  - A numeric quantile q or x, which corresponds to the desired value for which one wants to compute the intensity or probability for S<q.
  - A string corresponding to the desired SMoS which one wants to use to estimate the (near) crash intensity or probability. 

# Date
The "Date" columns of the dataframes correspond to the time points of the interactions. They have to include date and time in the following format:

y-m-d h:m:s

