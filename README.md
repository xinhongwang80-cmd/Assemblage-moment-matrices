# Assemblage-moment-matrices

To run the codes involving solving semi-definite programs, one has to use YALMIP as the modelling language and any solver compatible with YALMIP (e.g., SDPA, SDPT3, SDPNAL, etc.).

First, download all the codes and add them to the working path.

Filenames starts with "AB_C..." are examples.

For example, "AB_C_AMM_L_Mermin_SR_lowerbound.m" computes the lower bound of steering robustness for a given Mermin inequality violation, through AMM approach. In the filename "L" means "local" correlations.

Other example, "AB_C_AMM_NS_Mermin_SR_lowerbound.m" computes the lower bound of steering robustness for a given Mermin inequality violation, through AMM approach. In the filename "NS" means "no-signaling" correlations.
