# Assemblage-moment-matrices

This repository contains MATLAB code for computing lower bounds on tripartite steering robustness using the assemblage moment matrix (AMM) approach. We consider the steering scenario, in which Alice and Bob jointly steer Charlie.

Running the codes requires YALMIP together with a compatible SDP solver, such as SDPA, SDPT3, or SDPNAL.

Before running the examples, download all the files and add the relevant folders to the MATLAB search path.

Files whose names start with AB_C are examples.

For example, "AB_C_AMM_L_Mermin_SR_lowerbound.m" computes a lower bound on steering robustness for a given Mermin inequality violation using the AMM approach. In the filename, "L" indicates that the joint response of Alice and Bob in the LHS model is restricted to the local set.
Similarly, "AB_C_AMM_NS_Mermin_SR_lowerbound.m" uses the AMM approach to compute a lower bound on steering robustness for a given Mermin inequality violation. Here, "NS" indicates that the joint response of Alice and Bob is allowed to belong to the no-signalling set.

