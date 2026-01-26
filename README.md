# Dynamical-analysis-of-a-seasonally-forced-model-of-plankton-abundance-with-vertical-stratification
Repository of code used to generate Figures 3 and 4 in the same-named article for the StAMPS Summer 2026 Magazine issue. The driving script is `drivingResultsScript` which calls all other functions.

- To initialise parameters: run the 1$^{st}$ code-chunk of `drivingResultsScript`.
- To produce Figure 3: run the 2$^{nd}$ code-chunk.
- To produce Figure 4(a-c)(i): run the 3$^{rd}$ code-chunk.
- To produce Figure 4(a-c)(ii) and (d): run the 4$^{th}$ code-chunk.

Some adjustments may need to be made to the code and the output to account for nonconverging `ode45` time-series solutions and erroneous `fsolve` equilibria.

Please contact me if any questions arise.
