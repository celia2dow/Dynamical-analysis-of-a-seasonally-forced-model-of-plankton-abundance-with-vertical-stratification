# Dynamical-analysis-of-a-seasonally-forced-model-of-plankton-abundance-with-vertical-stratification
Repository of code used to generate Figures 3 and 4 in the same-named article for the StAMPS 2026 Magazine (Issue 2). The driving script is `drivingResultsScript` which calls all other functions.

- To initialise parameters: run the 1<sup>st</sup> code-chunk of `drivingResultsScript`.
- To produce Figure 3: run the 2<sup>nd</sup> code-chunk.
- To produce Figure 4(a-c)(i): run the 3<sup>rd</sup> code-chunk.
- To produce Figure 4(a-c)(ii) and (d): run the 4<sup>th</sup> code-chunk.

Some adjustments may need to be made to the code and the output to account for non-converging `ode45` time-series solutions and erroneous `fsolve` equilibria.

This work is an excerpt of my masters thesis entitled *Investigating the distribution of plankton using a seasonally forced model with ocean dynamics* (2023) which is included herein.

Please contact me if any questions arise.
