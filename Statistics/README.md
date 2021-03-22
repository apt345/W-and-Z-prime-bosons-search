# FYS5555 Statistics

Prepared by Magnar Bugge

This section allows the computation of "cut-and-count" ("single-bin") cross section limits on a new signal
process and production of standard limit plots including bands for the 1- and 2-sigma variations of the
limit around the expected (median) value.

Inputs are read from the file "inputs.txt", which is currently set up to reproduce the electron channel
limits (Fig. 2, top left) from https://arxiv.org/abs/1407.7494. The format is

[Number of mass points]
[W' mass] [mT threshold]
[signal efficiency] [signal eff. uncertainty] [number of predicted signal events] [background] [background uncertainty] [observed data]
[W' mass] [mT threshold]
[signal efficiency] [signal eff. uncertainty] [number of predicted signal events] [background] [background uncertainty] [observed data]
[W' mass] [mT threshold]
[signal efficiency] [signal eff. uncertainty] [number of predicted signal events] [background] [background uncertainty] [observed data]
...
...
...

The plot from the above publication can be reproduced with the commands:
- `root -l -q FYS5555/Statistics/runLimits.cpp`
- `root -l FYS5555/Statistics/limitPlot.cpp`
- `root -l FYS5555/Statistics/runLimits_combined.cpp`

The code depends on ROOT (which you have already installed in your Anaconda environment).
