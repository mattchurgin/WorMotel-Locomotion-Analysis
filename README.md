# WorMotel-Locomotion-Analysis
This repository includes the files required to quantify locomotion behavior state for C. elegans swimming in liquid in a WorMotel as described in the following paper: http://www.jneurosci.org/content/early/2017/07/11/JNEUROSCI.2636-16.2017

The following files are included.

1) WorMotel_GUI: Image analysis GUI.  Use this to process images acquired in the format 'Image0001.bmp', 'Image0002.bmp',...
A readme file is also included.

2) createMattFit2.m: This function is used to fit histogram data from a single worm's activity trace.  It requires two inputs: X and Y data for the histogram of activity values.  This fit function fits activity histogram data to the sum of two exponentials and a single Gaussian curve.  Dwelling is calculated as the area under the two exponential terms and roaming is calculated as the area under the Gaussian.  Quiescence is calculated as the total time spent at zero activity (and is therefore not included in the fit).
This function is called in Analysis_Example.m (below)

3) Analysis_Example.m: This matlab script shows how to analyze data.  An example data set is provided below.

4) An example data set of 48 worms imaged for 16 hours in the wormotel.
