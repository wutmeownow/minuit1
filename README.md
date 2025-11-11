# minuit1

- expFit.cpp: example using the more general fitting interface in minuit
- expFit.ipynb: equivalent example using lmfit
- SimultaneousExps(lm).ipynb: generation of histograms with correlated signals for simultaneous fit exercise
- rootExample.cpp: just another example of using ROOT classes in a C++ program
- *.root files: various input histograms for fitting exercises

-----

Team member names and computing IDs: Michael Lowry (mjl9ka)

-----

Exercise 1 comments:
--
The best fit of the Gumbel distribution has a smaller reduced chi square of 2.66 vs the 2.75 of the double gaussian fit. The p-value of the gumbel distribution fit is 1.5e-16 and 3.2e-17 for the double gaussian fit. The Gumbel distribution is favored over the double gaussian model.

-----
Exercise 2 comments:
--
I extracted the mean and sigma of the signal as 74.6 ± 0.5 and 4.60 ± 0.44. Both agree with the true values within one sigma. The number of degrees of freedom should be 2*50 bins - 8 parameters = 92, which results in a reduced chi squared of 0.83 and a chi squared probability of 0.88. I think this is a good fit to the data since the signal parameters agree with the true values and the chi squared probability is quite large. Since it is so large the fit may be a little too good and overfitting.

-----

Exercise 3 comments:
--


-----
