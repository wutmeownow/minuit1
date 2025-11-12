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
Best fit parameters and uncertainties:
- A      : 53.9973 +- 0.508397
- mux    : 3.51647 +- 0.00543451
- sigmax : 0.698328 +- 0.00484579
- muy    : 1.90636 +- 0.015083
- sigmay : 1.37966 +- 0.0124921
- Abkg   : 0.246314 +- 0.00228712

- Reduced chi sqr = 1.12
- p value = 4.29e-07

Note: For the lego plot for data minus fitted background, I set negative bins to zero to more easily compare to the data histogram.

I estimated the signal count as 32100 ± 200 by summing all the bins (including negative ones) in the histogram of data minus the fitted background. My rational is that while some bins may overcount the number of signal events, the bins that are negative will lead to an overall cancellation of these false signal counts. To find the error I used error propagation on the number of counts in each bin and then summed these in quadrature and took the square root. The error of the count in each bin in the final histogram should be the sqrt of the quadrature sum of root(N) of the data histogram bin and the error in the Abkg fitted parameter scaled by the value of the background model at that bin center.
-----
