James Morley and Benjamin Wong, "Estimating and Accounting for the
Output Gap with Large Bayesian Vector Autoregressions", Journal of
Applied Econometrics, Vol. 35, No. 1, 2020, pp. 1-18.

All files are zipped in mw-files.zip. Text files are in DOS format,
but there are also binary files, so not-Windows users should exercise
cuation.

MAIN.m will generate all the Figures in the paper. If you just want
the output gap in the benchmark 23 variable VAR, please run
Only_Benchmark_Gap.m

The entire Main.m script runs on Benjamin Wong's desktop in
approximately 15 mins in July 2019.

Figures A1 and A2 of the online appendix are also generated to
demonstrate how to do structural analysis with the framework. The oil
price and monetary policy shocks are identified using widely used and
standard identification restrictions. The individually named scripts
are meant to correspond with each figure in the paper. They should
each individually work as long as the benchmark output gap is
generated (Figure 2) and the 3 datasets of different coverage are
generated when setup_dataset.m is run.

The three datasets are in the y cell which will be generated as long
as the setup_data.m file is run.

y{1} contains the variables for the 8 variable VAR.

y{2} contains the variables for the benchmark 23 variable VAR.

y{3} contains the variables for the 138 variable VAR.

All data are sourced from FRED. The code will call out the relevant
variables within the spreadsheet.

Most of the analysis will be with the data in y{2} as this is the
benchmark set of variables.

Generating Figure_9 will require MATLAB's parrallel computing toolbox.
If you do not require Figure_9, you can just comment out Figure_9 in
Main.m. To do the parrallel computing, you might have to change the
number of processor cores. If you do not have parrallel computing, you
can comment out the sections calling and shutting down the multiple
cores (i.e. parpool and delete(gcp)), and change the parfor loop into
a for loop. While this will still do the rolling window exercise, you
can expect the calculations to slow down by a factor of 4 or 8
depending on how many cores are being used. 

If you use the code, please cite the paper.

If you find errors, please get in touch with us.

We cannot be held responsible for any use and mis-use of our code.

James Morley
james.morley [AT] sydney.edu.au
University of Sydney

Benjamin Wong
benjamin.wong [AT] monash.edu
Monash University

July 2019
