% Benjamin Wong
% Monash University
% July 2019
%
% The set of files replicates Morley & Wong,
% "Estimating and Accounting for the Output Gap with Large Bayesian Vector
% Autoregression" Journal of Applied Econometrics, forthcoming
%%
clear all
clc

addpath(genpath('_func'))
addpath('datasets')
addpath('other_scripts')

%% Preliminaries

% VAR lag order
p = 4;
% 1959Q3 to 2016Q4
dates = (1959.5:0.25:2016.75)';

setup_dataset
tic

%% Generate Figures

Figure_1
Figure_2
Figure_3
Figure_4
Figure_5
Figure_6
Figure_7
Figure_8
Figure_9

%% Appendix Figures

% We reproduce Figures A1 and A2 to just demonstrate how to do structural anaysis.
% The identification is a standard Cholesky Decomposition
% Oil Price shocks identified by ordering oil price first
% Monetary policy shock identified by ordering Federal Funds Rate after
% slow moving but before fast moving variables
% Other identification schemes are definitely possible as long as one
% alters the A0 matrix

Figure_A1_A2

disp('Total Time Taken')
toc