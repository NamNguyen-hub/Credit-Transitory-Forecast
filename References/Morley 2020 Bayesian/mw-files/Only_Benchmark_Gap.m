% Benjamin Wong
% Monash University
% July 2019
%
% The script only calls out the files needed to just get the output gap, and nothing more
% If you want the full set of results, please use MAIN.m
%
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

% Estimate Gap only
lambda = fminsearch(@(lambda0) BN_BVAR_oos_RMSE(y{2},p,lambda0,target_variable(2)),0,options.optimisation);
[BN_cycle,Info_decom,shock_decom,FEVD,Std_err] = BN_BVAR(y{2},p,lambda,target_variable(2),'Decomposition','Standard_Errors');


%%
figure
h1=NBERbc(dates,BN_cycle(:,target_variable(2)),{'-'},3,{'r'});
hold on
plot(dates,BN_cycle(:,target_variable(2))+1.645*Std_err(target_variable(2)),'--k','LineWidth',1);
hold on
plot(dates,BN_cycle(:,target_variable(2))-1.645*Std_err(target_variable(2)),'--k','LineWidth',1);
hold on
plot([dates(1) dates(end)],zeros(2,1),'-k','LineWidth',2)
ylim([-8 6])
set(gca,'FontSize',16)

toc

