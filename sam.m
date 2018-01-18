close all;
clear all;

sam_constants;
sam_parameters;

seed=1;
[Fe Fee] = mjlstd(lambda,J,T,K,epsilon,N,P,As,Bs,Cs,Ds,F_opt,X_opt,seed)
