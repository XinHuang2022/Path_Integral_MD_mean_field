%
% Tester for PIMC_BB_sampling_run.m
%

filename0='testrun_PIMC_BB';
lulog0='1';
datadir0='./';
dim0 = '3';
Ne0 = '6';
beta0 = '1';
lambda0 = '0';
dt0 = '0.025';
Mxmax0 = '12';
dMx0 = '0.5';
Nrep0 = '8';
Ncore0 = '4';
TestCase0 = "V1";

exitcode = PIMC_BB_sampling_run_improve(filename0, lulog0, datadir0, ...
    dim0, Ne0, beta0, lambda0, dt0, Mxmax0, dMx0, Nrep0, Ncore0, TestCase0)

%%%
% dim0:      the dimension of the space.
% Ne0:       the number of electrons.
% beta0:     the inverse temperature 1 / ( k_B * T ).
% lambda0:   the relative force parameter of Coulomb interactions. 
%            For the fermi gas case under Harmonic Oscillator potential 
%            without Coulomb interactions, lambda will be automatically 
%            set to 0.
% dt0:       the time step for the integration in the exponent of the 
%            Feymann-Kac representation.
% 2^Mxmax0:  the sample size for each independent estimator of partition 
%            function Z and mean-field h.
% dMx0:      the x-axis distance between two points in the figures for the 
%            convergence of estimators with increasing sample size.
% 2^Nrep0:   the total number of independent replicas of 
%            partition function Z and mean-field h, each replica generates
%            and stores 2^Mxmax0 samples for Monte Carlo integral.
% Ncore0:    the number of parallel workers for evaluating all the 
%            2^Nrep0 replicas.
% TestCase0: Case V1 denotes Harmonic Oscillator potential, Case V2 for 
%            Coulomb interaction under a harmonic trap.
