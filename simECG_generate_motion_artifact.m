function [simNoise] = simECG_generate_motion_artifact(successProb, sigmaBernoGauss, sigmaARinput, modalityFlag)
% INPUTS:
% - successProb : the probability of success, i.e., spikes 

% - sigmaBernoGauss : standard deviation of guassian noise that is to be multiplied
% by the bernoulli process

% - sigmaARinput : standard deviation of the gaussian noise used in the
% AR/ARIMA model input
% - modalityFlag 0 or 1 (switch between AR and ARIMA model; integrator)
% 0 - Holter recording  1 - Thumb-ECG 

%--------
% Hesam Halvaei, Lund University
%--------
% - Load the AR(4) coefficient estimated based muscle noise.
% - AR_MN is a library of AR(4) coefficients obtained from 25 3-lead muscle
% noise recordings.

load("DATA_AR_MN_Dictionary.mat");
% The coefficient are selected randomly.
ar_coeffs = AR_MN(randi([1, size(AR_MN, 1)]), :, randi([1, size(AR_MN, 3)]));

% N is the length of simulated noise, corresponds to a 30 s noise based on
% sampling frequecy of 1 Hz
N = 30000; 
[bernogauss] = simECG_Bernoulli_Gaussian(N, successProb, sigmaBernoGauss);
[bernogauss_conv] = simECG_Bernoulli_Gaussian_convolution(N, bernogauss);

white_noise = normrnd(0, sigmaARinput, [N, 1]);
ar_input = white_noise  + bernogauss_conv;
ar_input_decimated = decimate(ar_input, 5);
ar_model = idpoly(ar_coeffs, 'integratenoise', modalityFlag);
simNoise_decimated= sim(ar_model, ar_input_decimated);
simNoise = interp(simNoise_decimated, 5);
end

