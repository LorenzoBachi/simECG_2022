function [simuMA] = simECG_generate_motion_artifact(successProb, sigmaBernoGauss, modalityFlag, noiseRMS, N)
% INPUTS:
% - successProb : the probability of success, i.e., spikes 

% - sigmaBernoGauss : standard deviation of guassian noise that is to be multiplied
% by the bernoulli process

% - modalityFlag 0 or 1 (switch between AR and ARIMA model; integrator)
% 0 - Holter recording  1 - Thumb-ECG 

% -  N is the length of simulated noise,

%--------
% Hesam Halvaei, Lund University
%--------
% - Load the AR(4) coefficient estimated based muscle noise.
% - AR_MN is a library of AR(4) coefficients obtained from 25 3-lead muscle
% noise recordings.

load("DATA_AR_MN_Dictionary.mat"); %parameters were calculated with signals in uV

% The coefficient are selected randomly.
if modalityFlag %3 VCG leads or 1-lead
    ar_coeffs = squeeze(AR_MN(randi([1, size(AR_MN, 1)]), :, randi([1, size(AR_MN, 3)])));
    L = 1;
else
    ar_coeffs = squeeze(AR_MN(randi([1, size(AR_MN, 1)]), :, :));
    L = 3;
end


[bernogauss] = simECG_Bernoulli_Gaussian(N, successProb, sigmaBernoGauss);
[bernogauss_conv] = simECG_Bernoulli_Gaussian_convolution(N, bernogauss);

simuMA = zeros(L,N);
sigmaARinput = noiseRMS*0.3;

for Li = 1:L
    white_noise = normrnd(0, sigmaARinput, [N, 1]); %in mV
    ar_input = (white_noise  + bernogauss_conv).*1e3; %in uV
    ar_input_decimated = decimate(ar_input, 5);
    ar_model = idpoly(ar_coeffs(:,Li)', 'integratenoise', modalityFlag);
    simNoise_decimated = sim(ar_model, ar_input_decimated);
    simuMA(Li,:) = interp(simNoise_decimated, 5).*1e-3; %in mV
end

if ~modalityFlag % Transform to the 15 leads
    simuMA_8 = leadcalc(simuMA,'stan');% V1,V2,V3,V4,V5,V6,I,II
    simuMA_12 = leadcalc(simuMA_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
    
    simuMA_15 = vertcat(simuMA_8(7,:),simuMA_8(8,:),simuMA_12(12,:),...
        -simuMA_12(9,:),-simuMA_12(7,:),simuMA_12(11,:),simuMA_8(1:6,:),simuMA);
    
    simuMA = simuMA_15;
end
end

