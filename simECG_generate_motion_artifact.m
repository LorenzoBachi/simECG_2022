function [simuMA] = simECG_generate_motion_artifact(ecgLength,ecgParameters, noiseRMS)
% INPUTS:
% - successProb : the probability of success, i.e., spikes 

% - sigmaBernoGauss : standard deviation of guassian noise that is to be multiplied
% by the bernoulli process

% - modalityFlag 0 or 1 (switch between AR and ARIMA model; integrator)
% 0 - Holter recording  1 - Thumb-ECG 

% -  N is the length of simulated noise,

%--------
% Hesam Halvaei, Lund University
% Update: Cristina Perez, University of Zaragoza, 09/2022
%--------
% - Load the AR(4) coefficient estimated based muscle noise.
% - AR_MN is a library of AR(4) coefficients obtained from 25 3-lead muscle
% noise recordings.

if  ~isfield(ecgParameters,'sigmav') %no previous simulated muscular noise
    ecgParameters.sigmav = (0.03*noiseRMS)*1e3; %in uV;
end

if ecgParameters.MA_Flag %Thumb-ECG
    L = 1;
else %3 VCG leads or 1-lead
    L = 3;
end

for Li = 1:L
    bernogauss = simECG_Bernoulli_Gaussian(ecgLength, ecgParameters.MA_Prob, ecgParameters.sigmav);
    bernogauss_conv(Li,:) = simECG_Bernoulli_Gaussian_convolution(ecgLength, bernogauss)'; %in uV
end

simuMA = bernogauss_conv.*1e-3; %in mV

if ~ecgParameters.MA_Flag % Transform to the 15 leads
    simuMA_8 = leadcalc(simuMA,'stan');% V1,V2,V3,V4,V5,V6,I,II
    simuMA_12 = leadcalc(simuMA_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
    
    simuMA_15 = vertcat(simuMA_8(7,:),simuMA_8(8,:),simuMA_12(12,:),...
        -simuMA_12(9,:),-simuMA_12(7,:),simuMA_12(11,:),simuMA_8(1:6,:),simuMA);
    
    simuMA = simuMA_15;
end
end

