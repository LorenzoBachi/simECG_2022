function [simuMN_noise_15] = simECG_Muscular_Noise(ecgLength)
%function [simuMN_noise] = simECG_Muscular_Noise(ecgLength, noiseType, noiseRMS, ESTParameters)
% simuMN_noise = simECG_Muscular_Noise() returns a simulated muscular noise
% signal.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 07/2022

% 1) Load dictionary AR(p) model and select the parameters that model the quasy-stationay part of
% the simulated MN signal
load('AR_MN_Dictionary.mat')
ARp = squeeze(AR_MN(randi([1,25]),:,:));
fs = 1000;

% 2) AP(1) model to calculate the variance of the MN signal
clc
N = ecgLength;
W = 30;
M = ceil((N/fs)/W) +1; %time-varying 30seconds

scal = zeros(1,N);
v1All = [];
out_AP1All = [];
v1 = [];
out_AP1 = [];

%--> 1) Select the value of the pole and noise distribution
p = rand(1)*(0.9999-0.999) + 0.999;
b = 1;
a = [1 -p]; %according to me

%--> 2) Apply the filter and Spectra analysis of each filter
v1 = randn(3, ecgLength + 50000); %Frank leads
out_AP1 = filter(b,a,v1')';
v1(:,1:50000) = [];
out_AP1(:,1:50000) = [];

v1All = [v1All v1];
out_AP1 = (out_AP1 - 0)./sqrt(var(v1,[],2)/(1-p^2));%Normalize
scal(1) = rand(1)*(4-1) + 1;
out_AP1 = (scal(1) + out_AP1.^2); %variance
out_AP1All = [out_AP1All out_AP1];


% 3)AP(n)model to obtain the desired simulated muscular noise signal
%-->1) Resample to 200Hz
out_AP1All_200 = [];
for ii = 1:3
    out_AP1All_200(ii,:) = resample(out_AP1All(ii,:), 200, 1000);
end

%-->2) Apply AR(n) filter
v2_200 = [];
v2_200 = randn(size(out_AP1All_200,1),size(out_AP1All_200,2));
v2_200 = v2_200.*out_AP1All_200;

simuMN_noise_200 = [];
for Li = 1:3
    simuMN_noise_200(Li,:) = filter(1,ARp(:,Li),v2_200(Li,:)')';
end

%-->3) Resample to 1000Hz
simuMN_noise = [];
v2 = [];
for ii = 1:3
    simuMN_noise(ii,:) = resample(simuMN_noise_200(ii,:), 1000, 200);
    v2(ii,:) = resample(v2_200(ii,:), 1000, 200);
end
if size(simuMN_noise,2) > ecgLength
    simuMN_noise = simuMN_noise(:,1:ecgLength);
end

% Transform to the 15 leads
%1)Obtain augmented unipolar limb leads
simuMN_noise_8 = leadcalc(simuMN_noise,'stan');% V1,V2,V3,V4,V5,V6,I,II
simuMN_noise_12 = leadcalc(simuMN_noise_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

simuMN_noise_15=vertcat(simuMN_noise_8(7,:),simuMN_noise_8(8,:),simuMN_noise_12(12,:),...
    -simuMN_noise_12(9,:),-simuMN_noise_12(7,:),simuMN_noise_12(11,:),simuMN_noise_8(1:6,:),simuMN_noise);

end