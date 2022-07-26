function [simuMN_noise_15] = simECG_Muscular_Noise(ecgLength, ecgParameters)
% simuMN_noise = simECG_Muscular_Noise() returns a simulated muscular noise
% signal.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 07/2022

% 1) Load dictionary AR(p) model and select the parameters that model the quasy-stationay part of
% the simulated MN signal
load('infoMuscularNoise_EST.mat')
ARp = squeeze(AR_MN_Dictionary(randi([1,25]),:,:));
fs = ecgParameters.fs;

% 2) AP(1) model to calculate the variance of the MN signal

v1 = [];
out_AP1 = [];
%--> 1) Select the value of the pole and noise distribution
p(1) = rand(1)*(0.9999-0.999) + 0.999
% p(1) = 0.9999;
b = 1-p(1);
a = [1, -p(1)]; %according to me

%--> 2) Apply ARX model
u0 = 99;
u = (u0 + ones(3, ecgLength)).*b; %time-varying amplitude
v1 = randn(3, ecgLength); %Frank leads
v1in = v1 + u; 
out_AP1 = filter(1,a,v1in')';
out_AP1(out_AP1<0) = 0; %ReLU


% 3)AP(n)model to obtain the desired simulated muscular noise signal
%-->1) Resample to 200Hz
out_AP1_200 = [];
for ii = 1:3
    out_AP1_200(ii,:) = resample(out_AP1(ii,:), 200, 1000);
end

%-->2) Apply AR(n) filter
v2_200 = [];
v2_200 = randn(size(out_AP1_200,1),size(out_AP1_200,2));
v2_200 = v2_200.*out_AP1_200;

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

simuMN_noise = rescale(simuMN_noise,-10,10); %according to the simulator

% Transform to the 15 leads
%1)Obtain augmented unipolar limb leads
simuMN_noise_8 = leadcalc(simuMN_noise,'stan');% V1,V2,V3,V4,V5,V6,I,II
simuMN_noise_12 = leadcalc(simuMN_noise_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

simuMN_noise_15=vertcat(simuMN_noise_8(7,:),simuMN_noise_8(8,:),simuMN_noise_12(12,:),...
    -simuMN_noise_12(9,:),-simuMN_noise_12(7,:),simuMN_noise_12(11,:),simuMN_noise_8(1:6,:),simuMN_noise);

end