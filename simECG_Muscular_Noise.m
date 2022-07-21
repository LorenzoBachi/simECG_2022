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
p = rand(1)*(0.9999-0.999) + 0.999;
b = 1;
a = [1 -p]; %according to me

%--> 2) Apply the filter
v1 = randn(3, ecgLength + 50000); %Frank leads
out_AP1 = filter(b,a,v1')';
v1(:,1:50000) = [];
out_AP1(:,1:50000) = [];

out_AP1 = (out_AP1 - 0)./sqrt(var(v1,[],2)/(1-p^2));%Normalize

if ecgParameters.ESTflag
    for Li = 1:3 %frank leads
        tNew = linspace(0,ecgParameters.peak,size(patternMN.signal(Li,1:patternMN.peak),2));
        constantVar_e(Li,:)  = interp1(tNew,patternMN.signal(Li,1:patternMN.peak),(0:round(ecgParameters.peak*fs)-1)./fs);
        
        tNew = linspace(ecgParameters.peak,ecgLength/fs,size(patternMN.signal(Li,patternMN.peak+1:end),2));
        constantVar_r(Li,:)  = interp1(tNew,patternMN.signal(Li,patternMN.peak+1:end),(round(ecgParameters.peak*fs):ecgLength-1)./fs);
    end
    
    constantVar = [constantVar_e, constantVar_r];
else
    constantVar = min(patternMN.signal,[],2); %take a look! %Cris 07/2022
end


out_AP1 = (constantVar + out_AP1.^2); %variance
out_AP1 = fillmissing(out_AP1,'movmean',5000,2); %due to interp1 to calculate constantVar_e and constantVar_r

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

% Transform to the 15 leads
%1)Obtain augmented unipolar limb leads
simuMN_noise_8 = leadcalc(simuMN_noise,'stan');% V1,V2,V3,V4,V5,V6,I,II
simuMN_noise_12 = leadcalc(simuMN_noise_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

simuMN_noise_15=vertcat(simuMN_noise_8(7,:),simuMN_noise_8(8,:),simuMN_noise_12(12,:),...
    -simuMN_noise_12(9,:),-simuMN_noise_12(7,:),simuMN_noise_12(11,:),simuMN_noise_8(1:6,:),simuMN_noise);

end