function [simuMN_noise_15] = simECG_Muscular_Noise(ecgLength, ecgParameters, noiseRMS)
% simuMN_noise = simECG_Muscular_Noise() returns a simulated muscular noise
% signal.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 07/2022

% 1) Load dictionary AR(p) model and select the parameters that model the quasy-stationay part of
% the simulated MN signal
load('DATA_AR_MN_Dictionary.mat')
ARp = squeeze(AR_MN(randi([1,25]),:,:));
fs = ecgParameters.fs;

% 2) AP(1) model to calculate the variance of the MN signal

%--> 1) Select the value of the pole and noise distribution
p(1) = rand(1)*(0.999-0.999) + 0.999;
p(1) = 0.999
b = 1-p(1);

%--> 2) Apply ARX model
u0 = noiseRMS*1e3; %in uV
N200 = ceil(ecgLength/5);

v1 = [];
out_APX1_200 = zeros(3,N200);

%Define time-varying amplitude
if ecgParameters.ESTflag
    N1 = length(1:ecgParameters.peak*ecgParameters.fs);
    N2 = length(ecgParameters.peak*ecgParameters.fs:ecgLength);
    ut = [rescale(1.5.^((1:N1)./(100*fs)),0,100),...
        rescale(flip(1.5.^((1:N2)./(100*fs))),0,100)];%exponential pattern exercise stress test
else
    ut = zeros(3, N200);
end
u = (u0 + ut).*b;

v1 = randn(3, N200).*0.2; %Frank leads
out_APX1_200(:,1) = u0;

for ii=1:N200-1
    out_APX1_200(:,ii+1) = p*out_APX1_200(:,ii)+v1(:,ii)+u(:,ii);
end

out_APX1_200(out_APX1_200<0) = 0; %ReLU


% 3)AP(n)model to obtain the desired simulated muscular noise signal
%-->1) Resample to 200Hz
%-->2) Apply AR(n) filter 
v2_200 = [];
v2_200 = randn(size(out_APX1_200,1),size(out_APX1_200,2));
v2_200 = v2_200.*out_APX1_200;

simuMN_noise_200 = [];
for Li = 1:3
    simuMN_noise_200(Li,:) = filter(1,ARp(:,Li),v2_200(Li,:)')';
end

%-->3) Resample to 1000Hz
simuMN_noise = [];
v2 = [];
for ii = 1:3
    simuMN_noise(ii,:) = resample(simuMN_noise_200(ii,:), 1000, 200);
end
if size(simuMN_noise,2) > ecgLength
    simuMN_noise = simuMN_noise(:,1:ecgLength);
end

simuMN_noise = simuMN_noise.*1e-3; % in mV (accroding to the ECG signal)

% Transform to the 15 leads
%1)Obtain augmented unipolar limb leads
simuMN_noise_8 = leadcalc(simuMN_noise,'stan');% V1,V2,V3,V4,V5,V6,I,II
simuMN_noise_12 = leadcalc(simuMN_noise_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

simuMN_noise_15=vertcat(simuMN_noise_8(7,:),simuMN_noise_8(8,:),simuMN_noise_12(12,:),...
    -simuMN_noise_12(9,:),-simuMN_noise_12(7,:),simuMN_noise_12(11,:),simuMN_noise_8(1:6,:),simuMN_noise);

end