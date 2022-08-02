function [simuMN_15] = simECG_Muscular_Noise(ecgLength, ecgParameters, noiseRMS)
% simuMN_noise = simECG_Muscular_Noise() returns a simulated muscular noise
% signal.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 07/2022

% 1) Load dictionary AR(p) model and select the parameters that model the quasy-stationay part of
% the simulated MN signal
load('DATA_AR_MN_Dictionary.mat')
ARp = squeeze(AR_MN(randi([1,25]),:,:));
fs = ecgParameters.fs;
v1 = [];
N200 = ceil(ecgLength/5);

%--> 1) Select the value of the pole and time-varying model
nu = rand(1)*(0.9995-0.990) + 0.990
nu = 0.9999

%--> 2) Apply 1st model and then sum the different signals
u0 = noiseRMS*1e3; %in uV

if ecgParameters.ESTflag
    peak = (ecgParameters.peak*fs)/5;
    N1 = length(1:peak);
    N2 = length(peak:N200);
    ut = [rescale(2.^((1:N1)./(100*fs)),-u0/2,u0),...
        rescale(flip(2.^((1:N2)./(100*fs))),-u0/2,u0)];%exponential pattern exercise stress test
    ut = repmat(ut,3,1);
    stdw = (u0/4)*0.5;
    
else
    ut = zeros(3, N200);
    stdw = u0*0.5;
end

sigma_v1 = stdw*sqrt(1-nu^2); %remember that the var_out = var_in/(1 - nu^2)
v1 = randn(3, N200).*sigma_v1; %Frank leads

out_1st_200 = zeros(3,N200);

for ii=1:N200-1
    out_1st_200(:,ii+1) = nu*out_1st_200(:,ii)+v1(:,ii);
end

Allout_1st_200 = max(1,out_1st_200 + u0 + ut); %all sum and ReLU


% --> AR filter
v2_200 = [];
v2_200 = randn(size(Allout_1st_200,1),size(Allout_1st_200,2));
% v2_200 = normalize(v2_200','range',[-1 1])';
v2_200 = v2_200.*Allout_1st_200;

simuMN_200 = [];
for Li = 1:3
    simuMN_200(Li,:) = filter(1,ARp(:,Li),v2_200(Li,:)')';
end

% --> Resample to 1000Hz
simuMN = [];
v2 = [];
for ii = 1:3
    simuMN(ii,:) = resample(simuMN_200(ii,:), 1000, 200);
end
if size(simuMN,2) > ecgLength
     simuMN = simuMN(:,1:ecgLength);
end

simuMN = simuMN.*1e-3; %inmV

% Transform to the 15 leads
%1)Obtain augmented unipolar limb leads
simuMN_8 = leadcalc(simuMN,'stan');% V1,V2,V3,V4,V5,V6,I,II
simuMN_12 = leadcalc(simuMN_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

simuMN_15=vertcat(simuMN_8(7,:),simuMN_8(8,:),simuMN_12(12,:),...
    -simuMN_12(9,:),-simuMN_12(7,:),simuMN_12(11,:),simuMN_8(1:6,:),simuMN);

end