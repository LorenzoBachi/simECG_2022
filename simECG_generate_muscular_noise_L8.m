function [simuMN_15] = simECG_generate_muscular_noise(ecgLength, ecgParameters, noiseRMS)
% simuMN_noise = simECG_Muscular_Noise() returns a simulated muscular noise
% signal in mV.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 10/2022

% 1) Load dictionary AR(p) model and select the parameters that model the quasy-stationay part of
% the simulated MN signal
load('DATA_AR_MN_Dictionary_L8.mat');
fs = ecgParameters.fs;
v1 = [];
N200 = ceil(ecgLength/5);

if ecgParameters.MA_Flag %Thumb-ECG
    L = 1;
else %8 standard-leads (I,II,V1-V6)
    L = size(AR_MN,3);
end


%--> 1) Select the value of the pole and time-varying model
nu = rand(1)*(0.9995-0.99) + 0.99;

%--> 2) Apply 1st model and then sum the different signals
u0 = noiseRMS*1e3; %in uV

if ecgParameters.ESTflag
    peak = (ecgParameters.peak*fs)/5;
    N1 = length(1:peak);
    N2 = length(peak:N200);
    ut = [rescale(2.^((1:N1)./(100*fs)),u0/3,3*u0/3),...
        rescale(2.^(-(1:N2)./(100*fs)),u0/3,3*u0/3)];%exponential pattern exercise stress test
    ut = repmat(ut,L,1);
    ut = ut(:,1:N200);
    sigmax = [linspace((u0/3)*0.3, (3*u0/3)*0.3, N1),linspace((3*u0/3)*0.3, (u0/3)*0.3, N2)];
    sigmax = sigmax(:,1:N200);    
else
    ut = u0 + zeros(L, N200);
    sigmax = u0*0.3;
end

sigmav = sigmax.*sqrt(1-nu^2); %remember that the var_out = var_in/(1 - nu^2)
v1 = randn(L, N200).*sigmav; %8standard-leads

out_1st_200 = zeros(L,N200);

for ii=1:N200-1
    out_1st_200(:,ii+1) = nu*out_1st_200(:,ii)+v1(:,ii);
end

sigmaw = max(1,out_1st_200 + ut); %all sum and ReLU


% --> AR filter: with a random walk model
v2_200 = [];
v2_200 = randn(size(sigmaw,1),size(sigmaw,2));
% v2_200 = normalize(v2_200','range',[-1 1])';
v2_200 = v2_200.*sigmaw;

simuMN_200 = [];


aAR = squeeze(AR_MN(randi([1,25],1),:,:));
nSteps = ceil(N200/(10*200)); %each 10 seconds
pnew = zeros(nSteps,4,3);

for Li = 1:L
    [z,p,k] = tf2zpk(1,aAR(:,Li)); %poles of random-selected AR coefficients
    pnew(:,1,Li) = RandomWalk(nSteps,p(1),0.02, 0.01);
    pnew(:,3,Li) = RandomWalk(nSteps,p(3),0.02, 0.01);
    pnew(:,2,Li) = conj(pnew(:,1,Li));
    pnew(:,4,Li) = conj(pnew(:,3,Li));

    pIni = 1;
    pEnd = 10*200;
    for ii = 1:nSteps
        [b,aARnew] = zp2tf(z,pnew(ii,:,Li),k);
%         simuMN_200(Li,pIni:pEnd) = filter(1,aARnew,v2_200(Li,pIni:pEnd)')';
        ar_model = idpoly(aARnew, 'integratenoise', ecgParameters.MA_Flag);
        simuMN_200(Li,pIni:pEnd) = sim(ar_model, v2_200(Li,pIni:pEnd)')';
        
        if ii == nSteps-1
            pIni = pEnd+1;
            pEnd = length(v2_200);
        else
            pIni = pEnd+1;
            pEnd = pIni + 10*200 -1;
        end
    end
    
%     figure(3), subplot(2,2,Li), %to plot the different poles of each lead
%     plot([-1 1], [0 0],':k',[0 0], [-1 1],':k')
%     hold on, plot(p,'or'), plot(squeeze(pnew(:,:,Li)),'xb')
%     xlim([-1 1]), ylim([-1 1])
%     title({'Lead ' num2str(Li)})
end

% --> Resample to 1000Hz
simuMN = [];
v2 = [];
for ii = 1:L
    simuMN(ii,:) = resample(simuMN_200(ii,:), 1000, 200);
end

simuMN = simuMN(:,1:ecgLength).*1e-3; %in mV ->8-standard leads


simuMN_12 = leadcalc(simuMN,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
simuMN_xyz = leadcalc(simuMN,'hcuzst');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
simuMN_15=vertcat(simuMN(7:8,:),simuMN_12(12,:),-simuMN_12(9,:),simuMN_12(7,:),simuMN_12(11,:),...
    simuMN(1:6,:),simuMN_xyz);
end