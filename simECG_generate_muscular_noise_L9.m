function [simuMN_15] = simECG_generate_muscular_noise_L9(ecgLength, ecgParameters, noiseRMS, signal)
% simuMN_noise = simECG_Muscular_Noise() returns a simulated muscular noise
% signal in mV.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 10/2022

% 1) Load dictionary AR(p) model and select the parameters that model the quasy-stationay part of
% the simulated MN signal
load('DATA_AR_MN_Dictionary_L9_v2.mat');

fs = ecgParameters.fs;
v1 = [];
N200 = ceil(ecgLength/5);

% %8 standard-leads (I,II,V1-V6)
% L = size(AR_MN,3);
%9 standard-leads (V1-V6,I,II,III)
L = size(AR_MN,3);



%--> 1) Select the value of the pole and time-varying model
nu = rand(1)*(0.9995-0.99) + 0.99;

%--> 2) Apply 1st model and then sum the different signals

if noiseRMS > 1 %info in dB
    DefnoiseRMS = noiseRMS;
    Rpos = round(cumsum(ecgParameters.RR).*fs); %in samples
    noiseRMS = zeros(1,size(signal,1));ppQRS = zeros(1,size(signal,1));
    for l = 1:size(signal,1)
        qrsAll = zeros(100,length(Rpos));
        for ii = 2:length(Rpos)-1
            qrsAll(:,ii) = signal(l,Rpos(ii)-50:Rpos(ii)+50-1);
        end
        Aqrs = mean(qrsAll,2);
        ppQRS(l) = max(Aqrs)-min(Aqrs);
        noiseRMS(l) = ((max(Aqrs)-min(Aqrs))/(10^(DefnoiseRMS/20)));%in mV
    end  
else %info in mV
    noiseRMS = repmat(noiseRMS,1,L);
end


u0 = noiseRMS*1e3; %in uV
simuMN_200 = [];

for Li = 1:L
    if ecgParameters.ESTflag %take a look!
        peak = fix(ecgParameters.peak*200);
        N1 = length(1:peak);
        N2 = length(peak:N200);
%         ut = [rescale(2.^((1:N1)./(100*fs)),(u0(Li))/3,3*(u0(Li))/3),...
%             rescale(2.^(-(1:N2)./(100*fs)),(u0(Li))/3,3*(u0(Li))/3)];%exponential pattern exercise stress test
        ut = [linspace(u0(Li)/3,u0(Li),N1) linspace(u0(Li),u0(Li)/3,N2)];%linear pattern exercise stress test
        ut = ut(:,1:N200);
        sigmax = u0(Li)*0.1; %varianze
    else
        ut = u0(Li) + zeros(1, N200);
        sigmax = u0(Li)*0.1; %varianze
    end
    
    sigmav = sigmax.*(1-nu^2); %remember that the var_out = var_in/(1 - nu^2)
    v1 = randn(1, N200).*sqrt(sigmav); %9standard-leads
    
    out_1st_200 = zeros(Li,N200);
    
    for ii=1:N200-1
        out_1st_200(:,ii+1) = nu*out_1st_200(:,ii)+v1(:,ii);
    end
    
    sigmaw = max(0,out_1st_200 + ut); %all sum and ReLU
    
    
    % --> AR filter: with a random walk model
    v2_200 = [];
    v2_200 = randn(size(sigmaw,1),size(sigmaw,2));
    v2_200 = v2_200.*sigmaw;
        
    aAR = squeeze(AR_MN(randi([1,25],1),:,:));
    nSteps = ceil(N200/(10*200)); %each 10 seconds
    pnew = zeros(nSteps,4,3);
    
    
    pIni = 1;
    pEnd = 10*200;
    for ii = 1:nSteps
        [z,p,k] = tf2zpk(1,aAR(:,Li)); %poles of random-selected AR coefficients
        pnew(:,1,Li) = RandomWalk(nSteps,p(1),0.005, 0.01);
        pnew(:,3,Li) = RandomWalk(nSteps,p(3),0.005, 0.01);
        pnew(:,2,Li) = conj(pnew(:,1,Li));
        pnew(:,4,Li) = conj(pnew(:,3,Li));
        [b,aARnew] = zp2tf(z,pnew(ii,:,Li),k);
        simuMN_200(Li,pIni:pEnd) = filter(1,aARnew,v2_200(Li,pIni:pEnd)')';
        
        %Apply random walk
        if ii == nSteps-1
            pIni = pEnd+1;
            pEnd = length(v2_200);
        else
            pIni = pEnd+1;
            pEnd = pIni + 10*200 -1;
        end
    end
end

% --> Resample to 1000Hz
simuMN = [];
v2 = [];
for ii = 1:L
    simuMN(ii,:) = resample(simuMN_200(ii,:), 1000, 200);
end

simuMN = simuMN(:,1:ecgLength).*1e-3; %in mV ->9-standard leads

% Transform to the 15 leads
simuMN_12 = leadcalc(simuMN,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
simuMN_xyz = leadcalc(simuMN,'hcuzst');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
simuMN_15=vertcat(simuMN(7:9,:),-simuMN_12(9,:),simuMN_12(7,:),simuMN_12(11,:),...
    simuMN(1:6,:),simuMN_xyz);
end