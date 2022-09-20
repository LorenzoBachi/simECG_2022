function [simuMNMA_15, pnew] = simECG_generate_muscularN_motionA(ecgLength, ecgParameters, noiseRMS)
% simuMNMA_15 returns a simulated muscular noise with motion artifacts signal in mV.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 09/2022

% 1) Load dictionary AR(p) model and select the parameters that model the quasy-stationay part of
% the simulated MN signal
load('DATA_AR_MN_Dictionary.mat');
fs = ecgParameters.fs;
v1 = [];
N200 = ceil(ecgLength/5);

%--> 1) Select the value of the pole and time-varying model
nu = rand(1)*(0.9999-0.99) + 0.99;

%--> 2) Apply 1st model and then sum the different signals
u0 = noiseRMS*1e3; %in uV

if ecgParameters.ESTflag
    peak = (ecgParameters.peak*fs)/5;
    N1 = length(1:peak);
    N2 = length(peak:N200);
    ut = [rescale(2.^((1:N1)./(100*fs)),u0/3,3*u0/3),...
        rescale(2.^(-(1:N2)./(100*fs)),u0/3,3*u0/3)];%exponential pattern exercise stress test
    ut = repmat(ut,3,1);
    stdv = [linspace((u0/3)*0.3, (3*u0/3)*0.3, N1),linspace((3*u0/3)*0.3, (u0/3)*0.3, N2)];
else
    ut = u0 + zeros(3, N200);
    stdv = u0*0.3;
end

sigma_v1 = stdv.*sqrt(1-nu^2); %remember that the var_out = var_in/(1 - nu^2)
v1 = randn(3, N200).*sigma_v1; %Gaussian noise - Frank leads

out_1st_200 = zeros(3,N200);

for ii=1:N200-1
    out_1st_200(:,ii+1) = nu*out_1st_200(:,ii)+v1(:,ii);
end

stdw = max(1,out_1st_200 + ut); %all sum and ReLU


% --> Inputs AR filter: Gaussian noise + Bernoulli
v2_200 = [];
v2_200 = randn(size(stdw,1),size(stdw,2)); %Gaussian noise - Frank leads
v2_200 = v2_200.*stdw;
bernogauss = simECG_Bernoulli_Gaussian(N200, ecgParameters.MA_Prob, stdv);
bernogauss_conv = simECG_Bernoulli_Gaussian_convolution(N200, bernogauss);

v2B_200 = v2_200 + bernogauss_conv;

simuMNMA_200 = [];

% --> AR filter: with a random walk model
aAR = squeeze(AR_MN(randi([1,25],1),:,:));
nSteps = ceil(N200/(10*200)); %each 10 seconds
pnew = zeros(nSteps,4,3);

for Li = 1:3
    [z,p,k] = tf2zpk(1,aAR(:,Li)); %poles of random-selected AR coefficients
    pnew(:,1,Li) = RandomWalk(nSteps,p(1),0.02, 0.01);
    pnew(:,3,Li) = RandomWalk(nSteps,p(3),0.02, 0.01);
    pnew(:,2,Li) = conj(pnew(:,1,Li));
    pnew(:,4,Li) = conj(pnew(:,3,Li));

    pIni = 1;
    pEnd = 10*200;
    for ii = 1:nSteps
        [b,aARnew] = zp2tf(z,pnew(ii,:,Li),k);
        simuMNMA_200(Li,pIni:pEnd) = filter(1,aARnew,v2B_200(Li,pIni:pEnd)')';
        if ii == nSteps-1
            pIni = pEnd+1;
            pEnd = length(v2B_200);
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
simuMNMA = [];
v2 = [];
for ii = 1:3
    simuMNMA(ii,:) = resample(simuMNMA_200(ii,:), 1000, 200);
end
if size(simuMNMA,2) > ecgLength
     simuMNMA = simuMNMA(:,1:ecgLength);
end

simuMNMA = simuMNMA.*1e-3; %in mV

% Transform to the 15 leads
simuMNMA_8 = leadcalc(simuMNMA,'stan');% V1,V2,V3,V4,V5,V6,I,II
simuMNMA_12 = leadcalc(simuMNMA_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

simuMNMA_15 = vertcat(simuMNMA_8(7,:),simuMNMA_8(8,:),simuMNMA_12(12,:),...
    -simuMNMA_12(9,:),-simuMNMA_12(7,:),simuMNMA_12(11,:),simuMNMA_8(1:6,:),simuMNMA);


end