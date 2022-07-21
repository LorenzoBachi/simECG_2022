function [rr,hrmean,ecgParameters] = simECG_generate_sinus_rhythm(N,ecgParameters)
%
% rr = simPAF_gen_SR_RR_intervals() returns synthetic SR RR intervals. 
% Ventricular rhythm during SR is simulated using RR interval generator 
% proposed by McSharry et al. (2003) in which both the impact of parasympathetic 
% stimulation (respiratory sinus arrhythmia) and baroreflex regulation
% (Mayer waves) is modeled by a bimodal power spectrum.

%Last update: CPerez 05/2022

if N < 2  % something does not work when N = 1
    N = 2;
end

switch ecgParameters.ESTflag
    case 0 %Any case
        hrmean = randi([50,90]);         % Generate heart rate from an interval of [50-80] bpm
        hrstd = randi([5,30])/10;               % Generate SD of heart rate from an interval of [0.5-3] bpm
        rrstd = hrstd/60;
        lfhfratio = randi([5,20])/10;           % Generate LF/HF ratio [0.5 - 2]
        respRate = randi([20,50])/100;          % Respiratory rate [0.2 - 0.5]
        MayerFreq = 0.1;                        % Mayer waves
        % Define frequency parameters for RR process
        floStd = 0.1;
        fhiStd = 0.1;
        
        % Respiration pattern
        ecgParameters.Fr = ones(1,ceil(N*respRate)).*respRate; %number of cycles        
        
    case 1 %EST
        %RR information
        rrstd = randi([5,10])*1e-3; %seconds
        MayerFreq = 0.1; % Mayer waves 
        floStd = 0.1;
        fhiStd = 0.1;
        hrmean = [];
        lfhfratio = [];
        
        
        %2) Respiration pattern
        Fr_b = []; Fr_e = []; Fr_r = []; Fr_b2 = [];
        %--> Basal
        T = 0; bins = 1;
        while T < ecgParameters.Basal
            Fr_b = repmat(ecgParameters.Frini,1, bins);
            T = sum(1./Fr_b);
            bins = bins +1;
        end
        
        %--> Exercise
        T = 0; bins = 2;
        while T < ecgParameters.Exercise
            Fr_e = linspace(ecgParameters.Frini, ecgParameters.Frpeak, bins);
            T = sum(1./Fr_e);
            bins = bins +1;
        end

        %--> Recovery
        T = 0; bins = 2;
        while T < ecgParameters.Recovery
            Fr_r = linspace(ecgParameters.Frpeak - 0.001, ecgParameters.Frend, bins);
            T = sum(1./Fr_r);
            bins = bins +1;
        end
        
        %--> Basal2
        T = 0; bins = 1;
        while T < ecgParameters.Basal2
            Fr_b2 =  repmat(ecgParameters.Frend,1, bins);
            T = sum(1./Fr_b2);
            bins = bins +1;
        end
        
        ecgParameters.Fr = [Fr_b Fr_e Fr_r(2:end) Fr_b2];                     
end

% Compute rr process
[rr, ecgParameters] = simECG_RR_process(MayerFreq,floStd,fhiStd,lfhfratio,hrmean,rrstd,ecgParameters);
hrmean = 60./ecgParameters.RRmean;
end