function [rr,hrmean,simECGdata] = simECG_generate_sinus_rhythm(N,simECGdata)
% rr = simPAF_gen_SR_RR_intervals() returns synthetic SR RR intervals.
% Ventricular rhythm during SR is simulated using a novel RR interval
% generator based on McSharry et al. (2003) in which both the impact of
% parasympathetic stimulation (respiratory sinus arrhythmia) and baroreflex
% regulation (Mayer waves) is modeled by a bimodal power spectrum. In this
% model, time-verying respiration affects parasympathetic modulation. For
% additional details, please see "ECG Modeling for Simulation of
% Arrhythmias in Time-Varying Conditions" (2023).
% 
% Input arguments:
% N - number of normal, sinus RR interval to generate.
% simECGdata - struct of ECG simulation parameters defined in the main
% script.
% 
% Output arguments:
% rr - simulated sinus rhythm RR series
% hrmean - mean heart rate
% simECGdata - struct of ECG simulation parameters defined in the main
% script (updated).
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

%Last update: CPerez 05/2022

if N < 2  % something does not work when N = 1
    N = 2;
end

switch simECGdata.ESTflag
    case 0 %Any case
        hrmean = randi([50,80]);        % Generate heart rate from an interval of [50-80] bpm
        hrstd = randi([5,30])/10;       % Generate SD of heart rate from an interval of [0.5-3] bpm
        rrstd = hrstd/60;
        lfhfratio = randi([5,20])/10;   % Generate LF/HF ratio [0.5 - 2]
        respRate = randi([20,50])/100;	% Respiratory rate [0.2 - 0.5]
        MayerFreq = 0.1;                % Mayer waves
        % Define frequency parameters for RR process
        floStd = 0.1;
        fhiStd = 0.1;
        
        % Respiration pattern
        simECGdata.Fr = ones(1,ceil(N*respRate)).*respRate; %number of cycles
        simECGdata.Duration = N;
        
        
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
        while T < simECGdata.Basal
            Fr_b = repmat(simECGdata.Frini,1, bins);
            T = sum(1./Fr_b);
            bins = bins +1;
        end
        
        %--> Exercise
        T = 0; bins = 2;
        while T < simECGdata.Exercise
            Fr_e = linspace(simECGdata.Frini, simECGdata.Frpeak, bins);
            T = sum(1./Fr_e);
            bins = bins +1;
        end

        %--> Recovery
        T = 0; bins = 2;
        while T < simECGdata.Recovery
            Fr_r = linspace(simECGdata.Frpeak - 0.001, simECGdata.Frend, bins);
            T = sum(1./Fr_r);
            bins = bins +1;
        end
        
        %--> Basal2
        T = 0; bins = 1;
        while T < simECGdata.Basal2
            Fr_b2 =  repmat(simECGdata.Frend,1, bins);
            T = sum(1./Fr_b2);
            bins = bins +1;
        end
        
        simECGdata.Fr = [Fr_b Fr_e Fr_r(2:end) Fr_b2];                     
end

% Compute rr process
[rr, simECGdata] = simECG_RR_process(MayerFreq,floStd,fhiStd,lfhfratio,hrmean,rrstd,simECGdata);
hrmean = 60./simECGdata.RRmean;
end