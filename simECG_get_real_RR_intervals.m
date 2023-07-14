function [rr, ESTpeak, fa, ecgnr] = simECG_get_real_RR_intervals(rhythm, rrLength)
% [] = simECG_get_real_RR_intervals() returns a RR interval series
% constructed using real RR intevals obtained from the MIT-BIH Normal Sinus
% Rhythm database (SR) and the Long Term Atrial Fibrillation database (AF).
% In case of ECG stress test signals, returns RR interval series,
% constructed using real RR series obtained from the FINCAVAS Exercise
% Stress Test database (Normal Sinus Rhythm). %Cris 03/2022
%
% Input arguments:
% rhythm - type of rhythm: 0 for sinus rhythm, 1 for atrial fibrillation,
% 2 for stress test ECG.
% rrLength - number of desired real RR intervals
%
% Output arguments:
% rr - real RR series
% ESTpeak, fa, ecgnr - parameters of real stress test sampled from the
% FINCAVAS database
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

ESTpeak = 0; fa = 0; ecgnr = 0;
switch rhythm 
    case 0 % sinus rhythm
        sigNum = randi([1 18]);  % Select randomly a single SR RR signal from 18 possible
        load('DATA_RR_SR_real');
        rrSR = DATArrSR(sigNum).rrSR;
        clear DATArrSR
            
        if rrLength > length(rrSR)   % Concatenate if RR series is shorter than the desired length
            nCycles = ceil(rrLength/length(rrSR));
            rr = [];
            for i = 1:nCycles
                rr = [rr rrSR];
            end
            rr = rr(1:rrLength);  
        else
            rrStart = randi([1 (length(rrSR)-rrLength)]); % Randomly select start position
            rr = rrSR(rrStart+1:rrStart+rrLength);
        end   
    case 1 % atrial fibrillation
        sigNum = randi([1 69]);  % Select randomly a single AF RR signal from 69 possible
        load('DATA_RR_AF_real');
        rrAF = DATArrAF(sigNum).rrAF;
        clear DATArrAF
            
        if rrLength > length(rrAF)   % Concatenate if RR series is shorter than the desired length
            nCycles = ceil(rrLength/length(rrAF));
            rr = [];
            for i = 1:nCycles
                rr = [rr rrAF];
            end
            rr = rr(1:rrLength);  
        else
            rrStart = randi([1 (length(rrAF)-rrLength)]); % Randomly select start position
            rr = rrAF(rrStart+1:rrStart+rrLength);
        end   
        
    case 2 %CPerez 02/2022, stress test ecg
        clear rrLength
        sigNum = randi([1 42]);  % Select randomly a single SR RR signal from 42 possible
%         sigNum = 65
        load('DATA_RR_SR_EST_real_ECG-LR');
        rr = DATAestSR(sigNum).cRR';  %in seconds
        fa = DATAestSR(sigNum).fs;
        ecgnr = DATAestSR(sigNum).ecgnr;
        ESTpeak = DATAestSR(sigNum).peak; %in seconds
        clear DATAestSR
end

end
        
