function [rr ESTpeak fa ecgnr] = simECG_get_real_RR_intervals(AFrhythm, rrLength)
% [] = simPAF_construct_real_RR() returns RR interval series, constructed
% using real RR intevals obtained from the MIT-BIH Normal Sinus Rhythm database
% (SR) and the Long Term Atrial Fibrillation database (AF). 

%In case of EST signals, returns RR interval series, constructed
% using real RR series obtained from the FINCAVAS Exercise Stress Test
% database (Normal Sinus Rhythm). %Cris 03/2022

switch AFrhythm % If AFrhythm == 1 - AF, SR otherwise.
    case 0
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
    case 1
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
        
    case 2 %CPerez 02/2022
        clear rrLength
        sigNum = randi([1 42])  % Select randomly a single SR RR signal from 42 possible
%         sigNum = 65
        load('DATA_RR_SR_EST_real_ECG-LR');
        rr = DATAestSR(sigNum).cRR';  %in seconds
        fa = DATAestSR(sigNum).fs;
        ecgnr = DATAestSR(sigNum).ecgnr;
        ESTpeak = DATAestSR(sigNum).peak; %in seconds
        clear DATAestSR
end

end
        
