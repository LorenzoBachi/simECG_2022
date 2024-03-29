function [rr, simECGdata] = simECG_RR_process(flo, flostd, fhistd, lfhfratio, hrmean, rrstd,simECGdata)
% [] = simECG_RR_process() returns sinus rhythm RR intervals. RR intervals
% are generated according to the principle reported in the paper by
% McSharry et al. (2003) "A dynamical model for generating synthetic
% electrocardiogram signals" and expanded upon in "ECG Modeling for
% Simulation of Arrhythmias in Time-Varying Conditions" (2023). The same
% respiratory frequency using as calculating the RR interval as applying
% changes in the ECG morphology. 
% 
% Input arguments:
% flo - Mayer waves center frequency.
% flostd - standard deviation of Mayer waves frequency.
% fhistd - standard deviation of parasympathetic stimulation frequency.
% lfhfratio - low frequency to high frequency ratio.
% hrmean - mean heart rate.
% rrstd - desired standard deviation of RR intervals.
% simECGdata - struct of ECG simulation parameters defined in the main
% script.
% 
% Output arguments:
% rr - simulated sinus rhythm RR series.
% simECGdata - struct of ECG simulation parameters defined in the main
% script (updated).
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

%Last update: CPerez 05/2022

%Initial values
Fr = simECGdata.Fr;
fs = 2;

%1) Even-Sampled input filters and select corresponded value of Fr
FrTransient = ones(1,12).*0.2; %Transient period of 60 seconds
LTransient = sum(1./FrTransient);
Frall = [FrTransient Fr Fr(end)];
TFr = cumsum([0 1./Frall(1:end-1)]);
[tRR] = resamp(Frall', TFr', 1, fs);
Frn = []; %even-sampled respiratory signal
for ii = 1:length(tRR)
    nFr = find( TFr <= tRR(ii),1,'last');
    Frn(ii) = Frall(nFr);
end

v = randn(1,length(tRR)); %even-sampled white noise with std = 1

%Time-frequency parameters
w1 = 2*pi*flo/fs; %constant
w2 = 2*pi.*Frn./fs; %vector
c1 = flostd/fs; 
c2 = fhistd./fs;
P1 = 1;
if simECGdata.ESTflag
    %1) LFHF ratio pattern
        bias = 0.05;
        LF = rand(1,5)*2*bias +([0.74 0.82 0.42 0.20 0.30] - bias); %LF for basal, 60, 80, 100, recov for MaxT from Hernando2018
        HF = rand(1,5)*2*bias +([0.36 0.10 0.18 0.44 0.48] - bias);   %HF for basal, 60, 80, 100, recov for MaxT from Hernando2018
        RR60 = ceil(simECGdata.Exercise*0.6);
        RR80 = ceil(simECGdata.Exercise*0.8);
        Plfexercise = [linspace(LF(1),LF(2),RR60), linspace(LF(2),LF(3),(RR80-RR60)),linspace(LF(3),LF(4),(simECGdata.Exercise-RR80))];
        Phfexercise = [linspace(HF(1),HF(2),RR60), linspace(HF(2),HF(3),(RR80-RR60)),linspace(HF(3),HF(4),(simECGdata.Exercise-RR80))];
        Plf = [repmat(LF(1),1,sum(1./FrTransient)), repmat(LF(1),1,simECGdata.Basal), Plfexercise(1:simECGdata.Exercise), linspace(LF(4),LF(5),simECGdata.Recovery),repmat(LF(5),1,simECGdata.Basal2)];
        Phf = [repmat(HF(1),1,sum(1./FrTransient)), repmat(HF(1),1,simECGdata.Basal), Phfexercise(1:simECGdata.Exercise), linspace(HF(4),HF(5),simECGdata.Recovery),repmat(HF(5),1,simECGdata.Basal2)];     
        lfhfratio = Plf./Phf; %LFHF ratio  
        lfhfratio = interp1(1:length(lfhfratio),lfhfratio',1:length(tRR),'linear','extrap'); %to have the same length as RR
end
P2 = 1./lfhfratio;


%2) Calculate RRV: filters
x = []; %RRV

%LF
%--> cuadratic-logarithmic filter
texp = (-6000:6000)./fs;
fexp_LF = 1./(2*c1^2.*texp.^2 + 1);
limit = find(fexp_LF >= 0.05,1); % 5%
texp = -abs(texp(limit)):1/fs:abs(texp(limit));
fexp_LF = 1./(2*c1^2.*texp.^2 + 1);
tfilter = (0:length(texp)-1)./fs;

%--> Filter aplication 
A_LF = nthroot(32/pi^5,4)*sqrt(c1);
hLF = A_LF.*fexp_LF.*cos(w1.*tfilter);
hLFout = filter(sqrt(P1).*hLF,1,v);

%HF
Q = (0.2/fs)/c2; %Q factor = ratio between center frequency and bandwidth
if simECGdata.ESTflag
    for ii = 1:length(tRR)
        c2new(ii) = (Frn(ii)/fs)/Q; 
        %--> cuadratic-logarithmic filter
        texp = (-6000:6000)./fs;
        fexp_HF = 1./(2*c2new(ii)^2.*texp.^2 + 1);
        limit = find(fexp_HF >= 0.05,1); % 5%
        texp = -abs(texp(limit)):1/fs:abs(texp(limit));
        fexp_HF = 1./(2*c2new(ii)^2.*texp.^2 + 1);
        tfilter = (0:length(texp)-1)./fs;
        
        %--> Filter aplication 
        A_HF = nthroot(32/pi^5,4)*sqrt(c2new(ii)); 
        hHF = A_HF.*fexp_HF.*cos(w2(ii).*tfilter);
        hHFout = filter(sqrt(P2(ii)).*hHF,1,v);
        
        %--> Both filter
        x(ii) = hLFout(ii) + hHFout(ii);
    end
    
else
    fexp_HF = 1./(2*c2^2.*texp.^2 + 1);
    A_HF = nthroot(32/pi^5,4)*sqrt(c2); 
    hHF = A_HF.*fexp_HF.*cos(w2(1).*tfilter);
    hHFout = filter(sqrt(P2).*hHF,1,v);
    x = hLFout + hHFout;
end

x = x(sum(1./FrTransient)*fs-1:end); %reject the transient response


%3) RR interval
xstd = std(x);
ratio = rrstd/xstd;

%Pattern RRmean 
if simECGdata.ESTflag
    RR_b = repmat(simECGdata.RRini,1,simECGdata.Basal*fs);
    RR_e = linspace(simECGdata.RRini,simECGdata.RRpeak, simECGdata.Exercise*fs);
    RR_r = linspace(simECGdata.RRpeak,simECGdata.RRend, simECGdata.Recovery*fs);
    RR_b2 = repmat(simECGdata.RRend,1,simECGdata.Basal2*fs);
    rrMean = [RR_b RR_e RR_r RR_b2];
    
else
    rrMean = 60/hrmean;
end
L = simECGdata.Duration*fs;
trrMean = (0:simECGdata.Duration*fs-1)./fs;

%sum of RRv to obtain the desired RR series --> even-sampled
simECGdata.RRmean = rrMean;
rrSeries = rrMean + x(1:L).*ratio;


%4) Compute the RR interval
rr = [];
p = 1;
pos = p;
while sum(rr) <= simECGdata.Duration-5 %bias of 5 seconds because Fr signal has to be always larger
    rr = [rr rrSeries(p)];
    p = find(trrMean>=sum(rr),1);
    pos = [pos p];
end


%5) Real parameters
simECGdata.RR = rr;
simECGdata.Duration = sum(rr);

if simECGdata.ESTflag
    simECGdata.RRmean = rrMean(pos);
    simECGdata.peak = sum(rr(1:find(cumsum(rr) <= simECGdata.peak,1,'last')));
end
        
end