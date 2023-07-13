function [QRSindex, TendIndex, rr, multileadVA, ecgLength, simECGdata] = simECG_generate_multilead_VA(targets_beats, rr, realVAon, simECGdata,state_history)
% [] = simECG_gen_multilead_VA() returns 15-lead ventricular activity. A
% set of 100 15-lead ECGs with SR selected from the PTB Diagnostic ECG
% Database is used as a basis for modeling ventricular activity of
% supraventricular beats. The ECGs of the PTB database are first subjected
% to baseline removal and QRST delineation. The original T-waves are then
% resampled to a fixed width and, depending on the type of rhythm,
% width-adjusted to match prevailing heart rate. Since the original ECGs
% last just for about 2 min, QRST complexes are subjected to repeated
% concatenation until desired length of ECG is obtained. The TQ interval is
% interpolated using a cubic spline interpolation. The ventricular beats
% are modeled from ventricular beats extracted from the St Petersburg
% INCART 12-lead database, using Hermite functions as a basis.
%
% Input arguments:
% targets_beats - array of beat codes.
% rr - RR series of the simulated record.
% realVAon - 0 for fully simulated venricular activity, 1 for real atrial
% activity from the PTB database 
% simECGdata - struct of ECG simulation parameters defined in the main
% script.
% state_history - array of Markov chain states.
%
% Output arguments:
% QRSindex - array of heartbeat time indexes (in samples).
% TendIndex: array of T end time indexes (in samples).
% rr - RR series of the simulated record (updated).
% multileadVA - simulated multilead ventricular activity.
% ecgLength - length of simulated ECG, in samples.
% simECGdata - struct of ECG simulation parameters defined in the main
% script (updated).
%
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

disp('Generating ventricular activity ...');

rrLength = length(rr);
rr = (round(rr*1000));    % rr intervals in miliseconds
rr = [rr(1) rr]';
Fs = simECGdata.fs;                      % sampling rate
rr_sec = rr/Fs;

%% QT interval adaptation
nfreq = 4; %Interpolation to 4 Hz;
timebeats = cumsum(rr(2:end))/1000; %In seconds
timebeats = [0;timebeats];
% resample RR series to an evenly spaced series
[pos,RR_ori]=resamp(rr_sec,timebeats*Fs,Fs,nfreq);%pos in samples, RR_ori in sec
pos = pos/Fs; %sec
RR_exp = NaN(size(RR_ori));
L = 300; %window size in seconds.
N = min(L*nfreq,length(RR_ori));
j2=[0:1:N-1]';
tau_s = 25;
tau = tau_s*nfreq; % tau in samples, according to sampling frequency
alpha = exp(-1/tau); %de los valores de la tesis de Esther
k=(1-alpha)/(1-alpha^N);

for i=N:size(RR_ori,1)
    RR_exp (i)= k*sum(exp(-j2/tau).*RR_ori(i-j2)); %la exponencial es decreciente, dando mayor peso al primero, que es el más actual
end
% kk = cumsum(k*exp(j2/tau)); %the sum of the weights w is equal to 1
% Regresson model: hyperbolic
beta_h = 0.49; %medium risk from C.Perez, CinC 2020
alpha_h= -0.09;
qt_hyper_s = beta_h + alpha_h./RR_exp; %sampled at 4Hz
qt_hyper_s(1:N-1) = qt_hyper_s(N); %Values for the first N beats are fixed to the first QT value
qt_hyper = interp1(pos,qt_hyper_s,timebeats,'spline'); %QT series in beats
qt_hyper = round(qt_hyper*Fs);

%% Signal amplitude scaling factors
ecg_amp = simECGdata.ecg_amp;
if ecg_amp < 0
    ecg_amp = 5;
end
%maximum and minimum ECG QRST amplitude scaling factor
scale_factor = ( rand * 0.5 ) + ecg_amp*0.15;
A = 0.5 + ecg_amp*0.15 ;
simECGdata.scale_factor = scale_factor;

%% Normal PQRST complexes generation
switch realVAon
    case 0 % Generate synthetic PQRST complexes
        %[data.pqrst, data.Qind, data.Rind, data.Sind, data.Tind] = simPAF_gen_syn_VA(100);
        [data.pqrst, data.Qind, data.Rind, data.Sind, data.Tind] = simECG_generate_XYZ_VA(100,[scale_factor,A]);%% Changed in v012020
    case 1 % Load random patch of PQRST complexes extracted from PTB database
        sigNum = randi([1 100]);
        load('DATA_PQRST_real')
        data.pqrst = DATApqrst(sigNum).pqrst;
        data.Qind = DATApqrst(sigNum).Qind;
        data.Rind = DATApqrst(sigNum).Rind;
        data.Sind = DATApqrst(sigNum).Sind;
        data.Tind = DATApqrst(sigNum).Tind;
end
% Original PQRST is subjected to repeated concatenation until required
% number of beats is obtained
arraySize = length(data.pqrst(1,:,1));
pqrstLength = length(data.pqrst(1,1,:));
numLeads = size(data.pqrst,1);
pqrst = zeros(numLeads,rrLength,pqrstLength); %% Changed in v012020

j = 1;
for i = 1:rrLength+1 % Additonal is required to meet exact number of RR intervals
    pqrst(:,i,:) = data.pqrst(:,j,:);
    j = j + 1;
    if j > arraySize
        j = 1;
    end
end

%% Ventricular PQRST complexes generation
% Kors regression transformation
K = transpose([ 0.380,  -0.070,   0.110;...
    -0.070,   0.930,  -0.230;...
    -0.130,   0.060,  -0.430;...
    0.050,  -0.020,  -0.060;...
    -0.010,  -0.050,  -0.140;...
    0.140,   0.060,  -0.200;...
    0.060,  -0.170,  -0.110;...
    0.540,   0.130,   0.310]);
%length of a ventricular beat buffer
Lv = 800;
%ventricular QRS complexes increased width
W = rand * 0;
%ventricular beat amplitude factor
vaf = 0.75 + ecg_amp*0.05 + ((0.1*rand) - 0.05);
%multiform VPBs?
multiform_vpbs = simECGdata.multiform_vpbs;
switch realVAon
    case 0 % Generate ventricular complexes from the Hermite-logistic
        %function model
        load('DATA_ventricular_coefficients');
    case 1 % Load random patch of ventricular PQRST complexes from
        %Alcaraz's Database
        load('DATA_ventricular_real');
end
% number of ventricular beats that are present in the record
Nbvt = sum(targets_beats==4);
if Nbvt
    if multiform_vpbs ~= 0
        % multiform VPBs
        vpb = zeros(Nbvt,8,Lv);
        vpb_record_beat = zeros(2,Nbvt);
        for b = 1:Nbvt
            vpb_record_beat(1,b) = randi(10);
            if vpb_record_beat(1,b)==10
                vpb_record_beat(2,b) = randi(87);
            else
                vpb_record_beat(2,b) = randi(100);
            end
        end
        bb = 0;
        for k = 1:10
            vpb_to_add = vpb_record_beat(:,vpb_record_beat(1,:)==k);
            Nvpb_c = sum(vpb_record_beat(1,:)==k);
            if Nvpb_c
                switch realVAon
                    case 0
                        fitCoeff = DATAventricularCoefficients(k).fitCoeff;
                        fitCoeff(:,:,6+1) = fitCoeff(:,:,6+1)*(1+W);
                        fitCoeff(:,:,10+6) = round(fitCoeff(:,:,10+6)*(1+W));
                        vpb0 = simECG_generate_XYZ_ventricular_VA(fitCoeff,DATAventricularCoefficients(k).J);
                        vpb_R = DATAventricularCoefficients(k).R;
                    case 1
                        vpb0 = DATAventricularReal(k).qrst;
                        vpb_R = DATAventricularReal(k).R;
                end
                for b = 1:Nvpb_c
                    bb = bb + 1;
                    vpb(bb,:,:) = vpb0(vpb_to_add(2,b),:,:) .* vaf;
                end
            end
        end
        rp = randperm(Nbvt);
        vpb = vpb(rp,:,:);
    else
        % single VPB shape
        vpb = zeros(Nbvt,8,Lv);
        k = randi(10);
        if k ==10
            q = randi(87);
        else
            q = randi(100);
        end
        switch realVAon
            case 0
                fitCoeff = DATAventricularCoefficients(k).fitCoeff;
                fitCoeff(:,:,6+1) = fitCoeff(:,:,6+1)*(1+W);
                fitCoeff(:,:,10+6) = round(fitCoeff(:,:,10+6)*(1+W));
                vpb0 = simECG_generate_XYZ_ventricular_VA(fitCoeff,DATAventricularCoefficients(k).J);
                vpb_R = DATAventricularCoefficients(k).R;
            case 1
                vpb0 = DATAventricularReal(k).qrst;
                vpb_R = DATAventricularReal(k).R;
        end
        for k = 1:Nbvt
            vpb(k,:,:) = squeeze(vpb0(q,:,:)) .* vaf;
        end
    end
    vpb_15l = zeros(Nbvt,15,Lv);
    for k = 1:Nbvt
        vpb_15l(k,1:2,:) = vpb(k,1:2,:);
        vpb_15l(k,3,:) = vpb(k,2,:) - vpb(k,1,:);
        vpb_15l(k,4,:) = -vpb(k,1,:) + vpb(k,2,:)/2;
        vpb_15l(k,5,:) =  vpb(k,1,:) - vpb(k,2,:)/2;
        vpb_15l(k,6,:) =  vpb(k,2,:) - vpb(k,1,:)/2;
        vpb_15l(k,7:12,:) = vpb(k,3:8,:);
        vpb_15l(k,13:15,:) = K * squeeze(vpb(k,1:8,:));
    end
    v_idx = find(targets_beats == 4);
end

%% Activity generation

pqrstResampled=zeros(numLeads,1);
ecgSig = [];
%Load QRST complexes for all leads
for lead = 1:numLeads %% Changed in v012020
    %AFend = 0;
    % counters for supraventricular and ventricular beat numbers
    ks = 0;
    kv = 0;
    nVPB = 0;
    kv_history = ones(1,Nbvt);
    
    % extracting fiducial points location
    Qind = data.Qind; Rind0 = data.Rind - data.Qind;
    Sind = data.Sind - data.Qind; Tind = data.Tind - data.Qind;
    % preparing first beat
    if targets_beats(1)~=4
        % first beat is supraventricular
        QRSTtempNext = pqrst(lead,1,Qind+1:end);
        QRSTtempNext = simECG_correct_baseline(QRSTtempNext);
    else
        % first beat is ventricular
        vQRST = zeros(1,Lv);
        QRSTtempNext = vQRST;
    end
    for beatNr = 1:rrLength
        QRSTtemp = QRSTtempNext;
        if (beatNr < rrLength)
            if targets_beats(beatNr+1)~=4
                QRSTtempNext = pqrst(lead,ks+1,Qind+1:end);
                QRSTtempNext = simECG_correct_baseline(QRSTtempNext);
            else
                % ventricular beat
                vQRST = zeros(1,Lv);
                QRSTtempNext = vQRST;
            end
        end
        if targets_beats(beatNr)~=4
            % supraventricular beat (either normal or atrial
            % ectopic)
            ks = ks + 1;
            % R wave fiducial point is selected
            Rind = Rind0;
            % Divide PQRST into two parts
            QRST_PS = QRSTtemp(1,1:Sind); % The PQRS part
            QRST_ST = QRSTtemp(1,Sind+1:Tind); % The ST part
            % Resample T wave according to current RR interval
            % Signal is prolonged in order to avoid boundary effects of resampling
            QRST_ST  = [QRST_ST(1,1)*ones(1,10) QRST_ST QRST_ST(1,end)*ones(1,10)];
            % QRSTtemp  = [QRSTtemp(1,1)*ones(1,20) QRSTtemp QRSTtemp(1,end)*ones(1,20)]; %added by Alba
            
            % Alba: Hyperbolic model QT/RR regression
            ST_length = qt_hyper(beatNr);
            % Original version by Andrius:
            % Gradual prolongation of T wave after AF is terminated
            %if (targets_beats(beatNr) == 1)||(targets_beats(beatNr) == 3)
            %    if AFend == 1
            %        ST_length = round(length(QRST_ST)*sqrt((0.3+m*0.1)*rr(beatNr)/1000));
            %        m = m + 1;
            %        if m == 7
            %            AFend = 0;
            %        end
            %    else
            %        ST_length = round(length(QRST_ST)*sqrt(rr(beatNr)/1000));
            %    end
            %else
            %    ST_length = round(length(QRST_ST)*sqrt(0.35));
            %    % Alba: Bazzet's refers to the whole QT interval
            %    %  ST_length = round(QTc*sqrt(rr(beatNr)/1000));
            %    AFend = 1;
            %    m = 1;
            %end
            %QRST_ST = resample(QRST_ST, ST_length+round(20*sqrt(rr(beatNr)/1000))-length(QRST_PS), length(QRST_ST));  % T wave length correction QT = QTc*sqrt(RR) ~ T = Tc*sqrt(RR) %Bazzet's
            border = round(10*(ST_length-length(QRST_PS))/(length(QRST_ST)-20));
            QRST_ST = resample(QRST_ST, ST_length-length(QRST_PS), length(QRST_ST)-20);  % T wave length correction QT = QTc*sqrt(RR) ~ T = Tc*sqrt(RR)
            QRSTc = [QRST_PS QRST_ST(border:end-border)]; % Alba: adapted QRST
            %QRSTc = QRST_ST(round(20*sqrt(rr(beatNr)/1000))+1:end-round(20*sqrt(rr(beatNr)/1000))); %Alba: it is not the QTc but the original QT interva before correction
            
            % Find TQ length
            TQlength = rr(beatNr+1) - length(QRSTc);
            % Interpolate TQ interval only if the space between beats
            % is enough
            if  TQlength >= 2
                % Interpolate TQ interval
                x = [1 TQlength];
                xi = 1:1:TQlength;
                y = [QRSTc(1,end) QRSTtempNext(1, 1)];
                TQ = interp1(x,y,xi,'linear');
                QRSTcQ = [ QRSTc, TQ ];
            else
                QRSTcQ = QRSTc;
            end
            nn = length(QRSTcQ);
            
        else
            % Ventricular beat - activity placed after respiration
            % correction
            % in bigeminy / trigeminy, the same VPB is used
            if ~((state_history(max(beatNr-4,1))==4)&&(kv>0))
                if ~((state_history(beatNr-1)==state_history(beatNr))&&(state_history(beatNr)==5))
                    kv = kv + 1;
                end
            end
            % R wave fiducial point
            Rind = vpb_R;
            nVPB = nVPB + 1;
            kv_history(nVPB) = kv;
            vpb_temp = squeeze(vpb(kv_history(nVPB),:,:));
            vpb_module = vecnorm(vpb_temp);
            vpb_module = vpb_module / sum(vpb_module);
            vpb_end = length(vpb_temp);
            area = 0;
            for v = length(vpb_module):-1:round(length(vpb_module)*0.5)
                area = area + vpb_module(v);
                if area >= 0.05
                    vpb_end = v;
                end
            end
            QRSTcQ = vQRST;
            nn = length(QRSTcQ);
        end
        
        % first beat inserted at rr(2)-Rind
        if beatNr == 1
            if lead == 1
                rIndex = rr(2);
                TendIndex = size(ecgSig,2)+length(QRSTc);
            end
            ecgSig = zeros(1,rr(2)-Rind+nn+1);
            ecgSig(end-nn+1:end) = ecgSig(end-nn+1:end) + QRSTcQ;
        else
            if lead == 1
                %QIndex = [QIndex data.];
                TendIndex = [TendIndex size(ecgSig,2)+length(QRSTc)];
                %rIndex = [rIndex (rIndex(1,end) + TQlength0 + length(QRSTc) - Rind + RindNext)];
                if (state_history(beatNr-1)==state_history(beatNr))&&(state_history(beatNr)==5)
                    if rr(beatNr+1)<vpb_end
                        rr(beatNr+1) = vpb_end;
                    end
                end
                rIndex = [rIndex (rIndex(1,end) + rr(beatNr+1) )];
            end
            % adding extra samples to ecgSig to allow for partial
            % superposition between ecgSig and QRSTcQ
            A = length(ecgSig);
            B = rr(beatNr+1);
            C = A - rIndex(beatNr-1);
            D = B - C;
            E = length(QRSTcQ) - Rind;
            if (D+E) > 0
                newSamples = zeros(1,D+E);
                ecgSig = [ecgSig newSamples];
                ecgSig(end-nn+1:end) = ecgSig(end-nn+1:end) + QRSTcQ;
            else
                %beat should be place a bit before end of ecgSig
                d = D+E;
                ecgSig(end-nn+d+1:end+d) = ecgSig(end-nn+d+1:end+d) + QRSTcQ;
            end
            if ecgSig(rIndex(beatNr))<0.5
                debug = 1;
            end
        end
    end
    
    nn = length(ecgSig);
    ecgLength = length(pqrstResampled(1,:));
    if nn > ecgLength
        pqrstResampled = [pqrstResampled zeros(numLeads,nn-ecgLength)];
    end
    pqrstResampled(lead,1:nn) = ecgSig; %% Changed in v012020
end
rr = rr(2:end);
QRSindex = rIndex;
ecgLength = length(pqrstResampled(1,:));

%% Respiratory influence for altering QRST morphology (Changed in v012020)
%Last update: CPerez 05/2022

T = ecgLength/Fs;
Zr = [5 5 5]; % Maximum angular variation around X,Y,Z axis

switch realVAon %%
    case 0
        %CPerez 05/2022
        %Include respiration
        [Q] = simECG_rotation_angles_XYZ_v2( simECGdata.Fr, T, Fs, Zr);
        
        for t = 1 : ecgLength
            pqrstResampled(:,t) = Q(:,:,t)*pqrstResampled(:,t);
        end
        
        % Dower Transformation for synthetic PQRST complexes
        %Dower transform
        D = [
            -0.515  0.157   -0.917;
            0.044  0.164   -1.387;
            0.882  0.098   -1.277;
            1.213  0.127   -0.601;
            1.125  0.127   -0.086;
            0.831  0.076    0.230;
            0.632 -0.235    0.059;
            0.235  1.066   -0.132];
        
        pqrst_v1_v6_I_II = D*pqrstResampled;
        
        multileadVA(1:8,:) = pqrst_v1_v6_I_II(1:8,:); % V1, V2, V3, V4, V5, V6, I, II,
        multileadVA(9,:) = pqrst_v1_v6_I_II(8,:) - pqrst_v1_v6_I_II(7,:);     % III
        multileadVA(10,:) = -(pqrst_v1_v6_I_II(7,:) + pqrst_v1_v6_I_II(8,:))/2; % aVR
        multileadVA(11,:) =  pqrst_v1_v6_I_II(7,:) - pqrst_v1_v6_I_II(8,:)/2;   % aVL
        multileadVA(12,:) =  pqrst_v1_v6_I_II(8,:) - pqrst_v1_v6_I_II(7,:)/2;   % aVF
        
        multileadVA(13,:) = pqrstResampled(1,:); % X
        multileadVA(14,:) = pqrstResampled(2,:); % Y
        multileadVA(15,:) = pqrstResampled(3,:); % Z
        aux = multileadVA;
        multileadVA(1:6,:) = multileadVA(7:12,:);
        multileadVA(7:12,:) = aux(1:6,:);
        
    case 1 % 12-lead ECG was available
        multileadVA = pqrstResampled;
        
end

%% Ventricular activity
if Nbvt
    
    nVPB = 0;
    for k = 1:length(v_idx)
        
        nVPB = nVPB + 1;
        QRSTtemp = squeeze(vpb_15l(kv_history(nVPB),:,:));
        QRSTtemp = QRSTtemp * (0.8 + rand*0.2);
        
        if realVAon
            for l=1:15
                QRSTtemp(l,:) = simECG_correct_baseline(QRSTtemp(l,:));
            end
        end
        
        multileadVA(:,rIndex(v_idx(k))-vpb_R+1:rIndex(v_idx(k))-vpb_R+Lv) = multileadVA(:,rIndex(v_idx(k))-vpb_R+1:rIndex(v_idx(k))-vpb_R+Lv) + QRSTtemp;
        
    end
    
end
% store rr intervals in simECGdata
simECGdata.rr = rr./simECGdata.fs;

end

