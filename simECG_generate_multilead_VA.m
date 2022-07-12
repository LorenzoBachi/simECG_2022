function [QRSindex, TendIndex, rr, multileadVA, ecgLength] = simECG_generate_multilead_VA(rrLength, targets_beats, rr, realVAon, realAAon, realRRon, ecgParameters,state_history)
% [] = simECG_gen_multilead_VA() returns multilead (15 lead) ventricular
% activity. A set of 100 15-lead ECGs with SR selected from the PTB Diagnostic
% ECG Database is used as a basis for modeling ventricular activity. The ECGs
% of the PTB database are first subjected to baseline removal and QRST delineation.
% The original T-waves are then resampled to a fixed width and, depending on
% the type of rhythm, width-adjusted to match prevailing heart rate. Since
% the original ECGs last just for about 2 min, QRST complexes are subjected
% to repeated concatenation until desired length of ECG is obtained. The TQ
% interval is interpolated using a cubic spline interpolation.
%
%Generated leads:
% multileadVA(1,:) - I        multileadVA(7,:) - V1      multileadVA(13,:) - X
% multileadVA(2,:) - II       multileadVA(8,:) - V2      multileadVA(14,:) - Y
% multileadVA(3,:) - III      multileadVA(9,:) - V3      multileadVA(15,:) - Z
% multileadVA(4,:) - aVR      multileadVA(10,:) - V4
% multileadVA(5,:) - aVL      multileadVA(11,:) - V5
% multileadVA(6,:) - aVF      multileadVA(12,:) - V6

disp('Generating ventricular activity ...');

rr = (round(rr*1000));    % rr intervals in miliseconds
rr = [rr(1) rr]';
Fs = 1000;                      % sampling rate
rr_sec = rr/Fs;

%% QT interval adaptation
nfreq = 4; %Interpolation to 4 Hz;
timebeats = cumsum(rr(2:end))/1000; %In seconds
timebeats = [0;timebeats];
% resample RR series to an evenly spaced series
%[pos,RR_ori]=resamp(rr_sec,timebeats*Fs,Fs,nfreq);
pos = (timebeats(1)*Fs:Fs/nfreq:timebeats(length(rr_sec))*Fs)';
RR_ori = interp1(timebeats,rr_sec,pos,'spline'); %pos in samples, RR_ori in sec
pos = pos/Fs; %sec
RR_exp = NaN(size(RR_ori));
L = 300; %window size in seconds. Note: this parameter is used when real EST RR is loaded
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


%% Normal PQRST complexes generation
switch realVAon
    case 0 % Generate synthetic PQRST complexes
        %[data.pqrst, data.Qind, data.Rind, data.Sind, data.Tind] = simPAF_gen_syn_VA(100);
        [data.pqrst, data.Qind, data.Rind, data.Sind, data.Tind] = simECG_generate_XYZ_VA(100);%% Changed in v012020
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

%% Ectopic PQRST complexes generation
% loading ventricular complexes (only real) from the contribution of
% Alcaraz (2019)
load('DATA_ventricular_Frank');
% randomly select one of the 15 possible ectopic beat collections
k = randi([1 15]);
vpb = DATAventricularFrank(k).pqrst;
vpb_R = DATAventricularFrank(k).R;
% number of the ventricular beats of the selected collection
Nvp = length(vpb(1,:,1));
% number of ventricular beats that are present in the record
Nbvt = sum(targets_beats==4);
% amplitude factors for each VPB
vpb_amp = zeros(1,Nbvt);
f = (1 + rand*1.2) * ( max(sum(abs(squeeze(mean(pqrst,2))))) / max(sum(abs(squeeze(mean(vpb,2))))) );
for k = 1:Nbvt
    vpb_amp(k) =  f * (0.9 + rand*0.2);
    %(rand*0.25+0.25); %Alcaraz: (randi(10*[0.8 2.5])/10) * sign(1.3+randn(1,1)) ;
end
% if necessary, replicate the beats to reach the required number
if Nbvt > Nvp
    % number of required replications
    rep = ceil(Nbvt/Nvp);
    for k = 1 : rep
        vpb = [vpb vpb];
        vpb_R = [vpb_R vpb_R];
    end
end

%% Activity generation

pqrstResampled=zeros(numLeads,1);
ecgSig = [];
%Load QRST complexes for all leads
for lead = 1:numLeads %% Changed in v012020
    AFend = 0;
    % counters for supraventricular and ventricular beat numbers
    ks = 0;
    kv = 0;
    
    switch realAAon
        case 0 % Synthetic atrial activity
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
                vQRST = squeeze(vpb(lead,1,:))';
                vQRST = vpb_amp(1) * vQRST;
                vQRST = simECG_correct_baseline(vQRST);
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
                        vQRST = squeeze(vpb(lead,kv+1,:))';
                        vQRST = vpb_amp(kv+1) * vQRST;
                        vQRST = simECG_correct_baseline(vQRST);
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
                    
                    % Gradual prolongation of T wave after AF is terminated
                    if (targets_beats(beatNr) == 1)||(targets_beats(beatNr) == 3)
                        if AFend == 1
                            %                         ST_length =
                            %                         round(length(QRST_ST)*sqrt((0.3+m*0.1)*rr(beatNr)/1000));
                            %                         %Original version Andrius
                            ST_length = qt_hyper(beatNr); % Alba: Hyperbolic model QT/RR regression
                            m = m + 1;
                            if m == 7
                                AFend = 0;
                            end
                        else
                            %                         ST_length =
                            %                         round(length(QRST_ST)*sqrt(rr(beatNr)/1000)); %v0
                            %                         Andrius
                            ST_length = qt_hyper(beatNr); % Alba: Hyperbolic model QT/RR regression
                        end
                    else
                        % ST_length = round(length(QRST_ST)*sqrt(0.35)); %Original code Andrius
                        %                     ST_length = round(QTc*sqrt(rr(beatNr)/1000)); % Alba: Bazzet's refers to the whole QT interval
                        ST_length = qt_hyper(beatNr); % Alba: Hyperbolic model QT/RR regression
                        AFend = 1;
                        m = 1;
                    end
                    %                 QRST_ST = resample(QRST_ST, ST_length+round(20*sqrt(rr(beatNr)/1000))-length(QRST_PS), length(QRST_ST));  % T wave length correction QT = QTc*sqrt(RR) ~ T = Tc*sqrt(RR) %Bazzet's
                    border = round(10*(ST_length-length(QRST_PS))/(length(QRST_ST)-20));
                    QRST_ST = resample(QRST_ST, ST_length-length(QRST_PS), length(QRST_ST)-20);  % T wave length correction QT = QTc*sqrt(RR) ~ T = Tc*sqrt(RR)
                    QRSTc = [QRST_PS QRST_ST(border:end-border)]; % Alba: adapted QRST
                    %QRSTc = QRST_ST(round(20*sqrt(rr(beatNr)/1000))+1:end-round(20*sqrt(rr(beatNr)/1000))); %Alba: it is not the QTc but the original QT interva before correction
                    
                else
                    % Ventricular beat
                    % R wave fiducial point
                    Rind = vpb_R(kv+1);
                    % in bigeminy / trigeminy, the same VPB si used
                    if state_history(beatNr)~=4
                        kv = kv + 1;
                    else
                        QRSTtemp = QRSTtemp * (0.9 + rand*0.2);
                    end
                    % QT correction for VPBs turned off
                    QRSTc = QRSTtemp;
                end
                
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
                        %                     QIndex = [QIndex data.];
                        TendIndex = [TendIndex size(ecgSig,2)+length(QRSTc)];
                        %rIndex = [rIndex (rIndex(1,end) + TQlength0 + length(QRSTc) - Rind + RindNext)];
                        rIndex = [rIndex (rIndex(1,end) + rr(beatNr+1) )];
                    end
                    % adding extra samples to ecgSig to allow for partial
                    % compenetration between ecgSig and QRSTcQ (not supported
                    % in previous versions)
                    A = length(ecgSig);
                    B = rr(beatNr+1);
                    C = A - rIndex(1,beatNr-1);
                    D = B - C;
                    E = length(QRSTcQ) - Rind;
                    if (D+E) > 0
                        newSamples = zeros(1,D+E);
                        ecgSig = [ecgSig newSamples];
                    end
                    ecgSig(end-nn+1:end) = ecgSig(end-nn+1:end) + QRSTcQ;
                end
            end
            
        case 1 % Real atrial activity
            %{
            if AFburden == 1 % Entire signal is AF
                rIndex = data.Rind - data.Qind;
                TendIndex =[];
            else
                rIndex = data.Rind;
                TendIndex =[];
            end
            %}
            for beatNr = 1:rrLength
                if (targets_beats(beatNr) == 1)||(targets_beats(beatNr) == 3)
                    %{
                        if (targets_beats(beatNr-1) == 1)||(targets_beats(beatNr-1) == 3)
                            QRSTtemp = pqrst(lead,beatNr,1:end);
                            QRSTtemp = simECG_correct_baseline(QRSTtemp);
                            Qind = data.Qind; Rind = data.Rind;
                            Sind = data.Sind; Tind = data.Tind;
                        else % targets_SR_AF(beatNr-1) == 2
                            QRSTtemp = pqrst(lead,beatNr,Qind+1:end);
                            QRSTtemp = simECG_correct_baseline(QRSTtemp);
                            Qind = data.Qind; Rind = data.Rind - data.Qind;
                            Sind = data.Sind - data.Qind; Tind = data.Tind - data.Qind;
                        end
                    %}
                    Qind = data.Qind; Rind = data.Rind;
                    Sind = data.Sind; Tind = data.Tind;
                    QRSTtemp = pqrst(lead,beatNr,1:end);
                    QRSTtemp = simECG_correct_baseline(QRSTtemp);
                else % targets_SR_AF(beatNr) == 2
                    %{
                        if (targets_beats(beatNr-1) == 1)||(targets_beats(beatNr-1) == 3)
                            QRSTtemp = pqrst(lead,beatNr,Qind+1:end);
                            QRSTtemp = simECG_correct_baseline(QRSTtemp);
                            Qind = data.Qind; Rind = data.Rind - data.Qind;
                            Sind = data.Sind - data.Qind; Tind = data.Tind - data.Qind;
                        else % targets_SR_AF(beatNr-1) == 2
                            QRSTtemp = pqrst(lead,beatNr,1:end);
                            QRSTtemp = simECG_correct_baseline(QRSTtemp);
                            Qind = data.Qind; Rind = data.Rind;
                            Sind = data.Sind; Tind = data.Tind;
                        end
                    %}
                    Qind = data.Qind; Rind = data.Rind - data.Qind;
                    Sind = data.Sind - data.Qind; Tind = data.Tind - data.Qind;
                    QRSTtemp = pqrst(lead,beatNr,Qind+1:end);
                    QRSTtemp = simECG_correct_baseline(QRSTtemp);
                end
                
                if beatNr == length(targets_beats)
                    if (targets_beats(beatNr) == 1)||(targets_beats(beatNr) == 3)
                        %RindNext = data.Rind;
                        QRSTtempNext = pqrst(lead,beatNr,1:end);
                        QRSTtempNext = simECG_correct_baseline(QRSTtempNext);
                    else
                        %RindNext = data.Rind-data.Qind;
                        QRSTtempNext = pqrst(lead,beatNr,Qind+1:end);
                        QRSTtempNext = simECG_correct_baseline(QRSTtempNext);
                    end
                else
                    if (targets_beats(beatNr+1) == 1)||(targets_beats(beatNr+1) == 3)
                        if targets_beats(beatNr) == 2
                            %RindNext = data.Rind-data.Qind;
                            QRSTtempNext = pqrst(lead,beatNr+1,Qind+1:end);
                            QRSTtempNext = simECG_correct_baseline(QRSTtempNext);
                        else
                            %RindNext = data.Rind;
                            QRSTtempNext = pqrst(lead,beatNr+1,1:end);
                            QRSTtempNext = simECG_correct_baseline(QRSTtempNext);
                        end
                    else
                        if (targets_beats(beatNr) == 1)||(targets_beats(beatNr) == 3)
                            %RindNext = data.Rind;
                            QRSTtempNext = pqrst(lead,beatNr+1,1:end);
                            QRSTtempNext = simECG_correct_baseline(QRSTtempNext);
                        else
                            %RindNext = data.Rind-data.Qind;
                            QRSTtempNext = pqrst(lead,beatNr+1,Qind+1:end);
                            QRSTtempNext = simECG_correct_baseline(QRSTtempNext);
                        end
                    end
                end
                
                % Divide PQRST into two parts
                QRST_PS = QRSTtemp(1,1:Sind); % The PQRS part
                QRST_ST = QRSTtemp(1,Sind+1:Tind); % The ST part
                % Resample T wave according to current RR interval
                % Signal is prolonged in order to avoid boundary effects of resampling
                QRST_ST  = [QRST_ST(1,1)*ones(1,20) QRST_ST QRST_ST(1,end)*ones(1,20)];
                
                % Gradual prolongation of T wave after AF is terminated
                if (targets_beats(beatNr) == 1)||(targets_beats(beatNr) == 3)
                    if AFend == 1
                        %                         ST_length = round(length(QRST_ST)*sqrt((0.3+m*0.1)*rr(beatNr)/1000));
                        %                         m = m + 1;
                        %                         if m == 7
                        %                             AFend = 0;
                        %                         end
                        ST_length = qt_hyper(beatNr); % Alba: Hyperbolic model QT/RR regression
                    else
                        %                         ST_length = round(length(QRST_ST)*sqrt(rr(beatNr)/1000));
                        ST_length = qt_hyper(beatNr); % Alba: Hyperbolic model QT/RR regression
                        
                    end
                else
                    %                     ST_length = round(length(QRST_ST)*sqrt(0.35));
                    ST_length = qt_hyper(beatNr); % Alba: Hyperbolic model QT/RR regression
                    AFend = 1;
                    m = 1;
                end
                QRST_ST = resample(QRST_ST, ST_length+40, length(QRST_ST));  % T wave length correction QT = QTc*sqrt(RR) ~ T = Tc*sqrt(RR)
                QRSTc = [QRST_PS QRST_ST(21:end-20)]; % Corrected PQRST
                % Find TQ length
                TQlength = rr(beatNr+1) - length(QRSTc); % + Rind - RindNext;
                % Protect against the error of to low heart rate
                if  TQlength < 2
                    TQlength = 2;
                end
                % Interpolate TQ interval
                x = [1 TQlength];
                xi = 1:1:TQlength;
                y = [QRSTc(1,end) QRSTtempNext(1, 1)];
                TQ = interp1(x,y,xi,'linear');
                QRSTcQ = [ QRSTc TQ ]; % Lorenzo
                nn = length(QRSTcQ);
                
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
                        %                     QIndex = [QIndex data.];
                        TendIndex = [TendIndex size(ecgSig,2)+length(QRSTc)];
                        %rIndex = [rIndex (rIndex(1,end) + TQlength0 + length(QRSTc) - Rind + RindNext)];
                        rIndex = [rIndex (rIndex(1,end) + rr(beatNr+1) )];
                    end
                    % adding extra samples to ecgSig to allow for partial
                    % compenetration between ecgSig and QRSTcQ (not supported
                    % in previous versions)
                    A = length(ecgSig);
                    B = rr(beatNr+1);
                    C = A - rIndex(1,beatNr-1);
                    D = B - C;
                    E = length(QRSTcQ) - Rind;
                    if (D+E) > 0
                        newSamples = zeros(1,D+E);
                        ecgSig = [ecgSig newSamples];
                    end
                    ecgSig(end-nn+1:end) = ecgSig(end-nn+1:end) + QRSTcQ;
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

if ~realRRon %CPerez 05/2022
        %Include respiration
        [Q] = simECG_rotation_angles_XYZ_v2( ecgParameters.Fr, T, Fs, Zr, ecgParameters);

        for t = 1 : ecgLength
            pqrstResampled(:,t) = Q(:,:,t)*pqrstResampled(:,t);
        end
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


end

% function from Arcaraz used to fit the ventricular beat in the simulated
% record
function y = baseline(x)
t=[1 length(x)];
degree = 1;
[p, S, mu] = polyfit(t, x(t), degree);
xp = polyval(p, 1:length(x), [], mu);
y = x - xp;
end
