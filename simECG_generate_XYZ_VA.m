function [pqrstFrank, Qind, Rind, Sind, Tind] = simECG_generate_XYZ_VA(num)
%
% [] = simECG_generate_XYZ_VA() returns simulated QRST complexes.
% Ventricular activity is simulated by using a dynamic ECG model
% (Sameni et al. 2007) which is an extended version of a single channel
% ECG simulator originally proposed in (McSharry et al. 2003). By using this
% model, three orthogonal Frank leads of N samples length are generated.
% Ventricular activity is generated in 3 Frank leads

N = 1000;                       % Signal length in samples
fs = 1000;                      % sampling rate
pqrstFrank = zeros(3,num,540); %560 default;alba:650

% Dipole parameters
F = 1;                          % heart rate
teta0 = -pi/2;                  % initial phase of the ECG

%% QRS
Qw = simECG_random_number(0.05, 0.08);
Rw = simECG_random_number(0.05, 0.08);
Sw = simECG_random_number(0.05, 0.08);

% Increasing QRS complex maximum amplitude
A = 1.5; 
f = rand * A;
asf=zeros(1,6);
for k=1:6
    %amplitude scaling factor
    asf(k) = max(min(f + (1*(rand-0.5)),A),0); % 50% bounded variation around selected value
end

alphaiQRS.x = [ -0.05+((-0.4+0.05)*asf(1))      0.4+((1.5-0.4)*asf(2))      0 ];
biQRS.x     = [Qw    Rw    Sw];
tetaiQRS.x  = [-0.1    0    0.1];

alphaiQRS.y = [ 0                               0.1+((0.7-0.1)*asf(3))      -0.05+((-0.3+0.05)*asf(4)) ];
biQRS.y     = [Qw   Rw   Sw];
tetaiQRS.y  = [-0.1   0   0.1];

alphaiQRS.z = [ -0.05+((-0.4+0.05)*asf(5))     0                            0.1+((1-0.1)*asf(6)) ];
biQRS.z     = [Qw   Rw  Sw];
tetaiQRS.z  = [ -0.1   0   0.1];

%% T wave
Tw = simECG_random_number(0.5, 0.7);

% Txa = simECG_random_number(0.02, 0.12);
% Tya = simECG_random_number(0.01, 0.05);
% Tza = simECG_random_number(-0.02, -0.1);
Txa = simECG_random_number(0.02, 0.08);
Tya = simECG_random_number(0.01, 0.03);
Tza = simECG_random_number(-0.02, -0.06);

alphaiT.x = [Txa   2*Txa   3*Txa];
biT.x     = [Tw     Tw/2    Tw/4];
tetaiT.x  = [1.1     1.4    1.6];

alphaiT.y = [Tya   2*Tya   3*Tya];
biT.y     = [Tw     Tw/2    Tw/4];
tetaiT.y  = [1.1     1.4    1.6];

alphaiT.z = [Tza   2*Tza   3*Tza];
biT.z     = [Tw     Tw/2    Tw/4];
tetaiT.z  = [1.1     1.2    1.4];
%%
alphai.x = [alphaiQRS.x alphaiT.x];
alphai.y = [alphaiQRS.y alphaiT.y];
alphai.z = [alphaiQRS.z alphaiT.z];

bi.x = [biQRS.x biT.x];
bi.y = [biQRS.y biT.y];
bi.z = [biQRS.z biT.z];

tetai.x = [tetaiQRS.x tetaiT.x];
tetai.y = [tetaiQRS.y tetaiT.y];
tetai.z = [tetaiQRS.z tetaiT.z];

% Generate phase for each variability related frequency
if rand(1,1) > 0.5
    phaseX = 0;
else
    phaseX = pi;
end

if rand(1,1) > 0.5
    phaseY = 0;
else
    phaseY = pi;
end

if rand(1,1) > 0.5
    phaseZ = 0;
else
    phaseZ = pi;
end

[DIP, ~] = simECG_dipole_generator(N,fs,F,alphai,bi,tetai,teta0); 
for n = 1:num
    pqrstFrank(1:3,n,:) = [DIP.x(51:590) ; DIP.y(51:590) ; DIP.z(51:590)]; %default 610, Alba: 6000
end

[~,qrsMaxInd] = max(pqrstFrank(1,1,1:300),[],3);

Qind = qrsMaxInd-50;
Rind = qrsMaxInd;
Sind = qrsMaxInd+50;
Tind = size(pqrstFrank,3);
end

