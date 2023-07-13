function [ Q, phiX, phiY, phiZ] = simECG_rotation_angles_XYZ_v2( Fr, T, Fs, Zr)
% [] = simECG_rotation_angles_XYZ_v2() calculates the rotation matrix Q
% with time-varying angles (phiX(n), phiY(n), phiZ(n)).
% 
% Input arguments:
% Fr - respiratory Pattern, in Hz.
% Tr - time vector for Fr, in seconds.
% T  - length of ECG, in seconds.
% Fs - sampling Frequency, in Hz.
% Zr - maximum angular variation, in rad/sec.
% 
% Output arguments:
% Q - the rotation matrix.
% phiX, phiY, phiZ - the time varying angles.
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

%Last update: CPerez 05/2022

% Initialize the parameters of the sigmoidal functions
phiX = []; phiY = []; phiZ = [];

for nd = 1:numel(Fr)
    %Template respiratory cycle: random variables with a std of 10%
    sigIn = 3; sigEx = 3; deltIn = 2; deltEx = 7; Tr = 10;

    deltIn = 0.9*deltIn + (deltIn*(1.1-0.9)*rand(1));
    deltEx = 0.9*deltEx + (deltEx*(1.1-0.9)*rand(1));
    
    nTemp = 0:1/Fs:Tr;
    sigmoid_InTemp = 1./(1 + exp( -sigIn* (nTemp - deltIn)));
    sigmoid_ExTemp = 1./(1 + exp( sigEx* (nTemp - deltEx)));
    cycle_Temp = sigmoid_InTemp.*sigmoid_ExTemp;
    
%     % Rotation angles around X axis
%     Z0 = (1/(1+exp(-(deltIn - deltEx)/(1/(-sigIn) - 1/sigEx))))^2;
    
    % Amplitude modulation in the orthogonal leads 
    A = rand(1,3)*(1.1-0.9) + 0.9; %Cris 2022
    cycle = resample(cycle_Temp,round(Fs./Fr(nd)),Tr*Fs);
    
%     phiX = [phiX A(1)*(1/Z0).*cycle(1:round(Fs./Fr(nd)))];
    phiX = [phiX A(1).*cycle(1:round(Fs./Fr(nd)))];
    phiY = [phiY A(2).*cycle(1:round(Fs./Fr(nd)))];
    phiZ = [phiZ A(3).*cycle(1:round(Fs./Fr(nd)))];
end
phiX = Zr(1)*phiX(1:T*Fs); 
phiY = Zr(2)*phiY(1:T*Fs); 
phiZ = Zr(3)*phiZ(1:T*Fs);


h = pi/180;% Rotation angle in radians

for t = 1 : length(phiX)
    Qx = [1 0 0; 0 cos(h*phiX(t)) sin(h*phiX(t)); 0 -sin(h*phiX(t)) cos(h*phiX(t))];
    Qy = [cos(h*phiY(t)) 0 sin(h*phiY(t)); 0 1 0 ; -sin(h*phiY(t)) 0 cos(h*phiY(t))];
    %     Qz = [cos(h*phiZ(t)) sin(h*phiZ(t)) 0 ; -sin(h*phiZ(t)) 0 cos(h*phiZ(t)) ; 0 0 1];
    Qz = [cos(h*phiZ(t)) -sin(h*phiZ(t)) 0 ; sin(h*phiZ(t)) cos(h*phiZ(t)) 0; 0 0 1];
    Q(:,:,t) = Qx*Qy*Qz;
end

end

