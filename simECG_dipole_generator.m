function [DIP, teta]= simECG_dipole_generator(N,fs,f,alphai,bi,tetai,teta0)
% []= simECG_dipole_generator() Synthetic cardiac dipole generator using
% the 'differential form' of the dipole equations. For further details see
% a paper by Sameni et al. Multichannel ECG and noise modeling: application
% to maternal and fetal ECG signals. EURASIP Journal on Advances in Signal
% Processing.
%
% Input arguments:
% N - signal length, in samples.
% fs - sampling rate.
% f - heart rate, in Hz. 
% alphai - structure contaning the amplitudes of Gaussian functions used for 
%       modeling the x, y, and z coordinates of the cardiac dipole.
% bi - structure contaning the widths of Gaussian functions used for 
%       modeling the x, y, and z coordinates of the cardiac dipole.
% tetai - structure contaning the phase of Gaussian functions used for 
%       modeling the x, y, and z coordinates of the cardiac dipole.
% teta0 - initial phase of the synthetic dipole .
% 
% Output arguments:
% DIP - structure contaning the x, y, and z coordinates of the cardiac
% dipole. 
% teta - vector containing the dipole phase.
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html 

w = 2*pi*f; 
dt = 1/fs; 
 
teta = zeros(1,N); 
X = zeros(1,N); 
Y = zeros(1,N); 
Z = zeros(1,N); 
 
teta(1) = teta0; 
for i = 1:N-1
    teta(i+1) = mod(teta(i) + w*dt + pi , 2*pi) - pi; 
 
    dtetaix = mod(teta(i) - tetai.x + pi , 2*pi) - pi; 
    dtetaiy = mod(teta(i) - tetai.y + pi , 2*pi) - pi; 
    dtetaiz = mod(teta(i) - tetai.z + pi , 2*pi) - pi; 
 
    if(i==1)
        X(i) = sum(alphai.x .* exp(-dtetaix .^2 ./ (2*bi.x .^ 2))); 
        Y(i) = sum(alphai.y .* exp(-dtetaiy .^2 ./ (2*bi.y .^ 2))); 
        Z(i) = sum(alphai.z .* exp(-dtetaiz .^2 ./ (2*bi.z .^ 2))); 
    end 
 
    X(i+1) = X(i) - dt*sum(w*alphai.x ./ (bi.x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2* bi.x .^ 2)));   % x state variable 
    Y(i+1) = Y(i) - dt*sum(w*alphai.y ./ (bi.y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2* bi.y .^ 2)));   % y state variable 
    Z(i+1) = Z(i) - dt*sum(w*alphai.z ./ (bi.z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2* bi.z .^ 2)));   % z state variable 
 
end 

 
DIP.x = X; 
DIP.y = Y; 
DIP.z = Z;