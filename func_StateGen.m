function [states] = simECG_StateGen(l2h, h2l)
%% Generation of states
% l2h is the probability of transition to noisy segment 
% h2l is the probability of transition to clean segment 
% The function generates a markov chain with two states with 30 samples,
% which then upsampled by the factor of 1000 to match with a 30s thumb ECG with sampling frequency of 1000 Hz. 
% Thresholding is done to remove the ringing that otherwise would have occured. 
% Hesam Halvaei, June 01 2022, Lund University
%%
N = 30;
P = ([1-l2h, l2h; h2l, 1-h2l]); % transition matrix
mc = dtmc(P);
states = simulate(mc, N);
states = interp(states, 1000);

states(states > 1.5) = 2;
states(states < 1.5) = 1;

states(states == 1) = 0;
states(states == 2) = 1;
end

