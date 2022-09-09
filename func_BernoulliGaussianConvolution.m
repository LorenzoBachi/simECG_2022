function [bernogauss_conv] = func_BernoulliGaussianConvolution(len, bernogauss)
% func_BernoulliGaussianConvolution convolve the a realizaiton of bernoulli
% gaussian process with varying filter, i.e., 
% INPUT a bernoulli gaussian realization
% OUTPUT convolved bernoulli gaussian realization

inds = find(bernogauss ~= 0); % finds the non-zero indices from the process
bernogauss_conv = 0;
a1 = 0.95;
a2 = 0.99;
K = randi([250, 1250]);
N = 4000;
n = 1:K-1; 
h(n) = a1.^(-n);
n = K:N;
h(n) = a1.^(-(K-1)) * (a2.^(n-K));
h = h/(max(abs(h)));
if inds ~= 0
    for q = 1:length(inds)
        temp = zeros(len, 1);
        temp(inds(q)) = bernogauss(inds(q));
        temp = conv(temp, h, 'same');
        bernogauss_conv  = bernogauss_conv + temp;
    end
end

end

