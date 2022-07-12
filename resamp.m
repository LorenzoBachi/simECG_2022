function [ntimes,y,data]=resamp(x,times,freq,nfreq,time_segments,n_max_interp,InterpYes)

%Function
% [ntimes,y]=resamp(x,times,freq,nfreq,time_segments)
%
%Purpose
%  resamples a non-even sampled vector 
%
%Synopsis
% function resamp(x,times,freq,nfreq,time_segments)
%
%Description
%  resamples a non-even sampled series vector (RR, QT, KL, ...) with occurrence
%  times vector 'times' referred to a 'freq' sampling frequency into an even-sampled
%  vector with sampling frequency 'nfreq' Hz. It outputs the new even-sampled vector
%  and the time occurrences. It also works with an input matrix (L,no_leads).
%  time_segments = [t1 t2; t3 t4;...]
%  n_max_interp = maximum time between two samples (in seconds)
%
% Copyright (c), Jose Garcia, Zaragoza University, Spain 1997-08-08
%--------------------------------------------------------------------------------

% ----- series resampling at nfreq Hz ----------------------------------------
data = [];
if nargin < 5
    time_segments = [];
end
if nargin < 6
    n_max_interp = [];
end
[L,no_leads]=size(x);
ntimes=(times(1):freq/nfreq:times(L))';
if ~isempty(time_segments)   
    ntimes_all = [];
    for i = 1:size(time_segments,1)
        ntimes_new = ntimes;
        indx_f = find(ntimes_new > round(time_segments(i,2))); % Juan 16/07/2011
        indx_i = find(ntimes_new < round(time_segments(i,1))); % Juan 16/07/2011
        indx = [indx_i ;indx_f];% Juan 16/07/2011
        ntimes_new(indx) = [];%
        ntimes_all = [ntimes_all; ntimes_new];
    end
    ntimes = ntimes_all;
elseif ~isempty(n_max_interp)
    indxs = find(isnan(x));
    x(indxs) = [];
    times(indxs) = [];
    indxs = find(diff(times) >= n_max_interp*1000); 
    x([indxs+1 indxs+2]) = NaN;
%     indxs = find(diff(times) > n_max_interp*1000);   
%     if ~isempty(indxs)
%         for i = 1:length(indxs)+1
%             if i == 1
%                 data(i,:) = [1 times(indxs(i))];            
%             else
%                 try
%                     data(i,:) = [times(indxs(i-1)+1) times(indxs(i))];
%                 catch
%                     data(i,:) = [times(indxs(i-1)+1) times(end)];
%                 end
%             end
%         end
%         if InterpYes
%             x([indxs(i)+1 indxs(i)+2]) = NaN;
%         else
%             ntimes_all = [];
%             for i = 1:size(data,1)
%                 ntimes_new = ntimes;
%                 indx_f = find(ntimes_new > round(data(i,2))); % Juan 16/07/2011
%                 indx_i = find(ntimes_new < round(data(i,1))); % Juan 16/07/2011
%                 indx = [indx_i ;indx_f];% Juan 16/07/2011
%                 ntimes_new(indx) = [];%
%                 ntimes_all = [ntimes_all; ntimes_new];
%             end
%             ntimes = ntimes_all;
%         end
%     else
%         data = [times(1) times(end)];
%     end
end
y=zeros(length(ntimes),no_leads);
y = interp1(times,x,ntimes,'spline');
% y=interp1q(times,x,ntimes);
% aux=isnan(y);
% y(aux)=0;  