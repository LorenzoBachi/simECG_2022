function walk = RandomWalk(nSteps,center,movVal)
%nSteps is number of steps
%center is the value of real (imaginary) part
%movVal is the increment/decrement

% Copyright (c), Cristina Perez, University of Zaragoza, 08/2022

sig = sign(center);
center = abs(center);

sum = center;
walk(1) = sum;
prob = rand(nSteps,1);

for ii = 2:nSteps
    if (prob(ii) > 0.5 && (sum + movVal)<1) || (prob(ii) <= 0.5 && (sum - movVal)<0)
        sum = sum + movVal;
    else
        sum = sum - movVal;
    end
    walk(ii) = sum;
end
walk = walk.*sig;
end