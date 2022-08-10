function walkP = RandomWalk(nSteps,p,movMod, movAng)
%nSteps is number of steps
%p is the pole with real and imaginary part
%movMod is the increment/decrement of radious
%movAng is the increment/decrement of angle

% Copyright (c), Cristina Perez, University of Zaragoza, 08/2022

rad = abs(p);
ang = angle(p);
sig = 0;

if ang >= pi/2
    sig = 1;
    ang = ang - (pi/2);
end

%Random walk model for Radious
sum = rad;
newRad(1) = sum;
prob = rand(nSteps,1);

for ii = 2:nSteps
    if (prob(ii) > 0.5 && (sum + movMod)<1) || (prob(ii) <= 0.5 && (sum - movMod)<0)
        sum = sum + movMod;
    else
        sum = sum - movMod;
    end
    newRad(ii) = sum;
end

%Random walk model for angle
sum = ang;
newAng(1) = sum;
prob = rand(nSteps,1);

for ii = 2:nSteps
    if (prob(ii) > 0.5 && (sum + movAng)<(pi/2)) || (prob(ii) <= 0.5 && (sum - movAng)<0)
        sum = sum + movAng;
    else
        sum = sum - movAng;
    end
    newAng(ii) = sum;
end

%Polar to complex of all walk
if sig; newAng = newAng+pi/2; end

for ii = 1:nSteps
    walkP(ii,1) = newRad(ii)*exp(newAng(ii)*i);
end

end