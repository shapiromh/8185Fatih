function f = CHEBY(a,b,c,m,x)

% CHEBY Calculates m order Chebyshev polynomial approximation for input based on coefficients

% a and b are the same max and min on the grid for CHEBYC
% c will be the first m < n coefficients generated from CHEBYC
% m is the order of the polynomial
% x is the point on which the approximation will be evaluated

% Error checking

if ((x-a)*(x-b)>0)
    display 'X not in range [a,b]!'
    return
end

d=0;
dd=0;
y=(2*x-a-b)/(b-a); % Change variable to point on [-1,1]
y2=2*y;

% Calculate the function at y

for j=m:-1:2
    sv=d;
    d=y2*d-dd+c(j);
    dd=sv;
end
f=y*d-dd+.5*c(1);
