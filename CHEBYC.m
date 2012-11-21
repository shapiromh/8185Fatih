function f = CHEBYC(a,b,n,func)

% CHEBY Calculates n coefficients for Chebyshev polynomial approximation of a function

% a and b are upper and lower limits on a grid with n points
% func is the function; note for interpolation with func unknown we need 
% a condition for the func (think VF in DP problem). With an initial guess of the
% coefficients and a contraction mapping we could attempt to find
% the "true" parameters 

% Transform grid from points on [-1,1] and evaluate function for these points
bma = .5*(b-a);
bpa = .5*(b+a);
f = ones(1,n); % placeholder for func evals
c = ones(1,n); % placeholder for coefficients

for k=1:n
    y=cos(pi*(k-.5)/n);
    f(k)=func(y*bma+bpa);
end

% Calculate coefficients

fac=2/n;

for j=1:n
    sum=0;
    for k=1:n
        sum=sum+f(k)*cos((pi*(j-1))*((k-.5)/n));
    end
    c(j)=fac*sum;
end

f=c;
    
