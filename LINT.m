function f= LINT(func,x,c)

% LINT Linear interpolation for scalar function on R1

% x is ordered grid, func is the function, c is the desired point

% Error test for extrapolation
if c>max(x) || c<min(x)
    display 'Do not extrapolate!';
    return
end

% Error test for interpolating on grid point

[gridmin, I] = min(abs(c-x));

if gridmin==0
    display 'You know this!';
    f=func(x(I));
    return
end

y = func(x);

% find segment of interpolation and create evaluate for interpolation

if c<x(I)
    J = I-1;
    A = (x(I)-c)/(x(I)-x(J));
    B = 1-A;
    fc = A*y(J) + B*y(I);
elseif c>x(I)
    J = I+1;
    A = (x(J)-c)/(x(J)-x(I));
    B = 1-A;
    fc = A*y(I) + B*y(J);
end

f=fc;

