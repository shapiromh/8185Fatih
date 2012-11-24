function f= LINT(func,x,c)

% LINT Linear interpolation for scalar function on R1

% x is ordered grid, func is the function, c is the desired point

y = func(x);

% Test for extrapolation
if c>max(x)
    display 'This is outside of your range for interpolation! Be careful.';
    I = length(x);
    J = I-1;
    A = (x(I)-c)/(x(I)-x(J));
    B = 1-A;
    fc = A*y(J) + B*y(I);
end

if c<min(x)
    display 'This is outside of your range for interpolation! Be careful.';
    I = 2;
    J = I-1;
    A = (x(I)-c)/(x(I)-x(J));
    B = 1-A;
    fc = A*y(J) + B*y(I);
end

% Test for interpolating on grid point

[gridmin, I] = min(abs(c-x));

if gridmin==0
    f=func(x(I));
    return
end

% find segment of interpolation and create evaluate for interpolation

if x(1)<c && c<x(I)
    J = I-1;
    A = (x(I)-c)/(x(I)-x(J));
    B = 1-A;
    fc = A*y(J) + B*y(I);
elseif x(length(x))>c && c>x(I)
    J = I+1;
    A = (x(J)-c)/(x(J)-x(I));
    B = 1-A;
    fc = A*y(I) + B*y(J);
end

f=fc;

