function f = SPLT(x,y,m,y2,c)

% SPLT Function for cubic spline interpolation
% x is grid, y is fn eval on grid, m is grid points
% y2 is 2nd deriv output of SPL, c is desired point of interpolation

% Error check

if length(x)~=m || length(y)~=m
    display 'You lied about the grid size.'
    return
end

% Find location on grid of c

klo =1;
khi =m;
while (khi-klo)>1
    k = (khi+klo)/2;
    if x(k)>c 
        khi=k;
    else
        klo=k;
    end
end

% Error check 2

h=x(khi)-x(klo);
if h == 0
    display 'Grid points must be distinct.'
    return
end

% Interpolation for point c using cubic spline

a=(x(khi)-c)/h;
b=(c-x(klo))/h;
yc=a*y(klo)+b*y(khi)+((a^3-a)*y2(klo)+(b^3-b)*y2(khi))*(h^2)/6;

f=yc;
