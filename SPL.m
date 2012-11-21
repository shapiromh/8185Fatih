function f = SPL(x,y,m,yp1,ypm)

% SPL Generate second derivatives of cubic spline for the function
% x is the grid, y is function evaluated on grid, m is number of grid points
% yp1 and ypn are the first derivatives of the orig fn on the first and last nodes

% Error check

if length(x)~=m || length(y)~=m
    display 'You lied about the grid size.'
    return
end

% Lower boundary conditions on second derivative

y2 = zeros(1,m); % holder for second derivatives
u = zeros(1,m-1); % holder for something used in tridiagonal matrix solver

if yp1 > 1*10^30
    y2(1)=0;
    u(1)=0;
else
    y2(1)=-0.5;
    u(1)=(3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1);
end

% Part of algorithm for solving tridiagonal matrix problem

for j=2:m-1
    sig=(x(j)-x(j-1))/(x(j+1)-x(j-1));
    p=sig*y2(j-1)+2;
    y2(j)=(sig-1)/p;
    u(j)=(6*((y(j+1)-y(j))/(x(j+1)-x(j))-(y(j)-y(j-1))/(x(j)-x(j-1)))/ ...
        (x(j+1)-x(j-1))-sig*u(j-1))/p;
end

% Finish algorithm with upper boundary condition and backward substitution 

if ypm > 1*10^30 
    qm =0;
    um =0;
else
    qm =.5;
    um=(3/(x(m)-x(m-1)))*(ypm-(y(m)-y(m-1))/(x(m)-x(m-1)));
end

y2(m)=(um-qm*u(m-1))/(qm-y2(m-1)+1);

for i=m-1:-1:1
    y2(i)=y2(i)*y2(i+1)+u(i);
end

f=y2;