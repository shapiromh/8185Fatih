% Code: Written / Adopted (for Q2) by Matt Shapiro
% Date: November 16, 2012
% Macro 8185 - Prof. Fatih Guvenen
% Problem Set 2

clear,clc
pl = 0; % Set pl=1 to plot
%% Q1 Function Minimization

if pl==1
ezsurf('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',[-50:.1:50])
ezsurf('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',[-25:.1:25])
end

% options=optimset('Algorithm','interior-point');
x0 =[-25;25];
[q1, fvalq1] = fminunc('Q1',x0);

x1 =[15;15];
[q2, fvalq2] = fminunc('Q1',x1);

[q3, fvalq3] = fminsearch('Q1',x0);
[q4, fvalq4] = fminsearch('Q1',x1);


%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1d Fatih's Algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%

maxiter= 1000; % maximum iterations
xd = zeros(2,maxiter); % placeholder for starting points of N-M
z = zeros(2,maxiter); % placeholder for minimizers
Vz = ones(1,maxiter); % placeholder for minimum value
iter = 1; % Iteration Counter

theta = linspace(0,1,maxiter); % Increasing sequence of thetas
%theta = linspace(0,1,maxiter*2); % More gradual increasing sequence of thetas

xd(:,1)=[15;15];
[z(:,1) Vz(:,1)] = fminsearch('Q1',xd(:,1));
iter=iter+1;

A = sobolset(2); % Generate a 2D Sobol Sequence
B = net(A,maxiter); % Store the first n values of the transformed sequence generated
sobseq = (-50*B + (1-B)*50)'; % Transform elements of sequence into guess range (-50,50)

xd(:,iter)=theta(1)*z(:,1)+(1-theta(1))*sobseq(:,1);

while iter<maxiter
    [z(:,iter) Vz(:,iter)] = fminsearch2('Q1',xd(:,iter)); 
 
    % Update z* for algorithm and find new starting point for next iteration
    [temp1, temp2] = min(Vz); % find index for "best z" so far
    xd(:,iter+1)=theta(iter)*z(:,temp2)+(1-theta(iter))*sobseq(:,iter);
    iter = iter+1;
end

minVz = Vz(iter-1)

%% Q2a Utility Interpolation 

% Spline and Chebyshev code transformed for MATLAB from NR in F77

%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Interpolation %
%%%%%%%%%%%%%%%%%%%%%%%%

m = 20; % grid size
x = linspace(0.05,2,m); % create equal spaced grid 
c = .05:.05:2; % grid for interpolation
PW = [2 5 10];

% U = log(c)

linlog = zeros(1,length(c));
for j=1:length(c)
    linlog(j) = LINT(@log,x,c(j));
end

if pl==1
plot(c,linlog,c,log(c))
end
llogerror = abs(linlog-log(c));

% U = sqrt(c)

linsr = zeros(1,length(c));
for j=1:length(c)
    linsr(j) = LINT(@sqrt,x,c(j));
end

if pl==1
plot(c,linsr,c,log(c))
end
lsrerror = abs(linsr-log(c));

% U = c^(1-alf)/(1-alf)

linpwr=zeros(1,length(c),length(PW));
lpwrerror = zeros(1,length(c),length(PW));
for i=1:length(PW)
    pwr=@(c) (c.^(1-PW(i)))/(1-PW(i));
    for j=1:length(c)
        linpwr(1,j,i)=LINT(pwr,x,c(j));
    end
    lpwrerror(1,:,i)=abs(linpwr(1,:,i)-pwr(c));
end

if pl==1
    for i=1:length(PW)
        pwr=@(c) (c.^(1-PW(i)))/(1-PW(i));
        plot(c,linpwr(1,:,i),c,pwr)
    end    
end

%%%%%%%%%%%%%%%%%
% Cubic Splines %
%%%%%%%%%%%%%%%%%

n = 20; % grid size
x = linspace(0.05,2,n); % create equal spaced grid 
d = .05:.05:2; % grid for interpolation

% U = log(c)

ylog = log(x); % evaluate function on grid
yld1 = 1/x(1); % derivative on first grid point
yldn = 1/x(n); % derivative on last grid point
y2log = SPL(x,ylog,n,yld1,yldn); % Find second derivatives for interp fn

cublog = zeros(1,length(d));
for j=1:length(d)
    cublog(j) = SPLT(x,ylog,n,y2log,d(j));
end

if pl==1
plot(d,cublog,d,log(d))
end
clogerror = abs(cublog-log(d));

% U = sqrt(c)

ysr = sqrt(x); % evaluate function on grid
ysd1 = .5*x(1)^(-.5); % derivative on first grid point
ysdn = .5*x(n)^(-.5); % derivative on last grid point
y2sr = SPL(x,ysr,n,ysd1,ysdn); % Find second derivatives for interp fn

cubsr = zeros(1,length(d));
for j=1:length(d)
    cubsr(j) = SPLT(x,ysr,n,y2sr,d(j));
end

if pl==1
    plot(d,cubsr,d,sqrt(d))
end
csrerror = abs(cubsr-sqrt(d));


% U = c^(1-alf)/(1-alf)

% Placeholders for evaluations with each value of alpha 
ypr=zeros(1,length(x),length(PW));
ypd1=zeros(1,length(PW));
ypdn=zeros(1,length(PW));
y2pr=zeros(1,length(x),length(PW));
cubpwr = zeros(1,length(d),length(PW));
cpwrerror = zeros(1,length(d),length(PW));

for i=1:length(PW)
    ypr(1,:,i) = (x.^(1-PW(i)))/(1-PW(i)); % evaluate function on grid
    ypd1(1,i) = x(1)^(-PW(i)); % derivative on first grid point
    ypdn(1,i) = x(n)^(-PW(i)); % derivative on last grid point
    y2pr(1,:,i) = SPL(x,ypr(1,:,i),n,ypd1(1,i),ypdn(1,i)); % Find second derivatives for interp fn
    for j=1:length(d)
        cubpwr(1,j,i)=SPLT(x,ypr(1,:,i),n,y2pr(1,:,i),d(j));
    end
    if pl==1
        pwr=@(d) (d.^(1-PW(i)))/(1-PW(i));
        plot (d,cubpwr(1,:,i),d,pwr)
    end
    cpwrerror(1,:,i)=abs(cubpwr(1,:,i)-(d.^(1-PW(i)))/(1-PW(i)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Chebyshev Polynomials %
%%%%%%%%%%%%%%%%%%%%%%%%% 

max = 2; % max on grid
min = .05; % min on grid
h = 10; % number of grid points
o = 5; % order of Cheby polynomial, must be lower than h
e = .05:.05:2; % grid for interpolation

% U = log(c)

logcof = CHEBYC(min,max,h,@log); % Calculate coefficients for polynomial
chelog = zeros(1,length(e));
for j=1:length(e)
  chelog(j)= CHEBY(min,max,logcof,o,e(j)); % Evaluate for interpolation on 
end

if pl==1
plot(e,chelog,e,log(e))
end
chlogerror = abs(chelog-log(e));

% U = sqrt(c)

sqrcof = CHEBYC(min,max,h,@sqrt);
chesqr = zeros(1,length(e));
for j=1:n
    chesqr(j)=CHEBY(min,max,sqrcof,o,e(j));
end

if pl==1
plot(e,chelog,e,log(e))
end
chsrerror = abs(chelog-log(e));

% U = c^(1-alf)/(1-alf)

pwrcof = zeros(1,h,length(PW));
chepwr = zeros(1,length(e),length(PW));
chpwrerror = zeros(1,length(e),length(PW));

for i=1:length(PW)
    pwr=@(c) (c.^(1-PW(i)))/(1-PW(i));
    pwrcof(1,:,i)=CHEBYC(min,max,h,pwr);
    for j=1:n
        chepwr(1,j,i)=CHEBY(min,max,pwrcof(1,:,i),o,e(j));
    end
    chpwrerror(1,:,i)=abs(chepwr(1,:,i)-pwr(e));
end


if pl==1
    for i=1:length(PW)
        pwr=@(e) (e.^(1-PW(i)))/(1-PW(i));
        if pl==1
        plot(e,chepwr(1,:,i),e,pwr)
        end
    end    
end


%% Q2c Utility Extrapolation
