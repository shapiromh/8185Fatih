% Code: Written / Adopted (for Q2) by Matt Shapiro
% Date: November 16, 2012
% Macro 8185 - Prof. Fatih Guvenen
% Problem Set 2

clear,clc
format long
pl = 0; % Set pl=1 to plot

%% Q1 Function Minimization

if pl==1
figure(1)
subplot(1,2,1)
ezsurf('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',[-50:.1:50])
subplot(1,2,2)
ezsurf('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',[-50:.1:50])

figure(2)
subplot(1,2,1)
ezsurf('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',[-25:.1:25])
subplot(1,2,2)
ezsurf('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',[-25:.1:25])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton and Nelder-Mead Minimization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quasi-Newton Search

x0 =[-25;25];
[q1, fvalq1] = fminunc('Q1',x0);
% Minimum found at (-25,25) for value 0

x1 =[15;15];
[q2, fvalq2] = fminunc('Q1',x1);
% Minimum found at (17.583054987010893,0) for value -9.878399392795970 


% Nelder-Mead Search

[q3, fvalq3] = fminsearch('Q1',x0);
% Minimum found at (-25,25) for value 0
[q4, fvalq4] =fminsearch2('Q1',x0); % Simplex length 10 pts from start
% Minimum found at (-29.811596, 5.055848e(-10)) for value -6.103466

[q5, fvalq5] = fminsearch('Q1',x1);
% Minimum found at (17.220505542191308, -1.58e(-6)) for value -10.554137
[q6, fvalq6] = fminsearch2('Q1',x1); % Simplex length 10 pts from start
% Minimum found at same point with same value



%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1d Fatih's Algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%

maxiter= 10000; % iterations
xd = zeros(2,maxiter); % placeholder for starting points of N-M
z = zeros(2,maxiter); % placeholder for minimizers
VZ = ones(1,maxiter); % placeholder for minimum value
iter = 1; % Iteration Counter

theta = linspace(0,1,maxiter*1.5); % Increasing sequence of thetas
%theta = linspace(0,1,maxiter*2); % More gradual increasing sequence of thetas

xd(:,1)=[15;15];
[z(:,1) VZ(:,1)] = fminsearch('Q1',xd(:,1));
iter=iter+1;

A = sobolset(2); % Generate a 2D Sobol Sequence
B = net(A,maxiter); % Store the first n values of the transformed sequence generated
sobseq = (-50*B + (1-B)*50)'; % Transform elements of sequence into guess range (-50,50)

xd(:,iter)=theta(1)*z(:,1)+(1-theta(1))*sobseq(:,1);

while iter<maxiter
    [z(:,iter) VZ(:,iter)] = fminsearch2('Q1',xd(:,iter)); 
 
    % Update z* for algorithm and find new starting point for next iteration
    [temp1, temp2] = min(VZ); % find index for "best z" so far
    xd(:,iter+1)=theta(iter)*z(:,temp2)+(1-theta(iter))*sobseq(:,iter);
    iter = iter+1;
end

[minVZ, minzi] = min(VZ); % minimum value over all starting points
minz=z(:,minzi); % minimizer
% minimum value of -39.548776714419724 on iteration 363 at point (-4.4934, 2.8e(-9))
% for faster minimum


%% Q2a Utility Interpolation 

% Spline and Chebyshev code transformed for MATLAB from NR in F77

%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Interpolation %
%%%%%%%%%%%%%%%%%%%%%%%%

m = 4; % initial grid size
x = linspace(0.05,2,m); % create initial equal spaced grid 
c = .05:.05:2; % grid for interpolation
PW = [2 5 10];

% U = log(c)

linlog = zeros(1,length(c));

while max(abs(linlog-log(c))>abs(.01*log(c)))>0 % increase grid size
    % until 
    m=m+1; % increase the size of the grid to improve the fit
    x = linspace(0.05,2,m); 
    for j=1:length(c)
        linlog(j) = LINT(@log,x,c(j));
    end
end

linloggs = m; % grid size required to get appropriate fit for log(c)
% 40 points necessary (i.e. through every point)

if pl==1
figure(3)
plot(c,linlog,c,log(c))
end

% U = sqrt(c)

m = 4; % initial grid size
x = linspace(0.05,2,m); % create initial equal spaced grid 

linsr = zeros(1,length(c));
%

while max(abs(linsr-sqrt(c))>abs(.01*sqrt(c)))>0
    m=m+1;
    x=linspace(0.05,2,m);
    for j=1:length(c)
        linsr(j) = LINT(@sqrt,x,c(j));
    end
end

linsrgs = m; % grid size require to get appropriate fit for sqrt(c)
% 33 points necessary

if pl==1
figure(4)
plot(c,linsr,c,sqrt(c))
end

% U = c^(1-alf)/(1-alf)


linpwr=zeros(1,length(c),length(PW));
linpwrgs=zeros(1,length(PW)); % holder for grid sizes for each power

for i=1:length(PW)
    m = 4; % initial grid size
    x = linspace(0.05,2,m); % create initial equal spaced grid 
    pwr=@(c) (c.^(1-PW(i)))/(1-PW(i));
    while max(abs(linpwr(1,:,i)-pwr(c))>abs(.01*pwr(c)))>0
        m=m+1;
        x=linspace(0.05,2,m);
        for j=1:length(c)
            linpwr(1,j,i)=LINT(pwr,x,c(j));
        end
    end
    linpwrgs(1,i)=m;
    % for all powers 40 points necessary
    % clearly linear interpolation is no good for getting such fine estimates
end

if pl==1
    figure(5)
    for i=1:length(PW)
        pwr=@(c) (c.^(1-PW(i)))/(1-PW(i));
        subplot(2,2,i)
        plot(c,linpwr(1,:,i),c,pwr)
    end    
end

%%%%%%%%%%%%%%%%%
% Cubic Splines %
%%%%%%%%%%%%%%%%%

n = 4; % initial grid size
x = linspace(0.05,2,n); % create equal spaced grid 
d = .05:.05:2; % grid for interpolation

% U = log(c)

cublog = zeros(1,length(d));

while max(abs(cublog-log(d))>abs(.01*log(d)))>0
    n=n+1;
    x=linspace(.05,2,n);
    ylog = log(x); % evaluate function on grid
    yld1 = 1/x(1); % derivative on first grid point
    yldn = 1/x(n); % derivative on last grid point
    y2log = SPL(x,ylog,n,yld1,yldn); % Find second derivatives for interp fn
    for j=1:length(d)
        cublog(j) = SPLT(x,ylog,n,y2log,d(j));
    end
end

cubloggs=n; % grid points required to get appropriate fit for function
% 40 points necessary, again spline goes through every point

if pl==1
figure(6)
plot(d,cublog,d,log(d))
end

% U = sqrt(c)

n = 4; % reset initial grid size
x = linspace(0.05,2,n); % create equal spaced grid 
cubsr = zeros(1,length(d));

while max(abs(cubsr-sqrt(d))>abs(.01*sqrt(d)))>0
    n=n+1;
    x=linspace(.05,2,n);
    ysr = sqrt(x); % evaluate function on grid
    ysd1 = .5*x(1)^(-.5); % derivative on first grid point
    ysdn = .5*x(n)^(-.5); % derivative on last grid point
    y2sr = SPL(x,ysr,n,ysd1,ysdn); % Find second derivatives for interp fn
    for j=1:length(d)
        cubsr(j) = SPLT(x,ysr,n,y2sr,d(j));
    end
end

cubsrgs=n; % grid points required to get appropriate fit for function
% 15 points necessary

if pl==1
figure(7)
    plot(d,cubsr,d,sqrt(d))
end


% U = c^(1-alf)/(1-alf)

% Initial placeholders for evaluations with each value of alpha 
ypr=zeros(1,length(x),length(PW));
ypd1=zeros(1,length(PW));
ypdn=zeros(1,length(PW));
y2pr=zeros(1,length(x),length(PW));
cubpwr = zeros(1,length(d),length(PW));
cubpwrgs = zeros(1,length(PW));

n = 4; % reset initial grid size
x = linspace(0.05,2,n); % create equal spaced grid 

for i=1:length(PW)
    pwr=@(c) (c.^(1-PW(i)))/(1-PW(i));
    while max(abs(cubpwr(1,:,i)-pwr(d))>abs(.01*pwr(d)))>0
        n = n+1;
        x = linspace(0.05,2,n);
        ypr=zeros(1,length(x),length(PW));
        y2pr=zeros(1,length(x),length(PW));
        ypr(1,:,i) = pwr(x); % evaluate function on grid
        ypd1(1,i) = x(1)^(-PW(i)); % derivative on first grid point
        ypdn(1,i) = x(n)^(-PW(i)); % derivative on last grid point
        y2pr(1,:,i) = SPL(x,ypr(1,:,i),n,ypd1(1,i),ypdn(1,i)); % Find second derivatives for interp fn
        for j=1:length(d)
            cubpwr(1,j,i)=SPLT(x,ypr(1,:,i),n,y2pr(1,:,i),d(j));
        end
    end
    cubpwrgs(1,i)=n; % grid points required to get appropriate fit for function
    % 36, 40, 79 points necessary
    % 79 points necessary is curious since we are only estimating 40 points
end

if pl==1
    figure(8)
    for i=1:length(PW)
        pwr=@(d) (d.^(1-PW(i)))/(1-PW(i));
        subplot(2,2,i)
        plot (d,cubpwr(1,:,i),d,pwr)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% Chebyshev Polynomials %
%%%%%%%%%%%%%%%%%%%%%%%%% 

gmax = 2; % max on grid
gmin = .05; % min on grid
h = 50000; % max order of polynomial allowed
e = .05:.05:2; % grid for interpolation

% U = log(c)

logcof = CHEBYC(gmin,gmax,h,@log); % Calculate coefficients for polynomial
chelog = zeros(1,length(e));

o = 3; % initial order of Cheby polynomial, must be lower than h
while max(abs(chelog-log(e))>abs(.01*log(e)))>0
    o=o+1;
    for j=1:length(e)
      chelog(j)= CHEBY(gmin,gmax,logcof,o,e(j)); % Evaluate for interpolation on 
    end
end

cheloggs = o; % order of polynomials needed to get appropriate approximation
% 3622 order polynomial needed
    
if pl==1
figure(9)
plot(e,chelog,e,log(e))
end

% U = sqrt(c)

sqrcof = CHEBYC(gmin,gmax,h,@sqrt);
chesqr = zeros(1,length(e));

o = 1; % initial guess 
while max(abs(chesqr - sqrt(e))>abs(.01*sqrt(e)))>0
    o=o+1;
    for j=1:length(e)
        chesqr(j)=CHEBY(gmin,gmax,sqrcof,o,e(j));
    end
end

chesqrgs = o; % order of polynomials needed to get appropriate approximation
% 9th order polynomial needed

if pl==1
figure(10)
plot(e,chesqr,e,sqrt(e))
end

% U = c^(1-alf)/(1-alf)

pwrcof = zeros(1,h,length(PW));
chepwr = zeros(1,length(e),length(PW));
chepwrgs = zeros(1,length(PW));

for i=1:3
    o = 5; % initial guess
    pwr=@(c) (c.^(1-PW(i)))/(1-PW(i));
    pwrcof(1,:,i)=CHEBYC(gmin,gmax,h,pwr);
    while max(abs(chepwr(1,:,i)-pwr(e))>abs(.01*pwr(e)))>0
        o=o+1;
        for j=1:length(e)
            chepwr(1,j,i)=CHEBY(gmin,gmax,pwrcof(1,:,i),o,e(j));
        end
    end
    chepwrgs(i) = o; % order of polynomials needed to get appropriate approximation
    % 21 order polynomial needed for power 2, 74 for power 5
    % The order required for power 10 was over 50000
end


if pl==1
    figure(11)
    for i=1:length(PW)
        pwr=@(e) (e.^(1-PW(i)))/(1-PW(i));
        subplot(2,2,i)
        plot(e,chepwr(1,:,i),e,pwr(e))
    end    
end


%% Q2c Utility Extrapolation

% Extrapolation not feasible with Chebyshev polynomials
% The approximation is not defined outside of the original range set

% Extrapolation Grid

extrap = [.02 2.1 2.5 4]; 

%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Interpolation %
%%%%%%%%%%%%%%%%%%%%%%%%

% U = log(c)

llogextra=zeros(1,length(extrap)); % holder for extrapolation values 
llogextraerr=zeros(1,length(extrap)); % holder for errors from extrapolation
m = linloggs; % take the good linear model from part a
x = linspace(.05,2,m); 

for i=1:length(extrap)
    llogextra(i) = LINT(@log,x,extrap(i));
end

llogextraerr=abs(llogextra-log(extrap));
% errors .500402423538188 .001845451799148 .030034528528689 .319565138811648

% U = sqrt(c)

lsqrextra=zeros(1,length(extrap)); % holder for extrapolation values 
lsqrextraerr=zeros(1,length(extrap)); % holder for errors from extrapolation
m = linsrgs; % take the good linear model from part a
x = linspace(.05,2,m); 

for i=1:length(extrap)
    lsqrextra(i)=LINT(@sqrt,x,extrap(i));
end

lsqrextraerr=abs(lsqrextra-sqrt(extrap));
% errors .028294478242207 .000704717636790 .011218881702099 .126790160025870

% U = c^(1-alf)/(1-alf)

lpwrextra=zeros(1,length(extrap),length(PW)); % holder for extrapolation values 
lpwrextraerr=zeros(1,length(extrap),length(PW)); % holder for errors from extrapolation

for j=1:length(PW) 
    m = linpwrgs(1,j); % take the good linear model from part a
    x = linspace(.05,2,m);
    pwr=@(c) (c.^(1-PW(j)))/(1-PW(j));
    for i=1:length(extrap)
        lpwrextra(1,i,j)=LINT(pwr,x,extrap(i));
    end
    lpwrextraerr(1,:,j)=abs(lpwrextra(1,:,j)-pwr(extrap));
    % errors(2) 24 .001831601831502 .028205128205129 .262820512820515
    % errors(5) 1.5x10^(6) .000560249 .007427616 .051962025 
    % errors(10) 2.169x10^(14) .0000394717245425185 .00036747309 .0020048494
end


%%%%%%%%%%%%%%%%%
% Cubic Splines %
%%%%%%%%%%%%%%%%%

% U = log(c)

clogextra=zeros(1,length(extrap)); % holder for extrapolation values 
clogextraerr=zeros(1,length(extrap)); % holder for errors from extrapolation
m = cubloggs; % take the good linear model from part a
x = linspace(.02,2,m); 

for i=1:length(extrap)
    clogextra(i) = SPLT(x,ylog,m,y2log,extrap(i));
end

clogextraerr=abs(clogextra-log(extrap));
% errors .916290731874155 .000580904959675 .080266719211391 4.6366846615

% U = sqrt(c)

csrextra=zeros(1,length(extrap)); % holder for extrapolation values 
csrextraerr=zeros(1,length(extrap)); % holder for errors from extrapolation
m = cubsrgs; % take the good linear model from part a
x = linspace(.02,2,m); 

for i=1:length(extrap)
    csrextra(i) = SPLT(x,ysr,m,y2sr,extrap(i));
end

csrextraerr=abs(csrextra-sqrt(extrap));
% errors .082185441512669 .000077132454167 .012484446698173 .6606101902

% U = c^(1-alf)/(1-alf)

cpwrextra=zeros(1,length(extrap),length(PW)); % holder for extrapolation values 
cpwrextraerr=zeros(1,length(extrap),length(PW)); % holder for errors from extrapolation

for j=1:length(PW) 
    m = cubpwrgs(1,j);
    x = linspace(.05,2,m);
    yprex = pwr(x); % evaluate function on grid
    ypd1ex= x(1)^(-PW(j)); % derivative on first grid point
    ypdnex= x(m)^(-PW(j)); % derivative on last grid point
    y2prex= SPL(x,yprex,m,ypd1ex,ypdnex); % Find second derivatives for interp fn
    pwr=@(c) (c.^(1-PW(j)))/(1-PW(j));
    for i=1:length(extrap)
        cpwrextra(1,i,j)=SPLT(x,yprex,m,y2prex,extrap(i));
    end
    cpwrextraerr(1,:,j)=abs(cpwrextra(1,:,j)-pwr(extrap));
    % errors(2) 12.181147737 .00133306933 .079410917 4.4049939
    % errors(5) 1.3414x10^(6) .541180419 5.211100186 261.482167354
    % errors(10) 2.157x10^(14) .0412714 2.3928549 1.4010x10^(2)
end
