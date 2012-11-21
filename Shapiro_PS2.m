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



% Q1d Fatih's Algorithm %

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

%% Q2 Utility Interpolation 

% Spline code transformed for MATLAB from NR in F77

% Linear Interpolation

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

if plot==1
    for i=1:length(PW)
           
    
end

% Cubic Splines

% U = log(c)

% U = sqrt(c)

% U = c^(1-alf)/(1-alf)

% Chebyshev Polynomials


