% Code: Written by Matt Shapiro
% Date: September 23rd, 2012
% Macro 8185 - Prof. Fatih Guvenen
% Problem Set 1

clear, clc

%% Q1 First Method

% Quadratic Formula

quad1 = @(a,b,c) (-b+sqrt(b^2-4*a*c))/2*a;
quad2 = @(a,b,c) (-b-sqrt(b^2-4*a*c))/2*a;

% Parameters
n= -8:1:-1;

a = 1;
b = 100000;
c = ones(size(n));
for j=1:length(n)
    c(j)=10^(n(j));
end
quadm1 = ones(size(n));

% Formulas

quadm1(:) = max(quad1(a,b,c(:)),quad2(a,b,c(:)));

%% Q1 Second Method

quadm2 = ones(size(n));

q = @(a,b,c) (1/2)*(b+sign(b)*sqrt(b^2-4*a*c));
q1 = @(q,a) q/a;
q2 = @(c,q) c./q;

quadm2(:)= max(q1(q(a,b,c(:)),a),q2(c(:),q(a,b,c(:))));

quaddiff = quadm1 - quadm2; % difference in evaluations of two methods

%% Q2 Truncation Errors

m=2:1:20;
phim1=ones(1,21);
phim2=ones(1,21);

% Recursion Method

phim1(1)=1;
phim2(2)=.61803398;

for j=3:21
    phim1(j)=phim1(j-2)-phim1(j-1);
end

% Straight Exponentiation

phim2(1)=1;
phim2(2)=.61803398;

for j=3:21
    phim2(j)=phim2(2)^(j-1);
end

phimdiff = phim1 -phim2;

%% Q3 One sided derivative

g = @(p) .5*p^(-.5) + .5*p^(-.2);
k = 1:1:10;
gdest = ones(1,length(k));

for j=1:length(k)
    gdest(j)=(g(1.5+10^(-k(j)))-g(1.5))/(10^(-k(j)));
end

% True Derivative

gprime = @(p) -.25*p^(-1.5) - .10*p^(-1.2); 

gdiff = gprime(1.5)-gdest;
[C I] = min(abs(gdiff)); % I returns (-) exponent associated with epsilon returning most accurate result

%% Q4 Two Sided Derivative

gdest2=ones(1,length(k));

for j=1:length(k)
    gdest2(j)=(g(1.5+10^(-k(j)))-g(1.5-+10^(-k(j))))/(10^(-k(j)));
end

gdiff2 = gprime(1.5)-gdest2;
[D J] = min(abs(gdiff2));

%% Q7 Neoclassical Growth Model, Code originally for 8105

% Parameters and Grid

A = 2;
beta = .8;
alpha = .7; 
delta = .8; 
kss = ((1/(beta*alpha*A)-(1-delta)/(alpha*A)))^(1/(alpha-1)); % steady state capital
kmin = .1*kss; % minimum capital level on grid
kmax = 1.5*kss; % largest capital level on grid
kgrid = [kmin:.01:kmax]'; % grid, adjust increment
gridsize = size(kgrid,1);
iter = 100; % iterations allowed by the program
crit = 10^(-6); % critical value for convergence

optcap = zeros(length(kgrid),iter); % space for optimal capital policy for every cap level at each iter
optval = zeros(length(kgrid),iter); % same for value function


% Production possibilities and consumption possibilities

PPF = kron(ones(1,gridsize),(1-delta)*kgrid+A*kgrid.^alpha); % PPF using kgrid element capital
k1 = kgrid'; % choices for next state capital
cons = PPF - kron(ones(gridsize,1),k1); %(a,b) in this matrix is consumption with a capital level and b capital choice
util = log(cons); % utility for all capital and capital choice pairs
V = util+beta*kron(ones(gridsize,1),optval(:,1)'); % first iteration of the value fn

for j=2:iter
    [optval(:,j),kpolind] = max(V,[],2); % kpolind is index of policy function (found along columns)
    optcap(:,j) = kgrid(kpolind); % kpol is policy function
    error = norm(optval(:,j)-optval(:,j-1));
    if error < crit
        conind = j; % iteration at which we converged
        break
    end
    V= util + beta*kron(ones(gridsize,1),optval(:,j)'); % updating value function
end

optcap(:,conind)
optval(:,conind)