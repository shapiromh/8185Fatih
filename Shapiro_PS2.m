% Code: Written by Matt Shapiro
% Date: November 16, 2012
% Macro 8185 - Prof. Fatih Guvenen
% Problem Set 2

%% Q1 Function Minimization

ezsurf('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',[-50:.1:50])

x0 =[25;25];
[q1, fvalq1] = fmincon('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',x0);

x1 =[15;15];
[q2, fvalq2] = fmincon('(5*sin(x)/x)*(max(20-abs(y),0))^1.2',x0);

fminsearch

%% Q2 Utility Interpolation 