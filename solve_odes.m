% X ->[k] 2X
% dX/dt = kX
clear; clc;
syms X(t) X0 k;
conds = X(0) == X0;
odes = diff(X, t) == k*X;
dsolve(odes, conds)

% X + Y ->[k] 2Y
% dY/dt = -kXY + 2Y
clear;
syms Y(t) X(t) k C;
% conds = Y(0) == Y0;
odes = diff(X, t) == -k*X*(C-X);
dsolve(odes)

% X + Y ->[k] 2Y
% dY/dt = kX
clear;
syms Y(t) X(t) X0 Y0;
conds = [X(0) == X0; Y(0) == Y0];
odes = [diff(X, t) == -X*Y;  diff(Y, t) == -X*Y + 2*Y];
dsolve(odes, conds)