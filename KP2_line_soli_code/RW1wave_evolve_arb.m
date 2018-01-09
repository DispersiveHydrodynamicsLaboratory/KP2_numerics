function [R1] = RW1wave_evolve_arb(R1m,R1p,R2,m,choice)
% Calculates evolution of a jump in R1
% where R2 remains constant
% and the jump occurs over a line with
% slope m
% Choice determines which solution branch is taken

% Calculate slope of line perpendicular to m
mrec = -1/m;

% Calculate characteristic speeds
W1 = @(R1) -m/(2*sqrt(1+m^2)) * (R1 + 2/m) * (R1 - 2*R2);
W1m = W1(R1m);
W1p = W1(R1p);


