% Reference disc generator
% Used in 2019 publication


% CHOICE : PIXEL REFINEMENT


n = 100;


% INITS

% Initial point
X1 = 2;
X2 = 2;
% Length
L = 5;
% Angle
t = pi/10;
% Segment width
w = 1;
% Image size
F = ones(9);

% CALCULUS

% Fine image FF
xf1 = [1/(2*n):1/n:size(F,2)]+0.5;
xf2 = [1/(2*n):1/n:size(F,1)]+0.5;
[XF1,XF2] = meshgrid(xf1,xf2);
% Condition of appartenance
% A Faire