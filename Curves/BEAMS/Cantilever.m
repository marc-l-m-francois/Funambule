function [v,dvdL] = Cantilever(L,x);
% In adminesional coordinates:
% x = X/L reduced abscissa
% v = V/L reduced displacement along y
%
% Beam loaded by a concentrated force F at its end X=L
% Bending moment Mf = F(L-X)
% Physical problem V = FX^2(3L-X)/6EI
% Adim problem v = (FL^2/6EI) x^2(3-x)
% ==> Mechanical constant L = FL^2/6EI :

switch nargin
    
    case 0 % Initial value of the parameter
        v = 0;
        
    case 1 % Complementary points
        
        v = [];     % abscissa
        dvdL = [];  % name
        
    case 2 % Shape and derivatives

        v = L*x.^2.*(3-x);
        dvdL = x.^2.*(3-x);

end