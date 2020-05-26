function [v,dvdL] = PureFlexion(L,x);
% In adminesional coordinates:
% x = X/L reduced abscissa
% v = V/L reduced displacement along y
%
% Beam loaded by two momentums -M at X=0 and M at X=L
% Bending moment = M (for all x)
% Physical problem V = MX^2/2EI
% Adim problem v = (ML/2EI) x^2
% ==> Mechanical constant L = ML^2/2EI :

switch nargin
    
    case 0 % Initial value of the parameter
        v = 0;
                
    case 1 % Complementary points
        
        v = [];     % abscissa
        dvdL = [];  % name
        
    case 2 % Shape and derivatives

        % classical small perturbation approximation : parabolic shape 
        v  = L*x.^2;
        dvdL = x.^2;

        % exact large displacement circular solution is possible
end
