function [v,dvdL] = Cubic(L,x)
% In adminesional coordinates:
% x = X/L reduced abscissa
% v = V/L reduced displacement along y
%
% General cubic shape
% First and zero orders are included inthe the motion of the reference
% straight shape
% order zero  <=> location delta
% order one <=> angle theta
% orders two and three :

switch nargin
    
    case 0 % Initial value of the parameters
        v = [0 0];
        
    case 1 % Complementary points
        
        v = [];     % abscissa
        dvdL = [];  % name
                
    case 2 % Shape and derivatives

        v = L(2)*x.^3+L(1)*x.^2;
        dvdL = [x.^2;x.^3];
        
end
