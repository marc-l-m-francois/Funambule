function [v,dvdL] = Segment(L,x)
% In adminesional coordinates:
% x = X/L reduced abscissa
% v = V/L reduced displacement along y
%
% Straight segment - null function

switch nargin
    
    case 0 % Initial value of the parameters
        v = [];
        
    case 1 % Complementary points
        
        v = [];     % abscissa
        dvdL = [];  % name
                
    case 2 % Shape and derivatives

        v = zeros(size(x));
        dvdL = zeros(size(x));
        
end
