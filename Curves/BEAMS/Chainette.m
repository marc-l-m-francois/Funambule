function [v,dvdL] = Chainette(L,x)
% In adminesional coordinates:
% x = X/L reduced abscissa
% v = V/L reduced displacement along y
%
% The chainette curve is the shape of a string (beam without bending
% stiffness) under a horizontal constant tension T and a vertical
% distributed load p
% 
% Physical problem V = (T/p)(ch(pX/T)-1)
% Adim problem v = (N/pL)(ch(pLx/N)-1 ) plus shift : x ==> (x-x0)
% ==> Mechanical constant L(1) = pL/N 
% ==> and shift           L(2) = x0 :


global eps_eig
eps_eig = 1E-3;

switch nargin
    
    case 0 % Initial parameters
        v = [0 0.5];
        
    case 1 % Complementary points
        
        v = L(2);               % abscissa
        dvdL = num2str(L(2));   % name
        
    case 2 % Shape and derivatives
        
        K = L(1);  
        xo= L(2);

        if K~=0
            v = (cosh(K*(x-xo))-1)/K;
            dvdL = [ ( (x-xo).*sinh(K*(x-xo))/K - (cosh(K*(x-xo))-1)/K^2 );
                    -sinh(K*(x-xo))];
        else
           v = zeros(1,size(x,2));
           dvdL = [ (x-xo).^2;
                   -sinh(K*(x-xo))];
        end
       
end