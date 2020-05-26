function [Xc1,Xc2,dXc1dL,dXc2dL] = Cercle(L,x1,rien);
% x1 parametre abscisse curviligne entre 0 et 1
% Lk (lambda) parametres de la courbe
% L(1) L(2) centre
% L(3) rayon


switch nargin
    
    case 0 % Init manuelle cercle defini par 3 points ---------------------
        
        disp('   Please click 3 points');
        figure(1);[x,y] = ginput(3);
        % Centre / x
        L(1) = -( ( (x(3)^2-x(2)^2+y(3)^2-y(2)^2)/(2*(y(3)-y(2))) ) - ...
        ( (x(2)^2-x(1)^2+y(2)^2-y(1)^2)/(2*(y(2)-y(1))) ) ) / ...
        ( (x(2)-x(1))/(y(2)-y(1))-(x(3)-x(2))/(y(3)-y(2)) );
        % Centre / y
        L(2) = -L(1)*(x(2)-x(1))/(y(2)-y(1)) + (x(2)^2-x(1)^2+y(2)^2-y(1)^2)/(2*(y(2)-y(1)));
        % Rayon
        L(3) = sqrt((x(1)-L(1))^2+(y(1)-L(2))^2);
        Xc1 = L; % Xc2 = [];dXc1dL = [];dXc2dL = [];
        
    case 1 % Traces supplementaires ---------------------------------------
        
        plot(L(1),L(2),'r+');

    case 2 % Calcul simple ------------------------------------------------

        Xc1 = L(1) + L(3)*cos(2*pi*x1);
        Xc2 = L(2) + L(3)*sin(2*pi*x1);
        dXc1dL = [ ones(size(x1)) ; zeros(size(x1)) ; cos(2*pi*x1)];
        dXc2dL = [zeros(size(x1)) ;  ones(size(x1)) ; sin(2*pi*x1)];
        
    case 3 % Informations -------------------------------------------------
        
		disp(['    Circle of center [',num2str([L(1),L(2)]),'], and radius [',num2str([L(3)]),']']);


end