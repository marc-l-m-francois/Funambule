function [Xc1,Xc2,dXc1dL,dXc2dL] = Ellipse(L,x1,rien);
% x1 parametre abscisse curviligne entre 0 et 1
% Lk (lambda) parametres de la courbe
% L(1) L(2) centre
% L(3)/N angle (radians)
% L(4) rayon R1 
% L(5) rayon R2

N = 1; % Pour rendre la sensibilite de L(3) equivalente aux autres

switch nargin
    
    case 0 % Init manuelle cercle defini par 3 points ---------------------
        disp('    Click the great axis')
        [X1,X2] = ginput(2);plot(X1,X2,'r-');
        % Centre OC = OM1+M1M2/2
        C1 = X1(1)+(X1(2)-X1(1))/2;
        C2 = X2(1)+(X2(2)-X2(1))/2;
        % Rayon R1 = ||M1M2||/2
        R1 = sqrt( (X1(2)-X1(1))^2+(X2(2)-X2(1))^2 )/2;
        % Angle
        theta = -angle( (X1(2)-X1(1)) + i*(X2(2)-X2(1)) );
        % U Vecteur directeur de M1M2
        U1 = (X1(2)-X1(1))/(2*R1);
        U2 = (X2(2)-X2(1))/(2*R1);
        disp('    Click the small axis');
        [X1,X2] = ginput(2);plot(X1,X2,'r-');
        % Rayon R2 = || N1N2 v u ||/2
        R2 = abs( (X1(2)-X1(1))*U2 - (X2(2)-X2(1))*U1)/2;
        Xc1 = [C1,C2,theta*N,R1,R2];
 		
    case 1 % Traces supplementaires ---------------------------------------
        
        % Centres
% 		plot(L(1),L(2),'r+');
		% Axes
        x1=[0:1/1000:1];
   		X1 = L(4)*cos(2*pi*x1);
	    X2 = L(5)*sin(2*pi*x1);
		Xc1 = [  X1*cos(L(3)/N) + X2*sin(L(3)/N) + L(1)];
		Xc2 = [- X1*sin(L(3)/N) + X2*cos(L(3)/N) + L(2)];
		plot([L(1)-L(4)*cos(L(3)/N),L(1)+L(4)*cos(L(3)/N)],...
		     [L(2)+L(4)*sin(L(3)/N),L(2)-L(4)*sin(L(3)/N)],'r-');
		plot([L(1)+L(5)*sin(L(3)/N),L(1)-L(5)*sin(L(3)/N)],...
		     [L(2)+L(5)*cos(L(3)/N),L(2)-L(5)*cos(L(3)/N)],'r-');
%		text(L(1)-(L(4)*cos(L(3)/N)),L(2)+(L(4)*sin(L(3)/N)),...
%		'  \lambda_4','color','r');
%		text(L(1)+(L(5)*sin(L(3)/N)),L(2)+(L(5)*cos(L(3)/N)),...
%		'  \lambda_5','color','r');

    case 2 % Calcul simple ------------------------------------------------

         % En base locale
		X1 = L(4)*cos(2*pi*x1);
	    X2 = L(5)*sin(2*pi*x1);
		Xc1 = [  X1*cos(L(3)/N) + X2*sin(L(3)/N) + L(1)];
		Xc2 = [- X1*sin(L(3)/N) + X2*cos(L(3)/N) + L(2)];
        % Derivees
        dXc1dL = [  ones(size(x1)) ;
                    zeros(size(x1)) ;
                    - X1*sin(L(3)/N)/N + X2*cos(L(3)/N)/N;
                     cos(2*pi*x1)*cos(L(3)/N);
                     sin(2*pi*x1)*sin(L(3)/N)];
        dXc2dL = [  zeros(size(x1)) ;
                    ones(size(x1)) ;
                    - X1*cos(L(3)/N)/N - X2*sin(L(3)/N)/N;
                    -cos(2*pi*x1)*sin(L(3)/N);
                     sin(2*pi*x1)*cos(L(3)/N)];
                 
    case 3 % Informations- ------------------------------------------------
        
		disp(['    Ellipse of center [',num2str([L(1),L(2)]),']']);
        disp(['    angle (x1,great axis) ',num2str((L(3)/N)*180/pi),' degrees, ']);
        disp(['    radii R1 R2 [',num2str([L(4),L(5)]),']']);


end
