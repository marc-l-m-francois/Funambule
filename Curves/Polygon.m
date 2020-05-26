function [Xc1,Xc2,dXc1dL,dXc2dL] = Polygon(L,x1,rien);
% x1 parametre abscisse curviligne entre 0 et 1
% Lk (lambda) parametres de la courbe
% L(1) L(2) Centre
% L(3) rayon (centre-sommet)
% L(4) Angle

% Nombre de sommets
global number_of_summits 

switch nargin
    
    case 0 % Init manuelle cercle defini par 3 points ---------------------
        
        disp('    Please click the summits, enter afer the last one');
        figure(1);[x,y] = ginput;
        disp(['    Polygon with ',num2str(size(x,1)),' sides']);
        number_of_summits = size(x,1);
        plot(x,y,'r+');
        % Centre
        L(1) = mean(x);
        L(2) = mean(y);
        % Rayon
        L(3) = sqrt(mean( (x-L(1)).^2+(y-L(2)).^2 ));
        % Angle / abscisses
        L(4) = angle((x(1)-L(1)) + i*(y(1)-L(2))) + pi/size(x,1);
        Xc1 = L;
         
    case 1 % Traces supplementaires ---------------------------------------
        
        plot(L(1),L(2),'r+');
%        text(L(1),L(2),' C','FontSize',24,'Color','r');
        % Calcul des sommets en repere propre
		Sx = L(3)*cos(2*pi*[0:number_of_summits]/number_of_summits);
		Sy = L(3)*sin(2*pi*[0:number_of_summits]/number_of_summits);
   		% Mouvement de corps solide
		Xc1 = Sx*cos(L(4)+pi/number_of_summits) - Sy*sin(L(4)+pi/number_of_summits) + L(1);
		Xc2 = Sx*sin(L(4)+pi/number_of_summits) + Sy*cos(L(4)+pi/number_of_summits) + L(2);
        for p=1:number_of_summits
            text(Xc1(p),Xc2(p),[' P',num2str(p)],'FontSize',18,'Color','r')
        end

    case 2 % Calcul simple ------------------------------------------------

        % Calcul des sommets en repere propre
		Sx = L(3)*cos(2*pi*[0:number_of_summits]/number_of_summits);
		Sy = L(3)*sin(2*pi*[0:number_of_summits]/number_of_summits);
        
        x = zeros(size(x1));
        y = zeros(size(x1));
        for p=1:number_of_summits
            % Determination des indices dans un segment -x1(1) pour
            % recentrer en cas de derivee numerique...
            inds = x1-x1(1)>=(p-1)/number_of_summits & x1-x1(1)<p/number_of_summits;
            if p==number_of_summits;inds(end)=true;end
  			% Construction du poids relatif
			poids = number_of_summits*(x1-(p-1)/number_of_summits);
            noids = 1-poids;
            poids(~inds) = 0;
            noids(~inds) = 0;
            x = x+Sx(p+1)*poids + Sx(p)*noids;
            y = y+Sy(p+1)*poids + Sy(p)*noids;
         end
		% Mouvement de corps solide
		Xc1 = x*cos(L(4)+pi/number_of_summits) - y*sin(L(4)+pi/number_of_summits) + L(1);
		Xc2 = x*sin(L(4)+pi/number_of_summits) + y*cos(L(4)+pi/number_of_summits) + L(2);

        % Derivees       
        dXc1dL = [ ones(size(x1)) ;
                   zeros(size(x1)) ;
                   (Xc1-L(1))/L(3) ;
                   -x*sin(L(4)+pi/number_of_summits) - y*cos(L(4)+pi/number_of_summits)];
        dXc2dL = [ zeros(size(x1)) ;
                   ones(size(x1)) ;
                   (Xc2-L(2))/L(3) ;
                   x*cos(L(4)+pi/number_of_summits) - y*sin(L(4)+pi/number_of_summits)];
        
    case 3 % Informations- ------------------------------------------------

        disp(['    Regular polygon with ',num2str(number_of_summits),' sides']);
        
		disp(['    Center [',num2str([L(1),L(2)]),'], angle/x ',...
			       num2str((L(4))*180/pi),' degrees, ',...
			      ' radius [',num2str([L(3)]),']']);
       
        
end


% SOUCI Xc1(x1+1/1000) ne change pas le dernier terme !!!
% XcP1-Xc1