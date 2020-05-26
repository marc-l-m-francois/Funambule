function [Xc1,Xc2,dXc1dL,dXc2dL] = IBOX_cote(L,x1,rien);
% x1 parametre abscisse curviligne entre 0 et 1
% Lk (lambda) parametres de la courbe
% L(1) L(2) centre sphere
% L(3) angle / abscisses
% L(4) echelle = rayon sphere

global Rsph Rcyl Rpet LAo 
global theta_m theta_n theta_p Lmn Lno Lop Ltot houb

switch nargin
    
    case 0 % Init manuelle cercle defini par 3 points ---------------------
        
        % Dimensions principales mesurees sur le plan
        Rsph        = 98.88; % mm
        Rcyl        = 68.6;
        Rpet		= 30;
        LAP 		= 57;
        LAm 		= 51+6;
        LAo 		= 241.09;
        LAp 		= 265.62;
        entraxe = 350;
        % Zone des gravures "FMGC" situees 1 m plus bas
        x_d_grav = [122 150 160 178 190 210 228 263]*241.09/318;
        x_f_grav = [148 158 176 185 200 222 235 301]*241.09/318;;
        % Epaisseur gravure "FMGC"
        ep_grav = 1;
        % Zone des gravures de code coulee sur la sphere
        x_d_gras = -15;
        x_f_gras =  14;

		disp('   Click 3 points of the head circle and 1 point on the tail axis');
        [x1,x2] = ginput(4);
        
        % Centre du cercle
		L(1) = -( ( (x1(3)^2-x1(2)^2+x2(3)^2-x2(2)^2)/(2*(x2(3)-x2(2))) ) - ...
		      ( (x1(2)^2-x1(1)^2+x2(2)^2-x2(1)^2)/(2*(x2(2)-x2(1))) ) ) / ...
		      ( (x1(2)-x1(1))/(x2(2)-x2(1))-(x1(3)-x1(2))/(x2(3)-x2(2)) );
		L(2) = -L(1)*(x1(2)-x1(1))/(x2(2)-x2(1)) + (x1(2)^2-x1(1)^2+x2(2)^2-x2(1)^2)/(2*(x2(2)-x2(1)));
		% Angle axe queue / horizontale
		L(3) = angle( (x1(4)-L(1))+i*(x2(4)-L(2)) );
		% Echelle
        L(4) = Rsph/sqrt((x1(1)-L(1))^2+(x2(1)-L(2))^2);
        
        % Angles caracteristiques
        theta_m =-asin(LAm/Rsph);
        theta_n = acos(Rcyl/Rsph);
        theta_p = asin((LAp-LAo)/Rpet);
        
        % Longueur des parties du contour
        Lmn = Rsph*(-theta_m+theta_n);
        Lno = LAo-Rsph*sin(theta_n);
        Lop = Rpet*theta_p;
        Ltot= Lmn+Lno+Lop;
        
        % Determination partie haute ou basse par signe du produit vectoriel
        houb = -sign( (x1(2)-L(1))*(x2(4)-L(2))-(x2(2)-L(2))*(x1(4)-L(1)) );
        
        % Sortie
        Xc1 = L;

               
    case 1 % Traces supplementaires ---------------------------------------
        
        
        % Cenre de la sphere
        plot(L(1),L(2),'r+');

        
    case 2 % Calcul simple ------------------------------------------------
        
        
        if houb==1 % Demi-coquille haute
            
            % Partition de x1
            x1po = x1(x1<=Lop/Ltot);
            x1on = x1(x1> Lop/Ltot & x1<=(Lop+Lno)/Ltot);
            x1nm = x1(x1> (Lop+Lno)/Ltot);
            if ~all(x1 - [x1po, x1on, x1nm]==zeros(size(x1)));error;end

            % Coo locales allant de 0 a 1
            xl1po = x1po*Ltot/Lop;
            xl1on = (x1on-Lop/Ltot)*Ltot/Lno;
            xl1nm = (x1nm-(Lop+Lno)/Ltot)*Ltot/Lmn;

            % Conge
            X1po = LAo+Rpet*sin(theta_p*(1-xl1po));
            X2po = Rcyl-Rpet+Rpet*cos(theta_p*(1-xl1po));

            % Partie droite
            X1on = Rsph*sin(theta_n) + Lno*(1-xl1on);
            X2on = Rsph*cos(theta_n)*ones(size(xl1on));

            % Sphere de tete
            X1nm = Rsph*sin(theta_n+(theta_m-theta_n)*xl1nm);
            X2nm = Rsph*cos(theta_n+(theta_m-theta_n)*xl1nm);

            % Assemblage et echelle
            Xc10 = [X1po,X1on,X1nm]/L(4);
            Xc20 = [X2po,X2on,X2nm]/L(4);
            
        elseif houb==-1 % Demi coquille basse
            
            % Partition de x1
            x1mn = x1(x1<=Lmn/Ltot);
            x1no = x1(x1> Lmn/Ltot & x1<=(Lmn+Lno)/Ltot);
            x1op = x1(x1> (Lmn+Lno)/Ltot);
            if ~all(x1 - [x1mn, x1no, x1op]==zeros(size(x1)));error;end
            
            % Coo locales allant de 0 a 1
            xl1mn = x1mn*Ltot/Lmn;
            xl1no = (x1no-Lmn/Ltot)*Ltot/Lno;
            xl1op = (x1op-(Lmn+Lno)/Ltot)*Ltot/Lop;
            
            % Sphere de tete
            X1mn =  Rsph*sin(theta_m+(theta_n-theta_m)*xl1mn);
            X2mn = -Rsph*cos(theta_m+(theta_n-theta_m)*xl1mn);
            
            % Partie droite
            X1no =  Rsph*sin(theta_n) + Lno*xl1no;
            X2no = -Rsph*cos(theta_n)*ones(size(xl1no));
            
            % Conge
            X1op = LAo+Rpet*sin(theta_p*xl1op);
            X2op = -Rcyl+Rpet-Rpet*cos(theta_p*xl1op);
            
            % Assemblage et echelle
            Xc10 = [X1mn,X1no,X1op]/L(4);
            Xc20 = [X2mn,X2no,X2op]/L(4);
           
        end
        
        % Rotation + translation
		Xc1 = Xc10*cos(L(3)) - Xc20*sin(L(3)) + L(1);
		Xc2 = Xc10*sin(L(3)) + Xc20*cos(L(3)) + L(2);
        
        % Derivees
        dXc1dL = [ ones(size(x1)) ; zeros(size(x1)) ; -Xc10*sin(L(3)) - Xc20*cos(L(3)) ; (-Xc10*cos(L(3)) + Xc20*sin(L(3)))/(L(4)) ];
        dXc2dL = [zeros(size(x1)) ;  ones(size(x1)) ;  Xc10*cos(L(3)) - Xc20*sin(L(3)) ; (-Xc10*sin(L(3)) - Xc20*cos(L(3)))/(L(4)) ];
      
        
    case 3 % Informations- ------------------------------------------------
        
        disp('   IBOCS foundry');

end