function [x1,x2] = EBHOMHPP(choix,li,X1,X2,p)
% -------------------------------------------------------------------------
% Transfo de type Euler Bernoulli homogene en Petites Perturbations
% li(1) translation selon X repere meca / repere image F
% li(2) translation selon Y repere meca / repere image F
% li(3) rotation repere meca / repere image F
% li(4) courbure en repere meca (echelle pixel de F)
% li(5) coefficient de Poisson
% li(6) translation selon X repere image G / repere meca
% li(7) translation selon Y repere image G / repere meca
% -------------------------------------------------------------------------
% x1 x2	position image finale
% X1 X2	position image initiale (entiers) global
% li 	parametres du champ : On doit avoir li(n)=0 ==> x1=X1 et x2=X2
% p     indice  du lambda pour derivation
% -------------------------------------------------------------------------



global sizeF


% % LIGNES DE TEST
% clear
% M = 11;
%  [X1,X2] = meshgrid([1:M]);
%  li = [3 0 0 1/50 0 3 0];
%  X11=X1;
%  X22=X2;
%  sizeF = size(X1);
% choix = 'champ';
% % ------- FINTEST
% 

switch(choix)
    
    
    case 'nombre_param'
        
        
        x1 = 7;
        x2 = [];
        
        
    case 'infos'
        
        
        % Inits -----------------------------------------------------------
        
        X01 = li(1);
        X02 = li(2);
        TH0 = li(3);                % Rotation
        E   = 0;                    % Echelle 1
        rho = li(4)/max(sizeF);     % Courbure
        nu  = li(5);                % Poisson
        x01 = li(6);
        x02 = li(7);
        th0 = TH0;
        e   = E;
        
        % Displays
        infos = {...
             %['Courbure rho = ',num2str(rho),' pixels-1'];
             ['Rayon de courbure = ',num2str(abs(1/rho)),' pixels'];
             ['Coefficient de Poisson = ',num2str(nu)];
             };
        disp(infos);
        
        % Retour ligne moyenne pour trace
        X01A = X01+(sizeF(1)+1)/2; % Pour mettre par defaut le milieu au
        X02A = X02+(sizeF(2)+1)/2; % centre de l'image
        X1 = [0.5:sizeF(2)-0.5];
        X2 = tan(TH0)*(X1-X01A)+X02A;
        XI1 = ( (X1-X01A)*cos(TH0)+(X2-X02A)*sin(TH0) )/(1+E);
        XI2 = (-(X1-X01A)*sin(TH0)+(X2-X02A)*cos(TH0) )/(1+E);
        % champ meca etat final
        if abs(rho)<1E250
            xi1 =  XI1 - rho*XI1.*XI2;
            xi2 =  XI2 + rho*XI1.^2/2 + rho*nu*XI2.^2/2;
        else
            xi1 =  XI1;
            xi2 =  XI2;
        end
        % inits
        x01 = x01+(sizeF(1)+1)/2;
        x02 = x02+(sizeF(2)+1)/2;
        % Champ image G
        x1 = x01 + (1+e)*( xi1*cos(th0)-xi2*sin(th0) );
        x2 = x02 + (1+e)*( xi1*sin(th0)+xi2*cos(th0) );
        % Sortie en repere de F et en repere de g
        x1 = [X1;x1];
        x2 = [X2;x2];
        
        
    case 'champ'
        
        
        % Inits -----------------------------------------------------------
        
        X01 = li(1);
        X02 = li(2);
        TH0 = li(3);                % Rotation
        E   = 0;                    % Echelle 1
        rho = li(4)/max(sizeF);     % Courbure
        nu  = li(5);                % Poisson
        x01 = li(6);
        x02 = li(7);
        th0 = TH0;
        e   = E;
        
        % Plongement inverse du repere de F vers repere mecanique M -------
        
        X01A = X01+(sizeF(1)+1)/2; % Pour mettre par defaut le milieu au
        X02A = X02+(sizeF(2)+1)/2; % centre de l'image
        XI1 = ( (X1-X01A)*cos(TH0)+(X2-X02A)*sin(TH0) )/(1+E);
        XI2 = (-(X1-X01A)*sin(TH0)+(X2-X02A)*cos(TH0) )/(1+E);
        
        % Champ mecanique de M vers m =====================================
        
        % champ meca etat final
        if abs(rho)<1E250
            xi1 =  XI1 - rho*XI1.*XI2;
            xi2 =  XI2 + rho*XI1.^2/2 + rho*nu*XI2.^2/2;
        else
            xi1 =  XI1;
            xi2 =  XI2;
        end
        
        % Plongement du repere mecanique m vers le repere de G ------------
        
        % inits
        x01 = x01+(sizeF(1)+1)/2;
        x02 = x02+(sizeF(2)+1)/2;
        % Champ image G
        x1 = x01 + (1+e)*( xi1*cos(th0)-xi2*sin(th0) );
        x2 = x02 + (1+e)*( xi1*sin(th0)+xi2*cos(th0) );
        
        
    case 'derivees'

        
        % Inits -----------------------------------------------------------
        
        X01 = li(1);
        X02 = li(2);
        TH0 = li(3);                % Rotation
        E   = 0;                    % Echelle 1
        rho = li(4)/max(sizeF);     % Courbure
        nu  = li(5);                % Poisson
        x01 = li(6);
        x02 = li(7);
        th0 = TH0;
        e   = E;
        
        % Parametres du plongement X01 X02 TH0 E
        if p==1 % X01
            % Derivee champ mecanique dxidXI ==============================
            X01A = X01+(sizeF(1)+1)/2;
            X02A = X02+(sizeF(2)+1)/2;
            XI1 = ( (X1-X01A)*cos(TH0)+(X2-X02A)*sin(TH0) )/(1+E);
            XI2 = (-(X1-X01A)*sin(TH0)+(X2-X02A)*cos(TH0) )/(1+E);
            dxidXI11 = ones(size(XI1)) - rho*XI2;
            dxidXI12 = - rho*XI1;
            dxidXI21 =   rho*XI1;
            dxidXI22 = ones(size(XI2)) + rho*nu*XI2;
            % derivee plongement dXi / [dX01]
            dXi1dli = -cos(TH0)/(1+E);
            dXi2dli =  sin(TH0)/(1+E);
            % Derivee plongement inverse dx / dxi
            dxdxi11 =  (1+e)*cos(th0);
            dxdxi12 = -(1+e)*sin(th0);
            dxdxi21 =  (1+e)*sin(th0);
            dxdxi22 =  (1+e)*cos(th0);         
            % Assemblage
            vx1 = dxidXI11.*dXi1dli + dxidXI12.*dXi2dli;
            vx2 = dxidXI21.*dXi1dli + dxidXI22.*dXi2dli;
            x1 = dxdxi11.*vx1 + dxdxi12.*vx2;
            x2 = dxdxi21.*vx1 + dxdxi22.*vx2;
        elseif p==2 % X02
            % Derivee champ mecanique dxidXI ==============================
            X01A = X01+(sizeF(1)+1)/2;
            X02A = X02+(sizeF(2)+1)/2;
            XI1 = ( (X1-X01A)*cos(TH0)+(X2-X02A)*sin(TH0) )/(1+E);
            XI2 = (-(X1-X01A)*sin(TH0)+(X2-X02A)*cos(TH0) )/(1+E);
            dxidXI11 = ones(size(XI1)) - rho*XI2;
            dxidXI12 = - rho*XI1;
            dxidXI21 =   rho*XI1;
            dxidXI22 = ones(size(XI2)) + rho*nu*XI2;
            % derivee plongement dXi / [dX02]
            dXi1dli = -sin(TH0)/(1+E);
            dXi2dli = -cos(TH0)/(1+E);
            % Derivee plongement inverse dx / dxi
            dxdxi11 =  (1+e)*cos(th0);
            dxdxi12 = -(1+e)*sin(th0);
            dxdxi21 =  (1+e)*sin(th0);
            dxdxi22 =  (1+e)*cos(th0);         
            % Assemblage
            vx1 = dxidXI11.*dXi1dli + dxidXI12.*dXi2dli;
            vx2 = dxidXI21.*dXi1dli + dxidXI22.*dXi2dli;
            x1 = dxdxi11.*vx1 + dxdxi12.*vx2;
            x2 = dxdxi21.*vx1 + dxdxi22.*vx2;
        elseif p==3 % TH0
            % Derivee champ mecanique dxidXI ==============================
            X01A = X01+(sizeF(1)+1)/2;
            X02A = X02+(sizeF(2)+1)/2;
            XI1 = ( (X1-X01A)*cos(TH0)+(X2-X02A)*sin(TH0) )/(1+E);
            XI2 = (-(X1-X01A)*sin(TH0)+(X2-X02A)*cos(TH0) )/(1+E);
            dxidXI11 = ones(size(XI1)) - rho*XI2;
            dxidXI12 = - rho*XI1;
            dxidXI21 =   rho*XI1;
            dxidXI22 = ones(size(XI2)) + rho*nu*XI2;
            % derivee plongement dXi / [dTH0]
            dXi1dli = (-(X1-X01A)*sin(TH0)+(X2-X02A)*cos(TH0) )/(1+E);
            dXi2dli = (-(X1-X01A)*cos(TH0)-(X2-X02A)*sin(TH0) )/(1+E);
            % Derivee plongement inverse dx / dxi
            dxdxi11 =  (1+e)*cos(th0);
            dxdxi12 = -(1+e)*sin(th0);
            dxdxi21 =  (1+e)*sin(th0);
            dxdxi22 =  (1+e)*cos(th0);
            % A cause de l'egalite th0 = TH0
             if abs(rho)<1E250
                xi1 =  XI1 - rho*XI1.*XI2;
                xi2 =  XI2 + rho*XI1.^2/2 + rho*nu*XI2.^2/2;
            else
                xi1 =  XI1;
                xi2 =  XI2;
            end
            dx1dth0 = (-xi1*sin(th0)-xi2*cos(th0))*(1+e);
            dx2dth0 = ( xi1*cos(th0)-xi2*sin(th0))*(1+e);
            % Assemblage
            vx1 = dxidXI11.*dXi1dli + dxidXI12.*dXi2dli;
            vx2 = dxidXI21.*dXi1dli + dxidXI22.*dXi2dli;
            x1 = (dxdxi11.*vx1 + dxdxi12.*vx2 + dx1dth0) ;
            x2 = (dxdxi21.*vx1 + dxdxi22.*vx2 + dx2dth0) ;
        % derivee plongement dXi / [dE] non utilisee ici
            % dXi1dli   = (-(X1-X01)*cos(TH0)-(X2-X02)*sin(TH0) )/(1+E)^2;
            % dXi2dli   = ( (X1-X01)*sin(TH0)-(X2-X02)*cos(TH0) )/(1+E)^2;
        elseif p==4 % rho
            % Derivee champ mecanique dxidlp ==============================
            X01A = X01+(sizeF(1)+1)/2;
            X02A = X02+(sizeF(2)+1)/2;
            XI1 = ( (X1-X01A)*cos(TH0)+(X2-X02A)*sin(TH0) )/(1+E);
            XI2 = (-(X1-X01A)*sin(TH0)+(X2-X02A)*cos(TH0) )/(1+E);
            dxi1dlp = -XI1.*XI2/max(sizeF);
            dxi2dlp =  (XI1.^2/2 + nu*XI2.^2/2)/max(sizeF);
            % Derivee plongement inverse dx / dxi
            dxdxi11 =  (1+e)*cos(th0);
            dxdxi12 = -(1+e)*sin(th0);
            dxdxi21 =  (1+e)*sin(th0);
            dxdxi22 =  (1+e)*cos(th0);
            % Assemblage
            x1 = dxdxi11*dxi1dlp + dxdxi12*dxi2dlp;
            x2 = dxdxi21*dxi1dlp + dxdxi22*dxi2dlp;   
        elseif p==5 % nu
            % Derivee champ mecanique dxidlp ==============================
            X01A = X01+(sizeF(1)+1)/2;
            X02A = X02+(sizeF(2)+1)/2;
            XI1 = ( (X1-X01A)*cos(TH0)+(X2-X02A)*sin(TH0) )/(1+E);
            XI2 = (-(X1-X01A)*sin(TH0)+(X2-X02A)*cos(TH0) )/(1+E);
            dxi1dlp = zeros(size(XI1));
            dxi2dlp =  rho*XI2.^2/2;
            % Derivee plongement inverse dx / dxi
            dxdxi11 =  (1+e)*cos(th0);
            dxdxi12 = -(1+e)*sin(th0);
            dxdxi21 =  (1+e)*sin(th0);
            dxdxi22 =  (1+e)*cos(th0);
            % Assemblage
            x1 = dxdxi11*dxi1dlp + dxdxi12*dxi2dlp;
            x2 = dxdxi21*dxi1dlp + dxdxi22*dxi2dlp;   
        elseif p==6 % x01
            x1 = ones(size(X1));
            x2 = zeros(size(X1));
        elseif p==7 % x02
            x1 = zeros(size(X1));
            x2 = ones(size(X1)); 
        end
end
 

 
% TEST ----------------------------------
% figure(1)
% clf;plot(x1,x2,'k-');hold on;plot(x1',x2','k-');
% quiver(X11,X22,x1-X11,x2-X22,0);
% grid on;axis equal;
% 
% 
% F11 = gradient(x1 ) ;
% F12 = gradient(x1')'; 
% F21 = gradient(x2 ) ;
% F22 = gradient(x2')';
% % Green-Lagrange E = (F'*F - I)/2
% E11 = (F11.^2 + F21.^2 - 1)/2;
% E12 = (F11.*F12 + F21.*F22)/2;
% %E21 = (F12.*F11 + F22.*F21)/2;
% E22 = (F12.^2 + F22.^2 - 1)/2;
% % PETITES DEFOS H = F-I, EPS = (H+H')/2
% 
% figure(2);clf
% %pcolor(x1,x2,sqrt((x1-X11).^2+(x2-X22).^2))
% %pcolor(x1,x2,E11+E22)
% %pcolor(x1,x2,-E22./E11)
% %pcolor(x1,x2,(F12+F21)/2);%-(F22-1)/(F11-1))
% 
% % HPP
% %pcolor(x2-X22);title('u2')
% pcolor(x1,x2,F11+F22-2);title('E11+E22')
% shading interp
% colorbar
% axis equal
% colormap jet
% hold on
% xx1=x1;
% xx1(abs(E22)>1E-4)=NaN;
% plot(xx1(:),x2(:),'k+');
% 
% % Valeurs theorik
% figure(3);clf
% pcolor(x1,x2,-X22*rho+nu*(X22-(sizeF(2)+1)/2)*rho);title('E11+E22 tho')
% shading interp
% colorbar
% axis equal
% colormap jet
% hold on
% FIN TEST ----------------------------------


