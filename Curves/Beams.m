function [Xc1,Xc2,dXc1dL,dXc2dL] = Beams(L,x1,rien);
% Adresses the adimenstionnal problems in beam_equations folder
% in which you can add v = f(x) equations as you want
%
% Except linear displacement fields whcich are impossible 
% because they are similar to the placement of the reference straight line
%
% L(1)/beamlength = angle theta from x1 (image) to x (reference undeformed beam axis)
% L(2)            = distance delta along y from A0 to A
% L(3)... are the mechanical parameters used in the unitary beam equation


global A1 A2 B1 B2          % user fixed points
global beam_type            % type of mechanical problem
global beamlength           % Used for conditionning the matrix
global myfontsize

switch nargin
    
    case 0 % Manual settings
        disp('    Click first and last points')
        [AB1,AB2] = ginput(2);  
        A1 = AB1(1);
        A2 = AB2(1);
        B1 = AB1(2);
        B2 = AB2(2);
        % Plot
        plot(AB1,AB2,'r--+','MarkerSize',2*myfontsize);
        text(A1,A2,' 0','FontSize',2*myfontsize,'Color','r');
        text(B1,B2,' 1','FontSize',2*myfontsize,'Color','r');
         
        % Length of the beam at rest
        beamlength = sqrt((B2-A2)^2+(B1-A1)^2);
        % Initial angle t from X1 to the mechanical frame
        L(1) = angle( (B1-A1) + 1i*(B2-A2) ) * beamlength;
        % Initial distance d from A to O along y
        L(2) = 0;
        % User choice : type of curve
        disp('Choose the type of curve between:');
        cd Curves/BEAMS;a = dir;cpt=0;
        for p=1:size(a,1)
            if strcmp(a(p).name(end),'m') && ~a(p).isdir && size(a(p).name,2)>2
                if strcmp(a(p).name(end-1),'.')
                    cpt = cpt+1;
                    disp(['    ',num2str(cpt),' ',a(p).name(1:end-2)]);num(cpt) = p;
                end
            end
        end
        rep = input('Your choice: ');
        beam_type = a(num(rep)).name(1:end-2);cd ../..
        % Initial mechanical parameters
        L = [L,feval(beam_type) * beamlength];
        % Return the defaut values
        Xc1 = L;
        
       
    case 1 % Supplementary graphics ---------------------------------------
        
        % retriveing data
        t = L(1) / beamlength;
        d = L(2);
        % reference, undeformed beam equation
        x1 = [0 1];
        y1 = [0 0];
        c = cos(t);
        s = sin(t);
        l     = (B1-A1)*c + (B2-A2)*s;
        Xc1 = A1 + l*x1*c - (d + l*y1)*s;
        Xc2 = A2 + l*x1*s + (d + l*y1)*c;
        plot(Xc1,Xc2,'r--','LineWidth',1);
        % User cliked end points
        plot(A1,A2,'b+','MarkerSize',12);
        plot(B1,B2,'b+','MarkerSize',12);
        % Calculated end points
        text(Xc1(  1),Xc2(  1),' 0','Color','r','FontSize',myfontsize);
        text(Xc1(end),Xc2(end),' 1','Color','r','FontSize',myfontsize);
        
        % Internal supplementary graphics
        [x1,name] = feval(beam_type,L(3:end) / beamlength);
        y1 = 0;
        Xc1 = A1 + l*x1*c - (d + l*y1)*s;
        Xc2 = A2 + l*x1*s + (d + l*y1)*c;
        plot(Xc1,Xc2,'r+','LineWidth',1);
        text(Xc1,Xc2,[' ',name],'Color','r','FontSize',myfontsize);
        
       
    case 2 % Calculus -----------------------------------------------------
        
        % retriveing data
        t = L(1) / beamlength;
        d = L(2);
        
        % Curve equation
        [y1 ,yp1] = feval(beam_type,L(3:end) / beamlength,x1); % the unitary equation of the beam
        c = cos(t);
        s = sin(t);
        l   = (B1-A1)*c + (B2-A2)*s;
        Xc1 = A1 + l*x1*c - (d + l*y1)*s;
        Xc2 = A2 + l*x1*s + (d + l*y1)*c;
        
        % Derivatives / t, d
        dXc1dL = zeros(max(size(L)),size(x1,2));
        dXc2dL = zeros(max(size(L)),size(x1,2));
        dldt =-(B1-A1)*s + (B2-A2)*c;
        dXc1dL(1:2,:) = [ ((dldt*(x1*c -  y1*s) - l*x1*s - (d+l*y1)*c)) / beamlength;
                         -s*ones(size(x1))];
        dXc2dL(1:2,:) = [ ((dldt*(x1*s +  y1*c) + l*x1*c - (d+l*y1)*s)) / beamlength;
                          c*ones(size(x1))];
        
        % Derivatives / Ki
        for p=1:size(L,2)-2
            dXc1dL(p+2,:) =-l*yp1(p,:)*s / beamlength;
            dXc2dL(p+2,:) = l*yp1(p,:)*c / beamlength;
        end
        
    case 3 % Informations -------------------------------------------------
        
        eval(['help ',beam_type]);
        disp(L(3:end) / beamlength);

end