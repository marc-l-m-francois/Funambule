% Gencercle genere le cercle de test pour le papier de 2018
clear

x1 = [1:12];
x2 = [1:12];

x1c = 6.2345;
x2c = 6.6789;
rc   = 10/3;
L = [x1c x2c rc];
L0 = L; % Pour stats

% Raffinement
raff = 2500;

% Profondeur (nombre de bits)
prof = inf;

% Image fine
x1f = [x1(1)-0.5:1/raff:x1(end)+0.5]+1/(2*raff);
x2f = [x2(1)-0.5:1/raff:x2(end)+0.5]+1/(2*raff);
[X1f,X2f] = meshgrid(x1f,x2f);
FF = ones(size(X1f));
FF( (X1f-x1c).^2 + (X2f-x2c).^2 <= rc^2 ) = 0;
% imagesc(FF);colormap gray;axis equal

% Image grosse
[X1,X2] = meshgrid(x1,x2);
F = zeros(size(X1));
for p=1:size(x1,2)
    disp(['Etape ',num2str(p),' / ',num2str(size(x1,2))]);
    for q=1:size(x2,2)
        F(p,q) = mean(mean(FF(x1f>p-0.5 & x1f<=p+0.5 , x2f>q-0.5 & x2f<=q+0.5)));
    end
end


imagesc(F);colormap gray;axis equal
save Cercle_test F prof L0



