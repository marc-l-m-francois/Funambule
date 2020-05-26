function [speckle] = Speckle(F)

try
    F = imread(F);
    F = mean(double(F)/255,3);
end

[h,b] = histcounts(F);
cu = cumsum(h)/numel(F);
medianF = mean(diff(b))*sum(cu<.5);

mF = mean(F(:));

stats = regionprops(F>=mF,'EquivDiameter');
region = [(1:length(stats))',vertcat(stats.EquivDiameter)];

speckle.MoyenneTache = mean(region(region(:,2)<100,2));
speckle.EcarttypeTache = std(region(region(:,2)<100,2));
speckle.Densite = sum((F(:)>=mF))/numel(F);
speckle.Contraste = sqrt((1/numel(F))*sum((F(:)-mean(F(:))).^2));

speckle.MaxF = max(F(:));
speckle.MinF = min(F(:));

speckle.MedianF = medianF;
speckle.MoyenneF = mean(F(:));

figure(1)
subplot(1,2,1)
imagesc(F)
axis equal
colormap gray
colorbar
title 'Image'
xlabel 'X (pixel)'
ylabel 'Y (pixel)'
subplot(1,2,2)
imagesc(F>=mF)
axis equal
colormap gray
colorbar
title 'Threshold image'
xlabel 'X (pixel)'
ylabel 'Y (pixel)'

L = 2*10*ceil(speckle.MoyenneTache);

map.mst = zeros(size(F));
map.sst = map.mst;
map.max = map.mst;
map.min = map.mst;
map.dens = map.mst;
map.cont = map.mst;
for i=1:L:(size(F,1)-L)
    for j=(1:L:size(F,2)-L)
        id = i:(i+L-1);
        jd = j:(j+L-1);
        
        Fij = F(id,jd);
        
        stats = regionprops(Fij>=mF,'EquivDiameter');
        region = [(1:length(stats))',vertcat(stats.EquivDiameter)];

        try
            mstij = mean(region(region(:,2)<100,2));
            sstij = std(region(region(:,2)<100,2));
        catch
            mstij = 0;
            sstij = 0;
        end
        maxij = max(Fij(:));
        minij = min(Fij(:));
        densij = sum((Fij(:)>=mF))/numel(Fij);
        contij = sqrt((1/numel(Fij))*sum((Fij(:)-mean(Fij(:))).^2));
        
        map.mst(id,jd) = mstij;
        map.sst(id,jd) = sstij;
        map.max(id,jd) = maxij;
        map.min(id,jd) = minij;
        map.dens(id,jd) = densij;
        map.cont(id,jd) = contij;
        
    end
end

figure(2)
subplot(3,2,1)
imagesc(map.mst)
axis equal
colormap jet
colorbar
caxis([0 9]);
title 'Moyenne des taches'
xlabel 'X (pixel)'
ylabel 'Y (pixel)'
subplot(3,2,3)
imagesc(map.sst)
axis equal
colormap jet
colorbar 
caxis([0 9]);
title 'Ecart-type des taches'
xlabel 'X (pixel)'
ylabel 'Y (pixel)'
subplot(3,2,2)
imagesc(map.max)
axis equal
colormap jet
colorbar
caxis([0 1]);
title 'Niveau de gris max'
xlabel 'X (pixel)'
ylabel 'Y (pixel)'
subplot(3,2,4)
imagesc(map.min)
axis equal
colormap jet
colorbar
title 'Niveau de gris min'
xlabel 'X (pixel)'
ylabel 'Y (pixel)'
caxis([0 1]);
subplot(3,2,5)
imagesc(map.dens)
axis equal
colormap jet
colorbar
caxis([0 1]);
title 'Densité'
xlabel 'X (pixel)'
ylabel 'Y (pixel)'
subplot(3,2,6)
imagesc(map.cont)
axis equal
colormap jet
colorbar
caxis([0 .5]);
title 'Contraste'
xlabel 'X (pixel)'
ylabel 'Y (pixel)'
pause(.1)

path(path,'Routines');
path(path,'Routines_graphiques');
path(path,'Champs');

precision = 0.001;
cptmax=200;
type_interp = 'linear';
champ='DEFHOM';
F0 = F(10:end,10:end);
g0 = F(1:end-9,1:end-9);
DF = ones(size(F0))==1;
[X1,X2] = meshgrid([1:size(F0,2)],[1:size(F0,1)]);
rgs = [16 8 0];
[x1,x2,infos,li] = feval(champ,[],X1,X2,0);
li = zeros(size(li));
for rg = rgs
    F = F0;
    g = g0;
    if rg>0
        H = fspecial('gaussian',6*ceil(rg)+1,rg/2);
        F = filter2(H,F);
        g = filter2(H,g);
    end
    
    % INITS CORRELATION ---------------------------------------------------
    % Gradients de l'image F
    [gradF10,gradF20]  = gradient(F);
    % Init compteur de boucles et ecart
    cpt = 1;
    % Champ actuel
    [x1,x2] = feval(champ,li,X1,X2,0);
    % Image G = image g interpolee
    G = interp2(X1,X2,g,x1,x2,type_interp,NaN);
    % Son domaine de definition intersection avec celui de F
    DG = isfinite(G) & DF;
    % Surface = nombre de pixels correles
    S = sum(sum(DG));
    % Ecart GG-F
    G_F = G-F;
    % On annule l'effet des points hors domaine
    G_F(~DG) = 0;
    % Erreur actuelle
    Psi(cpt) = sum(sum(G_F.^2))/S;
    % Variables internes au Newton
    % gradF_dxdl = zeros(size(F,1),size(F,2),size(li,1));
    gradF_dxdlp = zeros(size(F,1),size(F,2));
    gradF_dxdlq = zeros(size(F,1),size(F,2));
    M = zeros(max(size(li)));
    V = zeros(max(size(li)),1);
    
    dli = inf*ones(size(li));
    cpt = 0;
    while max(abs(dli))>precision && cpt<cptmax && (sum(DG(:))>0)
        cpt = cpt+1;
        % Elimination des points hors domaine
        gradF1 = gradF10;
        gradF2 = gradF20;
        gradF1(~DG) = 0;
        gradF2(~DG) = 0;
        % Calcul matrice et vecteur
        for p=1:max(size(li))
            % Calcul de la p-ieme derivee en li
            [dx1dL,dx2dL] = feval(champ,li,X1,X2,p);
            gradF_dxdlp  =  gradF1.*dx1dL +  gradF2.*dx2dL;
            % Vecteur
            V(p) = -sum(sum(gradF_dxdlp.*G_F));
            % Matrice
            for q=1:p
                if q~=p
                    [dx1dL,dx2dL] = feval(champ,li,X1,X2,q);
                    gradF_dxdlq  =  gradF1.*dx1dL +  gradF2.*dx2dL;
                    M(p,q) = sum(sum( gradF_dxdlp.*gradF_dxdlq ));
                    M(q,p) = M(p,q);
                else
                    M(p,q) = sum(sum( gradF_dxdlp.*gradF_dxdlp));
                end
            end
        end
        
        % Verification de la determination de la matrice
        detM = abs(det(M));
        
        dli = (M\V)';
        li = li + dli;
        
        [x1,x2] = feval(champ,li,X1,X2,0);
        % Image G = image g interpolee
        G = interp2(X1,X2,g,x1,x2,type_interp,NaN);
        % Son domaine de definition intersection avec celui de F
        DG = isfinite(G) & DF;
        % Surface = nombre de pixels correles
        S = sum(sum(DG));
        % Ecart GG-F
        G_F = G-F;
        % On annule l'effet des points hors domaine
        G_F(~DG) = 0;
        % Erreur actuelle
        Psi(cpt) = sum(sum(G_F.^2))/S;
    end
    Psi = [];
end


speckle.MoyenneAutoCorrelation = mean(abs(G_F(:)));
speckle.MaxAutoCorrelation = max(abs(G_F(:)));
end

