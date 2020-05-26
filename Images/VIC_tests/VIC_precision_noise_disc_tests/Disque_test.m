% Reference disc generator
% Used in 2019 publication


% CHOICE : PIXEL REFINEMENT


n = 100;


% INITS


% Disc center coordinates
C1 = 7.2345;
C2 = 7.6789;
RA  = 10/3;
% Image size
F = ones(14);
% Noise %
noise = 10;    


% CALCULUS


% Fine image FF
xf1 = [1/(2*n):1/n:size(F,2)]+0.5;
xf2 = [1/(2*n):1/n:size(F,1)]+0.5;
[XF1,XF2] = meshgrid(xf1,xf2);
% Condition of appartenance
R2 = RA^2;
FF = (XF1-C1).^2+(XF2-C2).^2 < R2;
n2 = n^2;
for p=1:size(F,1)
    for q=1:size(F,2)
        F(p,q) = 1-nnz(FF(n*(p-1)+[1:n],n*(q-1)+[1:n]))/n2;
    end
end


% ADDITIVE NOISE


F = F + noise*randn(size(F))/100;


% DEBUG


if false
    figure(1);clf;
    imagesc(F,[0 1]);hold on
    plot(XF1(FF),XF2(FF),'ro','MarkerSize',12,'MarkerFaceColor','r')
    plot(XF1(~FF),XF2(~FF),'r+')
    P = 2*pi*RA;
    theta = 2*pi*[0:1/n:P]/P;
    % x1 = [0:1/nraff:P]/P;
    Xc1 = C1+RA*cos(theta);
    Xc2 = C2+RA*sin(theta);
    plot(Xc1,Xc2,'b-')
    [XG1,XG2] = meshgrid([0:size(F,1)]+0.5,[0:size(F,2)]+0.5);
    plot(XG2,XG1,'g-');plot(XG2',XG1','g-');
    colormap gray
    axis equal tight

    for p=1:size(F,1)
        for q=1:size(F,2)
            text(q,p,num2str(nnz(FF(n*(p-1)+[1:n],n*(q-1)+[1:n]))),'FontSize',24)
        end
    end
end


% RECORDING


comp_name = ['X1=',num2str(C1),',X2=',num2str(C2),',R=',num2str(RA),',noise=',num2str(noise)]

imwrite(F,['DISC_',comp_name,'.png'],'PNG','BitDepth',16);


% STATISTICS


if 0==1
    % Distance signee
    delta = sqrt((Xc1-C1).^2 + (Xc2-C2).^2) - RA;
    mean(delta)
    std(delta)
end