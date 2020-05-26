%function testim(F);
%
if 0==1
    %close(1)
    figure(1);clf;
    Range_mon_trace(1,1,1,2,4);
    set(1,'Name','Image etat initial F','MenuBar','none',...
        'GraphicsSmoothing','off','ToolBar','none');
    set(gca,'Position',[0 0 1 1]);
    colormap(gray);hold on
end

tic;
imagesc(F);axis equal tight off;hold on
plot([-50 100],[-50 100],'r-')

toc

% tic;imshow(F8);toc 10x plus lent !


% InnerPosition: [12 511 408 495]