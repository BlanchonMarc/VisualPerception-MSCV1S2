%plotimagepoints    plot image points
%    [M,nomia,nomib]=loadtestimages(image)
%
% by X. Armangue
% (c) Mr3D - University of Girona, September 2002
%
function plotimagepoints(M,nomia,nomib)

if ~isempty(nomia),
    figure;
    map=colormap;
    imga=imread(nomia);
    imshow(imga);
    hold on;
    plot(M(1,:),M(2,:),'ro');
    zoom on;
    figure;
    colormap(map);
    imgb=imread(nomib);
    imshow(imgb);
    hold on;
    plot(M(3,:),M(4,:),'ro');
    zoom on;
    
    figure;
    colormap(map);
    imga=imread(nomia);
    imshow(imga);
    hold on;
    h=quiver(M(1,:),M(2,:),M(3,:)-M(1,:),M(4,:)-M(2,:),0);
    set(h,'color','w','LineWidth',1);
    zoom on;
else 
    figure;
    hold on;
    plot(M(1,:),M(2,:),'ro');
    zoom on;
    figure;
    hold on;
    plot(M(3,:),M(4,:),'ro');
    zoom on;
    
    figure;
    hold on;
    h=quiver(M(1,:),M(2,:),M(3,:)-M(1,:),M(4,:)-M(2,:),0);
    set(h,'color','k','LineWidth',1);
    
    %	plot(M(1,:),M(2,:),'r+');
    %	plot(M(3,:),M(4,:),'g+');
    %	for i=1:size(M,2),
    %		plot([M(1,i) M(3,i)],[M(2,i) M(4,i)],'b');
    %	end
    zoom on;
end
