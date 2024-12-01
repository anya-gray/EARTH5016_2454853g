% create model setup from TIFF image

function [units, D, Nz] = ModelFromImage(filename, n_units, W, Nx)

% read in RGB image from TIFF file
img = double(importdata(filename));

% get size of image
[p,q,k] = size(img);

% identify units by k-means clustering analysis
[Ic,Fc] = kmeans(reshape(img,p*q,k),9,'MaxIter',1e3,'Replicates',5);

% sort clusters for ascending SiO2 content
[~,isort]   = sort(Fc(:,1)-sum(Fc,2),'descend');
Ics = zeros(size(Ic));
for ic = 1:n_units
    Ics(Ic==isort(ic)) = ic;
end
imgc = reshape(Ics,p,q);

% interpolate from original dimensions to target model size
D  = W*p/q;
Nz = floor(Nx*p/q);

ho  = W/q;
xco = ho/2:ho:W-ho/2;
zco = ho/2:ho:D-ho/2;
[Xco,Zco] = meshgrid(xco,zco);

h   = W/Nx;
xc  = h/2:h:W-h/2;
zc  = h/2:h:D-h/2;
[Xc,Zc] = meshgrid(xc,zc);

imgi = interp2(Xco,Zco,imgc,Xc,Zc);
% figure(1); clf
% 
% 
% subplot(2,1,1)
% imagesc(xco,zco,imgc); axis equal tight; colorbar
% subplot(2,1,2)
% im = imagesc(xc,zc,imgi); axis equal tight; colorbar


units = uint8(imgi);

end