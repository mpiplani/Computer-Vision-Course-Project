% Demo script for estimating surface gradient and albedo by PS

clear
close all

addpath('Toolbox/') % Some useful functions for integration of normals

% Data :
% I : nrows x ncols x nimgs
% mask : nrows x ncols
% S (light matrix) : nimgs x 3
load Beethoven

% Note : axis coordinates are Matlab style (with z oriented so that Oxyz is direct)
% 0 -- y
% |
% x


% Show data
figure(1)
subplot(1,4,1)
imagesc(I(:,:,1),[0 255]);
colormap gray;
axis equal
axis off
title('$I^1$','Interpreter','Latex','Fontsize',18)
subplot(1,4,2)
imagesc(I(:,:,2),[0 255]);
axis equal
axis off
title('$I^2$','Interpreter','Latex','Fontsize',18)
subplot(1,4,3)
imagesc(I(:,:,3),[0 255]);
axis equal
axis off
title('$I^3$','Interpreter','Latex','Fontsize',18)
subplot(1,4,4)
imagesc(mask);
axis equal
axis off
title('mask','Interpreter','Latex','Fontsize',18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate normals and albedo by photometric stereo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get masked pixels
indices_mask = find(mask>0);
nb_pixels = length(indices_mask);

% Stack all intensities in a nimgs x (card(mask)) matrix I_vect
[nrows,ncols,nimgs] = size(I);
I_vect = zeros(nimgs,nb_pixels);
for i = 1:nimgs
	Ii = I(:,:,i);
	I_vect(i,:) = transpose(Ii(indices_mask));
end

% Lambert's law reads : I(u,v) = rho(u,v) s.n(u,v)
% Call m(u,v) = rho(u,v) n(u,v) -> I(u,v) = s.m(u,v)
% In matrix form : I_vect = S*M_vect with M_vect : 3 x nb_pixels

% Estimate m = rho * N in each pixel
M_vect = S\I_vect;

% Albedo is rho = \|m\| in each pixel
rho_vect = sqrt(sum(M_vect.^2,1));

% Normal vector is N = m / \|m\| in each pixel
N_vect = M_vect./repmat(rho_vect,[3 1]);

% Form the estimated matrices
rho = zeros(nrows,ncols);
rho(indices_mask) = rho_vect;
N = zeros(nrows,ncols,3);
for i = 1:3
	Ni = N(:,:,i);
	Ni(indices_mask) = N_vect(i,:);
	N(:,:,i) = Ni;
end

% Show the estimated albedo and normal map
figure(2)
subplot(1,2,1)
imagesc(rho,[0 255])
axis equal
axis off
colormap gray
title('Albedo','Interpreter','Latex','Fontsize',18)
subplot(1,2,2)
imshow(uint8(128*(N+1)))
axis equal
axis off
title('Normals','Interpreter','Latex','Fontsize',18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normals to gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assuming orthographic projection, we have
% N = [-z_x,-z_y,1]/sqrt(z_x^2+z_y^2+1)
% -> z_x = -N1/N3 and z_y = -N2/N3
% We denote this estimation of the gradient : (p,q)
p = -N(:,:,1)./N(:,:,3);
q = -N(:,:,2)./N(:,:,3);

% Show the estimated albedo and normal map
figure(3)
subplot(1,2,1)
imagesc(p)
axis equal
axis off
colormap jet
title('$p$','Interpreter','Latex','Fontsize',18)
subplot(1,2,2)
imagesc(q)
axis equal
axis off
title('$q$','Interpreter','Latex','Fontsize',18)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-squares integration on masked data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calling integrate_with_prior without arguments comes to integrating with least-squares

z1 = integrate_with_prior(p,q,[],[],mask);

% Show result
figure(5)
surfl(fliplr(z1),[0 90])
shading flat
axis equal
view(130,30)
colormap gray
camlight headlight
camlight(0,90)
axis off
title('Least-squares integration','Interpreter','Latex','Fontsize',18)

save results_PS_Beethoven
