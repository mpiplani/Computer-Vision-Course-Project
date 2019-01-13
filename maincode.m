%% Tried to implement a part of paper : Enforcing Integrability by Error Correction using l1 minimisation
% Authors : Dikpal Reddy, Amit Agrawal and Rama chellappa
% Published in: 2009 IEEE Conference on Computer Vision and Pattern Recognition

tic
clear all;
close all;
clc;

global RMSE_TH;
global maxZ;

ADD_OUTLIERS = 1;
ADD_NOISE = 1;
RMSE_TH = 0.01;

%% generate ramp peak surface (im)

H = 128;
W = 64;

im = zeros(64,64);
[x,y] = meshgrid(-8:8,-8:8);

tt = (abs(x)-8).*(abs(y)-8);
[h,w] = size(tt);

st = 5;
im(st:st+h-1,st:st+w-1) = im(st:st+h-1,st:st+w-1) + 0.5*tt;

st = 10;
im(st:st+h-1,30+st:30+st+w-1) = im(st:st+h-1,30+st:30+st+w-1) + 0.5*tt;

st = 35;
im(st:st+h-1,st:st+w-1) = im(st:st+h-1,st:st+w-1) + 0.15*tt;

for j = 10:50
    im(50:60,j) = (j-10)/4;
end

im = imfilter(im,fspecial('gaussian',6,1),'symmetric');

maxZ = max(im(:));
[H,W] = size(im);
[ogx,ogy] = calculate_gradients(im,0,0);

%% add noise in gradients
if(ADD_NOISE)
    tt = sqrt(ogx.^2 + ogy.^2);
    sigma = 5*max(tt(:))/100;
    clear tt
else
    sigma = 0;
end

gx = ogx + sigma*randn(H,W);
gy = ogy + sigma*randn(H,W);

%% add uniformly distributed outliers in gradients
if(ADD_OUTLIERS)
    fac = 3;
    outlier_x = rand(H,W) > 0.90;
    outlier_x(:,end) = 0;
    
    gx = gx + fac*outlier_x.*(2*(rand(H,W)>0.5)-1);
    
    outlier_y = rand(H,W) > 0.90;
    outlier_y(end,:) = 0;
    
    gy = gy + fac*outlier_y.*(2*(rand(H,W)>0.5)-1);
    outlier_x = double(outlier_x);
    outlier_y = double(outlier_y);
    disp(sprintf('Gx outliers = %d',sum(outlier_x(:))));
    disp(sprintf('Gy outliers = %d',sum(outlier_y(:))));
end
gx(:,end) = 0;
gy(end,:) = 0;

disp('============================================');
disp('Least squares solution by solving Poisson Equation')
r_ls = poisson_solver(gx,gy);
r_ls = r_ls - min(r_ls(:));

disp('============================================');
disp('l1 minimisation using l1 magic')
% g=(gx+gy)./2;
% j = 1:H+1;
% k = 1:W+1;
% % gyy(j+1,k) = gy(j+1,k) - gy(j,k);
% % gxx(j,k+1) = gx(j,k+1) - gx(j,k);
f = gx + gy;
% f = f(2:end-1,2:end-1);
% r_ls1 = poisson_solverl1(gx,gy);
% r_ls1 = r_ls1 - min(r_ls1(:));

xpx = l1eq_pd(im, gx, [], r_ls, 1e-3);
xpy = l1eq_pd(im, gy, [], r_ls, 1e-3);
xp=(xpx+xpy)./2;
close all;

mydisplay(im);
title('original');
axis on;

mydisplay(r_ls);
axis on;
title('Least Squares');

mydisplay(xp);
axis on;
title('L1 minimisation');
toc