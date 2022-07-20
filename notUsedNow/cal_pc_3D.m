function eig3 = cal_pc_3D(im3D,sigma3D) 
imSmoothed = zeros(size(im3D));
for i = 1:size(im3D,3)  
   imSmoothed(:,:,i) = imgaussfilt(im3D(:,:,i),sigma3D); 
end
% imSmoothed(maskAll(:) == 0) = nan;
[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(imSmoothed,0);
eig3 = zeros(size(im3D));
a = -(Dxx + Dyy + Dzz);
a = a(:);
b = -(Dyz.^2 + Dxy.^2 + Dxz.^2 - Dxx.*Dzz -Dyy.*Dzz -Dxx.*Dyy);
b = b(:);
c = -(Dxx.*Dyy.*Dzz + 2*Dxy.*Dyz.*Dxz - Dxx.*(Dyz.^2) - Dzz.*(Dxy.^2) - Dyy.*(Dxz.^2));
c = c(:);
p = b - a.^2/3;
q = 2*(a.^3)/27 - (a.*b)/3 + c;
result3 = zeros(length(im3D(:)),3);
result3(:,1) = 2/sqrt(3)*(sqrt(-p)).*cos(1/3.*acos(3*q./(2*p).*sqrt(-3./p)) - 2*pi/3) - a/3;
result3(:,2) = 2/sqrt(3)*(sqrt(-p)).*cos(1/3.*acos(3*q./(2*p).*sqrt(-3./p)) - 4*pi/3) - a/3;
result3(:,3) = 2/sqrt(3)*(sqrt(-p)).*cos(1/3.*acos(3*q./(2*p).*sqrt(-3./p))) - a/3;
eig3(:) = max(result3,[],2);
end


function [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(Volume,Sigma)
%  This function Hessian3D filters the image with an Gaussian kernel
%  followed by calculation of 2nd order gradients, which aprroximates the
%  2nd order derivatives of the image.
% 
% [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,Sigma)
% 
% inputs,
%   I : The image volume, class preferable double or single
%   Sigma : The sigma of the gaussian kernel used. If sigma is zero
%           no gaussian filtering.
%
% outputs,
%   Dxx, Dyy, Dzz, Dxy, Dxz, Dyz: The 2nd derivatives
%
% Function is written by D.Kroon University of Twente (June 2009)
% defaults
if nargin < 2, Sigma = 1; end

if(Sigma>0)
    F=imgaussian(Volume,Sigma);
else
    F=Volume;
end

% Create first and second order diferentiations
Dz=gradient3(F,'z');
Dzz=(gradient3(Dz,'z'));
clear Dz;

Dy=gradient3(F,'y');
Dyy=(gradient3(Dy,'y'));
Dyz=(gradient3(Dy,'z'));
clear Dy;

Dx=gradient3(F,'x');
Dxx=(gradient3(Dx,'x'));
Dxy=(gradient3(Dx,'y'));
Dxz=(gradient3(Dx,'z'));
clear Dx;
end

function D = gradient3(F,option)
% This function does the same as the default matlab "gradient" function
% but with one direction at the time, less cpu and less memory usage.
%
% Example:
%
% Fx = gradient3(F,'x');

[k,l,m] = size(F);
D  = zeros(size(F),class(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:,:) = (F(2,:,:) - F(1,:,:));
    D(k,:,:) = (F(k,:,:) - F(k-1,:,:));
    % Take centered differences on interior points
    D(2:k-1,:,:) = (F(3:k,:,:)-F(1:k-2,:,:))/2;
case 'y'
    D(:,1,:) = (F(:,2,:) - F(:,1,:));
    D(:,l,:) = (F(:,l,:) - F(:,l-1,:));
    D(:,2:l-1,:) = (F(:,3:l,:)-F(:,1:l-2,:))/2;
case 'z'
    D(:,:,1) = (F(:,:,2) - F(:,:,1));
    D(:,:,m) = (F(:,:,m) - F(:,:,m-1));
    D(:,:,2:m-1) = (F(:,:,3:m)-F(:,:,1:m-2))/2;
otherwise
    disp('Unknown option')
end
end