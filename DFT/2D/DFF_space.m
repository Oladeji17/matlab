function [var_out] = DFF_space(var,theta1,theta2,land_mask)
%% Oladeji Siyanbola, USM, 2025-10-23
%% modification of Gong et al. 2020 DFF_HRET function
%% generates waves whose wavenumber vector vectors fall between theta1 and theta2
%% input
% var represents the field of interest (3-D); 2D in space [nx,ny] and 1D in time
% theta1 and theta2 are the lower and the upper angular limits of the
% angular sector respectively
%
% ensure nx and ny to be odd
%% output
% var_out is the filtered field



%% example

% a0 = 2; T = 742; sig = 2*pi/T;
% x = 0:1:160; nx = length(x);
% y = 0:1:320; ny = length(y);
% Lx = 161; k = 2*pi/Lx;
% Ly = 161; l = 2*pi/Ly;
% t = 0:1:741;nt = length(t);
%
% a_1q = zeros(nx,ny,nt);a_2q = a_1q;a_3q = a_1q;a_4q = a_1q; % q stands
% for quadrant
%
% for ii = 1:nx
%     for jj = 1:ny
%         for kk = 1:nt
%     a_1q(ii,jj,kk) = 0.5*a0*cos(sig*4*t(kk) - 4*k*x(ii) - 4*l*y(jj));
%     a_2q(ii,jj,kk) = 0.5*a0*cos(sig*4*t(kk) - 4*k*x(ii) + 4*l*y(jj));
%     a_3q(ii,jj,kk) = a0*cos(sig*4*t(kk) + 4*k*x(ii) + 4*l*y(jj));
%     a_4q(ii,jj,kk) = 0.5*a0*cos(sig*4*t(kk) + 4*k*x(ii) - 4*l*y(jj));
%         end
%     end
% end
%
% atot = a_1q + a_2q + a_3q + a_4q;
%
% theta1 = 0; theta2 = 90; % for first quadrant
% land_mask = ones(nx,ny);
% [var_theta] = space_filt(atot,theta1,theta2,land_mask);
%
% [X,Y] = meshgrid(x,y);
%
% figh = figure;
% for ii = 1:nt
% pcolor(X,Y,var_theta(:,:,ii)')
% shading flat
% caxis([-4,4])
% colormap(bluewhitered(256))
% colorbar
% movieVector(ii) = getframe(figh);
% end

%%
%%
var(isnan(var)) = 0;
[nx,ny,nt] = size(var);

% fft in time
var_t = fft(var,[],3); var_t = var_t(:,:,1:ceil(nt/2));var_t(:,:,2:end) = 2.*var_t(:,:,2:end);

% fft in space

var_xy = zeros(size(var_t));
for i = 1:ceil(nt/2)
    var_xy(:,:,i) = fftshift(fft2(var_t(:,:,i)));
end

% mask function

kx = fliplr(fftshift(2*pi./nx.*[0:floor(nx/2), -floor(nx/2):-1]));
ky = fliplr(fftshift(2*pi./ny.*[0:floor(ny/2), -floor(ny/2):-1]));

[KY, KX] = meshgrid(ky,kx);

[theta, radii]=cart2pol(KX,KY); % cartesian to polar coordinate transformation
theta=rad2deg(theta); % radian to degrees conversion
theta(theta<0)=360+theta(theta<0);

mask = theta>=theta1 & theta<=theta2;

% inverse transform in space_2D

var_theta = zeros(size(var_xy));
for i = 1:ceil(nt/2)
    var_theta(:,:,i) = ifft2(ifftshift(var_xy(:,:,i).*mask));
end
var_theta = cat(3,var_theta,zeros(nx,ny,floor(nt/2)));

% inverse transform in time

var_theta = ifft(var_theta,[],3);var_theta = real(var_theta);

% masking out land regions
var_out = zeros(size(var_theta));
for i = 1:nt
    var1 = var_theta(:,:,i);
    var1(land_mask == 0) = NaN;
    var_out(:,:,i) = var1;
end

end