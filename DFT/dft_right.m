function[var_right] = dft_right(var)
%% Oladeji Siyanbola, USM, 2022-06-15
%% separates the rightward propagating wave
%% input
% var represents the input field (2-D); 1D in space and time
% 1st dimension is space (delta_x = const.), 2nd dimension is time (delta_t = const.)
%% output
% var_right is the rightward propagating wave 

[ll,mm] = size(var);

% fourier transform in time
var_om = fft(var,[],2);var_om = var_om(:,1:ceil(mm/2));var_om(:,2:end) = 2.*var_om(:,2:end);var_om = conj(var_om);

% fourier transform in space
var_k = fft(var_om,[],1);var_k = var_k(1:ceil(ll/2),:);var_k(1,:) = 0.5*var_k(1,:); 

% inverse fourier transform in space
var_k = cat(1,var_k,zeros(floor(ll/2),ceil(mm/2)));
var_x = ifft(var_k,[],1); 

% inverse fourier transform in time
var_x = cat(2,var_x,zeros(ll,floor(mm/2)));
var_t = ifft(conj(var_x),[],2);var_right = real(var_t); 
end