function[varfftl] = df_1D13_out(var)
%% [varfftl] = DF_1D13_OUT(VAR) 
%% Oladeji Siyanbola, USM, 2020-12-11
%% computes leftward propagating decomposed fields
%% input
% var represents the field of interest (2-D); 1D in space and time
%% output
% varfftl is the rleftward propagating filtered field 
[ll,mm] = size(var);
ab = fft(var,[],2); varft = ab(:,1:ceil(mm/2));varft(:,2:end) = 2.*varft(:,2:end);% fourier transform in time
ac = fft(varft,[],1);u1fxl = ac(1:ceil(ll/2),:); u1fxl(1,:) = 0.*u1fxl(1,:); % fourier transform in space
u1fxl = cat(1,u1fxl,zeros(floor(ll/2),ceil(mm/2)));%u1fxr = conj(u1fxr);
u1iftl = ifft(u1fxl,[],1); % inverse fourier transform in space
u1iftl = cat(2,u1iftl,zeros(ll,floor(mm/2)));
varfftl = ifft(u1iftl,[],2);varfftl = real(varfftl); % inverse fourier transform in time
end

% [ll,mm] = size(var);
% ab = fft(var,[],2); varft = ab(:,1:ceil(mm/2));varft(:,1) = 0; varft(:,2:end) = 2.*varft(:,2:end);% fourier transform in time
% ac = fft(varft,[],1);u1fxl = ac(1:ceil(ll/2),:); u1fxl(1,:) = 0.5*u1fxl(1,:); % fourier transform in space
% u1fxl = cat(1,u1fxl,zeros(floor(ll/2),ceil(mm/2)));%u1fxr = conj(u1fxr);
% u1iftl = ifft(u1fxl,[],1); % inverse fourier transform in space
% u1iftl = cat(2,u1iftl,zeros(ll,floor(mm/2)));
% varfftl = ifft(u1iftl,[],2);varfftl = real(varfftl); % inverse fourier transform in time
% end




