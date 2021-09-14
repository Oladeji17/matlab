function[varfftr] = df_1D13(var)
%% [varfftr] = DF_1D22 
%% Oladeji Siyanbola, USM, 2020-12-11
%% computes rightward propagating decomposed fields
%% input
% var represents the field of interest (2-D); 1D in space and time
%% output
% varfftr is the rightward propagating filtered field 
[ll,mm] = size(var);
ab = fft(var,[],2); varft = ab(:,1:ceil(mm/2));varft(:,2:end) = 2.*varft(:,2:end); varft = conj(varft);% fourier transform in time
ac = fft(varft,[],1);u1fxr = ac(1:ceil(ll/2),:); u1fxr(1,:) = 0.*u1fxr(1,:); % fourier transform in space
u1fxr = cat(1,u1fxr,zeros(floor(ll/2),ceil(mm/2)));%u1fxr = conj(u1fxr);
u1iftr = ifft(u1fxr,[],1); % inverse fourier transform in space
u1iftr = cat(2,u1iftr,zeros(ll,floor(mm/2)));
varfftr = ifft(conj(u1iftr),[],2);varfftr = real(varfftr); % inverse fourier transform in time
end


%% leftward component
% [ll,mm] = size(var);
% ab = fft(var,[],2); varft = ab(:,1:mm/2);varft(:,2:end) = 2.*varft(:,2:end); % fourier transform in time
% ac = fft(varft,[],1);u1fxr = ac(1:ceil(ll/2),:); u1fxr(1,:) = 0.*u1fxr(1,:); % fourier transform in space
% u1fxr = cat(1,u1fxr,zeros(floor(ll/2),mm/2));%u1fxr = conj(u1fxr);
% u1iftr = ifft(u1fxr,[],1); % inverse fourier transform in space
% u1iftr = cat(2,u1iftr,zeros(ll,mm/2));
% varfftr = ifft(u1iftr,[],2);varfftr = real(varfftr); % inverse fourier transform in time