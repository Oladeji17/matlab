function [Ce,phi_w,phi_p,phi_w_surf,phi_p_surf] = SLeig_solver(N2,zf,fs)

%% [Ce,phi_w,phi_p,phi_w_surf,phi_p_surf] = SLEIG_SOLVER computes eigenvalues using Ashok & Bhaduria (2009) numerical scheme
%% with either a rigid lid or a free surface boundary condition (Griffiths & Grimshaw 2007, JPO, eq. 4c and Kelly 2016, JPO, eq. 3b)

% Oladeji Siyanbola, USM
% 

%% input
% N2 represents the Brunt Vaisala frequency vector (top to bottom) - mapped
% to layer faces
% zf represents layer faces' z-level vector (approaches -inf with increase in depth). zf(1) = 0;
% fs surface boundary condition: 0 = rigid lid (default), 1 = free surface
% ** for accuracy, N2 and zf should be double-precision **
%% output
% Ce, phi_w and phi_p represent eigen speed, vert. eigen and hor. eigen functions respectively
% phi_w_surf represents the vert. eigen function at the surface (z = 0); zero for rigid lid
% phi_p_surf represents the hor. eigen function at the surface (z = 0)

% Boundary conditions:
%   rigid lid   : phi_w = 0                     at z = 0
%   free surface: phi_w = (Ce^2/g) d(phi_w)/dz  at z = 0   (Griffiths & Grimshaw 2007, eq. 4c; Kelly 2016, eq. 3b)
%   flat bottom : phi_w = 0                     at z = -H
% surface (zeroth) mode included: added aposteriori for rigid lid, obtained
% directly from the EVP for free surface

if nargin<3
    fs = 0;
end

g = 9.81;

dz = -diff(zf); % layer thickness vector (top to bottom)

% Compute parameters
H = nansum_ca(dz); 
N = numel(dz);

% Main solver
A = [];
if fs
    % Unknowns are W at faces 1 (surface) to N (bottom face N+1 has W = 0)
    % Free surface row: same 3-point stencil applied at the surface face using a
    % mirror ghost face at height dz(1) above the surface. The ghost value is
    % eliminated with the centered form of eq. 3b,
    %   (W_ghost - W_2)/(2*dz(1)) = W_1/(g*invCe2)
    % which gives
    %   2/dz(1)^2*(W_2 - W_1) = -invCe2*(N2(1) + 2*g/dz(1))*W_1
    B = diag (-[N2(1) + 2*g/dz(1); N2(2:end-1)]);

    A(1,1) = -2*1/(dz(1)*dz(1));
    A(1,2) =  2*1/(dz(1)*dz(1));
    for i = 2:N-1
        A(i,-1+i) =  2*1/(dz(i-1)*(dz(i-1) + dz(i)));
        A(i, 0+i) = -2*1/(dz(i-1)*dz(i));
        A(i, 1+i) =  2*1/(dz(i)*(dz(i-1) + dz(i)));
    end
    i = i+1;
    A(i,-1+i) =  2*1/(dz(i-1)*(dz(i-1) + dz(i)));
    A(i, 0+i) = -2*1/(dz(i-1)*dz(i));
else
    % Unknowns are W at faces 2 to N (surface and bottom faces have W = 0)
    B = diag (-N2(2:end-1));

    A(1,1) = -2*1/(dz(1)*dz(2));
    A(1,2) = 2*1/(dz(2)*(dz(1) + dz(2)));
    for i = 2:N-2
        A(i,-1+i) =  2*1/(dz(i)*(dz(i) + dz(i+1)));
        A(i, 0+i) = -2*1/(dz(i)*dz(i+1));
        A(i, 1+i) =  2*1/(dz(i+1)*(dz(i) + dz(i+1)));
    end
    i = i+1;
    A(i,-1+i) =  2*1/(dz(i)*(dz(i) + dz(i+1)));
    A(i, 0+i) = -2*1/(dz(i)*dz(i+1));
end

ll = size(A,1);

% Solve the EVP
[W1,invCe2] = eig(A,B);
Ce2 = 1./(diag(invCe2));
Ce = sqrt(Ce2); [Ce,Is]=sort(Ce,'descend');

W1 = W1(:,Is);

%% Get W structure functions
if fs
    phi_w = [W1; zeros(1,ll)];                % W at cell faces
else
    phi_w = [zeros(1,ll); W1; zeros(1,ll)];   % W at cell faces
end
phi_w_surf = phi_w(1,:);                      % W at the surface

%% Get U structure

dph_w = phi_w(2:end,:)-phi_w(1:end-1,:);
dzu = repmat(dz,[1,ll]);
phi_p = dph_w./dzu;

% U structure at the surface (z = 0)
if fs
    % centered derivative at the surface face from eq. 3b: dW/dz = g*W/Ce^2
    % (minus sign since dph_w is taken downward)
    phi_p_surf = -g.*phi_w_surf./(Ce.^2)';
else
    % linear extrapolation of the top two layer centers (equivalent to a
    % 2nd-order one-sided derivative using W = 0 at the surface)
    phi_p_surf = phi_p(1,:) + (phi_p(1,:) - phi_p(2,:)).*dz(1)./(dz(1) + dz(2));
end

%% Normalization (N2-weighted integral plus free surface term; the surface term vanishes for rigid lid)

N2 = N2(1:end-1)./2 + N2(2:end)./2;
phi_w = phi_w(1:end-1,:)./2 + phi_w(2:end,:)./2;

AA = ((sum(phi_w.^2.*repmat(N2,[1,ll]).*dzu,1) + g.*phi_w_surf.^2)./(Ce.^2)')./H;
AA = repmat(AA,[N 1]).^(1/2);
AA(AA==0) = Inf;
phi_w = phi_w./AA;
phi_w_surf = phi_w_surf./AA(1,:);

AA = repmat(sum(phi_p.^2.*dzu,1)./H,[N 1]).^(1/2);
AA(AA==0) = Inf;
phi_p = phi_p./AA;
phi_p_surf = phi_p_surf./AA(1,:);

%% make eig shapes positive at the bottom

flp = phi_p(N,:)<0;
phi_w(:,flp) = -phi_w(:,flp);
phi_w_surf(flp) = -phi_w_surf(flp);
phi_p(:,flp) = -phi_p(:,flp);
phi_p_surf(flp) = -phi_p_surf(flp);

phi_w = -phi_w; % this corrects for z being more negative with depth so that phi_p = dphi_w/dz;
phi_w_surf = -phi_w_surf;

%% including surface (zeroth) mode aposteriori for rigid lid

if ~fs
    phi_w=[zeros([N 1]) phi_w];
    phi_p=[ones([N 1]) phi_p];
    phi_w_surf=[0 phi_w_surf];
    phi_p_surf=[1 phi_p_surf];

    Ce_surf = sqrt(g*H);
    Ce=[Ce_surf; Ce];
end
