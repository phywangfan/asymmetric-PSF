%% This program is for computing the asymmetric PSF 
%% It computes for one point in Fig.2 and 5*4 points in Fig.4
% mm is unit
clear all
%% Basic parameters
pp = 3.8e-3;      % pixel pitch
wlen = 532e-6;    % wavelength
k  = 2*pi/wlen;
Nx0 = 3840;         % resolution of the global plane
Ny0 = 2160; 
Lx0 = Nx0*pp;  % the canvas size of the global plane
Ly0 = Ny0*pp;
%% rotation matrix
% ----obtain by rotate angle----
% alpha_y = deg2rad(36); % rotate around y to generate Fig.2
% alpha_x = deg2rad(45); % rotate around x to generate Fig.2
alpha_y = deg2rad(22); % rotate around y to generate Fig.4
alpha_x = deg2rad(0); % rotate around x to generate Fig.4

[R_angle] = rotmatrix([alpha_x alpha_y 0]); % rotation from global to local
n = [0 0 1]*R_angle; % the normal vector of the tilted plane
% ----obtain by the normal vector---- 
% one can user-define a normal vector
% n = [ 0. 5 0.2 0.8]; % for example
% n = n./norm(n); % normalization

R_g2l = RotMatrix_Rodrigues(n'); % rotation matrix obtained by Rodrigues
%% sampling in the local system
alpha_x2 = [n(1),n(3)]./norm([n(1),n(3)]);
alpha_y2 = [n(2),n(3)]./norm([n(2),n(3)]);
Lx = abs(Lx0./alpha_x2(2)); if isinf(Lx) Lx =Lx0; end  % the x size of the local convas 
Ly = abs(Ly0./alpha_y2(2)); if isinf(Ly) Ly =Ly0; end  % the y size of the local convas 
Nx = ceil(Lx/pp); Nx = Nx+mod(Nx,2); % ensure an even number
Ny = ceil(Ly/pp); Ny = Ny+mod(Ny,2); % ensure an even number
x = pp*(-Nx/2+1:Nx/2); 
y = pp*(-Ny/2+1:Ny/2); 
[x, y] =meshgrid(x,-y); % the sampling grid in the local system (tilted plane)  

%%  computing wavefield for single points in Fig. 2
% X = -0.8;
% Y = 0.2;
% z0 = 20;
%% computing wavefield for multiplie points in Fig. 4
X = linspace(-Lx0/6,Lx0/6,5);
Y = linspace(-Lx0/4,Lx0/4,3);
z0 = 195;

wf_loc = zeros(Ny,Nx);  %% wavefield in the local plane
for n = 1:length(X)  %  5 points is defined
    for m = 1:length(Y)
%% source point
x0 = X(n);   % coordinate of point;
y0 = Y(m);    % coordinate of point;
z0 = z0 + 5;  % propagation depth
%% Filter mask of aPSF in the tilted plane
theta = asin(wlen/2/pp);  % Limited diffraction angle by the pixel pitch
interP = intersection_LinePlane([x0,y0,z0], theta,R_g2l);
N_smp = ceil(max([abs(interP(1:2,:)),[Lx;Ly]],[],2)/pp); 
mask_int = round(abs(interP(1:2,:)/pp+[N_smp(1); -N_smp(2)]));
mask = poly2mask(mask_int(1,:),mask_int(2,:),N_smp(2)*2,N_smp(1)*2);
mask = mask(N_smp(2)-Ny/2+1:N_smp(2)+Ny/2,N_smp(1)-Nx/2+1:N_smp(1)+Nx/2);
% figure(2),imshow(mask,[])
index = find(mask==1);
%% proposed aPSF formula
dx = R_g2l(1)*x(mask) + R_g2l(2)*y(mask);
dy = R_g2l(4)*x(mask) + R_g2l(5)*y(mask);
dz = R_g2l(7)*x(mask) + R_g2l(8)*y(mask);
r = sqrt((x0 - dx).^2 + (y0 - dy).^2 + (z0 - dz).^2)*(-1)^(z0<0);
aPSF = exp(1i*k * (r + R_g2l(7)*x(mask) + R_g2l(8)*y(mask)))./(1j*wlen*r); % proposed aPSF 
wf_loc(index) = wf_loc(index) + aPSF; 

% imwrite(aPSF_phs,['.\aPSF_image\','phi_y=',num2str(rad2deg(alpha_y)),' phi_x=',num2str(rad2deg(alpha_y)),'.bmp'])
figure(1), imshow(angle(wf_loc),[]),% axis off,colormap(parula(10)), title('Proposed input local plane')
    end
end
aPSF_phs = (angle(wf_loc)+pi)/2/pi;
figure(2), imshow(aPSF_phs),  title('Proposed input local plane')
% imwrite(aPSF_phs,['E:\OneDrive\8 Tilted PSF\Experiments\UPO SLM\multi_points\','5pointArray angY=',num2str(rad2deg(alpha_y)),' angX=',num2str(rad2deg(alpha_x)),'_200mm_Spacing5mm.png'])



%% rotate aPSF of the local system to the global system
FT = fftshift(fft2(fftshift(wf_loc)));
FS = rot_poly2D(FT, R_g2l', [Nx,Ny]*pp, wlen, false);  %  rotmatrix_w(0, rad2deg(alpha))
wf_glb = ifftshift(ifft2(ifftshift(FS))); % wavefield in the global plane
wf_glb = wf_glb(Ny/2-Ny0/2+1:Ny/2+Ny0/2,Nx/2-Nx0/2+1:Nx/2+Nx0/2); 

figure(2), imshow(real(wf_glb), []), %title('Output resulting proposed psf');

for d= z0-5:z0+5 %linspace(0.1e-2,z0*1.5,30)  %z0-0.005:0.001:z0+0.005
    [rH,~] = ASM_prop_W(wf_glb,-d,pp,wlen);
    figure(101),imshow(abs(rH),[]),  title(['focus at',num2str(d),'mm'])  
end
