%% This program is for computing the asymmetric PSF used for the point cloud object: two pieces of pawns
% mm is unit
clear all
%% Basic parameters
pp = 3.8e-3;      % pixel pitch
wlen = 532e-6;    % wavelength
k  = 2*pi/wlen;
Nx0 = 1920*2;         % resolution of the global plane
Ny0 = 1080*2;
Lx0 = Nx0*pp;  % the canvas size of the global plane
Ly0 = Ny0*pp;
%% rotation matrix: define a tilted plane
% ----obtain by rotate angle----
alpha_y = deg2rad(-22.2); % rotate around y
alpha_x = deg2rad(-0.1); % rotate around x
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
Nx = round(Lx/pp); Nx = Nx+mod(Nx,2); % ensure an even number
Ny = round(Ly/pp); Ny = Ny+mod(Ny,2); % ensure an even number
x = pp*(-Nx/2+1:Nx/2); 
y = pp*(-Ny/2+1:Ny/2); 
[x, y] =meshgrid(x,-y); % th
%% 
dx = R_g2l(1)*x + R_g2l(2)*y;
dy = R_g2l(4)*x + R_g2l(5)*y;
dz = R_g2l(7)*x + R_g2l(8)*y;
%% load the point-cloud object
model_name = '.\chess_two_2';
% [vertex,vertexNor,vertex_tex, faces,index_nor, index_tex] = readOBJ([model_name,'.obj']);
data = readply([model_name '.ply']);
vertex = getfields(data.vertex, 2,'x','y','z');
% [vertex, ~] = rotateOBJ(vertex, 1,[ 90 0 0]);
[vertex,~] = resizeOBJ(vertex,[],Lx0,Ly0);

vertex(:,3) = vertex(:,3) +200;
figure(1),plot3(vertex(:,1),vertex(:,2),vertex(:,3),'k.','markersize',5), xlabel('x'),ylabel('y'),zlabel('z'), %view([0 0 1])
axis off


%% plot 3D object 
% TR=triangulation(double(faces(:,1:3)),vertex);
% figure,trimesh(TR,'EdgeColor','r'), axis([-Lx0/2 Lx0/2 -Ly0/2 Ly0/2]),  xlabel('x'),ylabel('y'),zlabel('z')
%% computing wavefield for each points
wavefield_loc = zeros(Ny,Nx);   %% wavefield in the local plane
tic
for np = 1:size(vertex,1) 
    %% source point
    x0 = vertex(np,1);
    y0 = vertex(np,2);
    z0 = vertex(np,3);     % propagation depth
    %% Filter mask of aPSF in the tilted plane
    theta = asin(wlen/2/pp);  % Limited diffraction angle by the pixel pitch
    interP = intersection_LinePlane([x0,y0,z0], theta,R_g2l);
    N_smp = ceil(max([abs(interP(1:2,:)),[Lx;Ly]],[],2)/pp); N_smp = N_smp+mod(N_smp,2);
    mask_int = round(abs(interP(1:2,:)/pp+[N_smp(1); -N_smp(2)]));
    if ~isequal(mask_int(2,1),mask_int(2,2),mask_int(2,3),mask_int(2,4))&& ~isequal(mask_int(1,1),mask_int(1,2),mask_int(1,3),mask_int(1,4))  % make sure four vertexes of mask
        mask = poly2mask(mask_int(1,:),mask_int(2,:),N_smp(2)*2,N_smp(1)*2);
        mask = mask(N_smp(2)-Ny/2+1:N_smp(2)+Ny/2,N_smp(1)-Nx/2+1:N_smp(1)+Nx/2);
    else
        continue
    end
    % figure(2),imshow(mask,[])
    index = find(mask==1);
    %% proposed aPSF formula
    dx = R_g2l(1)*x(mask) + R_g2l(2)*y(mask);
    dy = R_g2l(4)*x(mask) + R_g2l(5)*y(mask);
    dz = R_g2l(7)*x(mask) + R_g2l(8)*y(mask);
    r = sqrt((x0 - dx).^2 + (y0 - dy).^2 + (z0 - dz).^2)*(-1)^(z0<0);
    aPSF = exp(1i*k * (r + R_g2l(7)*x(mask) + R_g2l(8)*y(mask)))./(1j*wlen*r);

    wavefield_loc(index) = wavefield_loc(index) + aPSF;
    % figure(1), imagesc(angle(wf_loc)), axis off,colormap(parula(10)), title('Proposed input local plane')
end
figure,imshow(angle(wavefield_loc),[])
SLM = wavefield_loc(Ny/2-Ny0/2+1:Ny/2+Ny0/2,Nx/2-Nx0/2+1:Nx/2+Nx0/2);
phs_SLM = (angle(SLM)+pi)/2/pi;

%% off-axis
% [H_off,R_cpx] = offaxis(SLM,8,wlen,pp,'~0');
% H_off= rescale(angle(H_off),0,1);
imwrite(phs_SLM,['.\aPSF_image','angY=',...
    num2str(rad2deg(alpha_y)),' angX=',num2str(rad2deg(alpha_x)),'D200mm_big2','.bmp'])

% SLM_phs = rescale(angle(SLM),0,1);
% imwrite(SLM_phs,['E:\OneDrive\8 Tilted PSF\Experiments\UPO SLM\Objects with multi points\phiy_-22° phix_-0.7° D-6mmTo6mm','.bmp'])

% SLM_amp = rescale(SLM.*conj(SLM),0,1);
% imwrite(SLM_amp,['E:\OneDrive\8 Tilted PSF\Experiments\UPO SLM\Objects with multi points\amp_phiy_-22° phix_-0.7° D-6mmTo6mm','.bmp'])
%% rotate aPSF of the local system to the global system
FT = fftshift(fft2(fftshift(wavefield_loc)));
FS = rot_poly2D(FT, R_g2l', [Nx,Ny]*pp, wlen, false);  %  rotmatrix_w(0, rad2deg(alpha))
wavefile_glb = ifftshift(ifft2(ifftshift(FS))); % wavefield in the global plane
wavefile_glb = wavefile_glb(Ny/2-Ny0/2+1:Ny/2+Ny0/2,Nx/2-Nx0/2+1:Nx/2+Nx0/2); 
toc
figure(2), imshow(real(wavefile_glb), []), %title('Output resulting proposed psf');

for d= z0-10:z0+10 %linspace(0.1e-2,z0*1.5,30)  %z0-0.005:0.001:z0+0.005
    [rH,~] = ASM_prop_W(wavefile_glb,-d,pp,wlen);
    reco = abs(rH);
    reco = reco./Ny0/Nx0;
    figure(101),imshow(reco,[]),  title(['focus at',num2str(d),'mm'])  
    
%     [rH,~] = ASM_prop_W(wavefield_loc,-d,pp,wlen);
%     reco = abs(rH);
%     reco = reco./Ny0/Nx0;
%     figure(102),imshow(reco,[]),  title(['focus at',num2str(d),'mm']) 
end

