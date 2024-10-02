function [rH, iF] = ASM_prop_W(H,di,pp,wlen,radius)
% radius : set a ciryle filter to reconstruct， radius is the value of 0 to 1
if nargin==5
    FH=fftshift(fft2(H));
%     mu=[0,0];
%     Sigma=[10,0;0,10];
%     filter = gauss2D_W(size(FH,2),size(FH,1),pp,mu,Sigma,1); %% gauss filter
    filter = circleW(zeros(size(FH)),radius); % radius is the value of 0 to 1
%     figure(100),imshow(filter);
    FH=FH.*filter;
    H = ifft2(fftshift(FH));
%  figure,surfl(FH.*conj(FH)),shading interp,colormap(gray); xlabel('u');ylabel('v');title('diffractive spectrum')
end
[Ny0,Nx0]=size(H);
H=padarray(H,[Ny0/2,Nx0/2],0,'both');
FH=fftshift(fft2(fftshift(H)));


k=2*pi/wlen;
[Ny,Nx]=size(FH);
Lx=Nx*pp;
Ly=Ny*pp;
fx=linspace(-1/2/pp,1/2/pp-1/Lx,Nx);% 
fy=linspace(1/2/pp,-1/2/pp+1/Ly,Ny);
[fx,fy]=meshgrid(fx,fy');
fx_BL=Nx/2*pp/wlen/abs(di);
fy_BL=Ny/2*pp/wlen/abs(di);  % fy_max;%
fx(abs(fx)>fx_BL)=0;
fy(abs(fy)>fy_BL)=0;
iF=FH.*(exp(1j*k*di*sqrt(1-(wlen*fx).^2-(wlen*fy).^2)));
rH=fftshift(ifft2(fftshift(iF)));
rH=rH(Ny/2-Ny0/2+1:Ny/2+Ny0/2,Nx/2-Nx0/2+1:Nx/2+Nx0/2);
iF=fftshift(fft2(fftshift(rH)));

end

function cir = circleW(ori,radi)
% radi is a radio to y-axis
    [r,c]=size(ori);
    y = 1:r;
    x = 1:c;
    [x,y] = meshgrid(x,y);
    radi = c*radi;
    cir = (x-c/2).^2+(y-r/2).^2<=radi^2;
end