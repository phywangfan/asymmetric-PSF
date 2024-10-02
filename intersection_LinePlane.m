function P_inlocal = intersection_LinePlane(pnt, theta,R)
% pnt: source point
% theta: diffraction angle
% R : rotation matrix

%% plane equation
n = [0 0 1]*R;  % n=[nx,ny,nz] colum vector
%% intersection in the orthogonal plane
x0 = pnt(1);
y0 = pnt(2);
z0 = pnt(3);

inter_orth = [z0*tan(theta) z0*tan(theta) 0]';
inter_orths = inter_orth.*[1 1 1; -1 1 1; -1 -1 1; 1 -1 1]'+[x0 y0 0]';  % four intersection points [first second third forth] quartile
coff = inter_orths-[x0 y0 z0]';
%% t
X = n * coff;
Y = -n * [x0 y0 z0]';
t = Y./X;

ff = @(i) coff(:,i).*t(i) + [x0,y0,z0]';
P1 = ff(1); %% intersection point
P2 = ff(2); 
P3 = ff(3); 
P4 = ff(4); 

P_inGlobal = [P1, P2, P3, P4];  % intersection points in the global system

% conv = [1 -1 1; 1 1 -1;1 1 1];
P_inlocal = R*P_inGlobal;  % intersection points in the local system



 

end