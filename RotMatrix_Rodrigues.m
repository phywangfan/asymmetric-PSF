function [R2] = RotMatrix_Rodrigues(n)
% Rodriues rotation matrix
% input: n is the normal vector (column vector)
n=n./norm(n);
ez = [0,0,1]';
xn = n(1);
yn = n(2);
zn = n(3);

%% based on Matsushima's book
% s = sqrt(xn^2+yn^2);
% R11 = zn - (zn-1)*yn*yn/s;
% R12 = (zn-1)*yn*xn/s;
% R21 = (zn-1)*xn*yn/s;
% R22 = zn-(zn-1)*xn*xn/s;
% 
% R1 = [R11 R12 xn; 
%      R21 R22 yn;
%      -xn -yn zn];
% R1 = R1'; % Rotated global system to local system: loc = R1*glo 
%% based on Rodridues formula
k = cross(ez,n); % rotation axis
if norm(k)~=0
    k = k./norm(k);
end
theta = acos(ez'*n);
k_mat = [0 -k(3) k(2);
        k(3) 0  -k(1);
        -k(2) k(1) 0] ;
R2 = eye(3)+k_mat*sin(theta)+k_mat*k_mat*(1-cos(theta));
R2 = R2'; 
end

