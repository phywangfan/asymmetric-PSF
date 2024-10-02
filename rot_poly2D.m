function [Y, D, J] = rot_poly2D(X, R, L, wlen, fw, enable_shift)
% fw = true or false
% when fw is true, perform the rotation from the global to the local
% when fw is false, perform the rotation from the lcoal to the global

% Rotate 2D plane spectrum (input/output in fourier domain)
% assumptions: fftshifted (central frequency is zero), rotating around midpoint

% Computes the spectrum of an arbitrary rotated plane using rotation matrix R
% Returns the shifting amount and the jacobian as well if desired.
% Optional frequency shift (on by default)

[r,c] = size(X);
U = (-c/(2*L(1)):1/L(1):(c-1)/(2*L(1)));
V = (-r/(2*L(2)):1/L(2):(r-1)/(2*L(2)));
[U, V] = meshgrid(U,-V);  % sampled frequencies
if ~exist('enable_shift', 'var'), enable_shift = true; end

%%
if fw  %%   rotating from the global to the local
    D = enable_shift * [R(7)/wlen, R(8)/wlen];
    
    W = real(sqrt(wlen^(-2) -  (U-D(1)).^2 - (V-D(2)).^2));           % z-component k-vector
%     W0 = real(sqrt(wlen^(-2) -  (U).^2 - (V).^2));           % z-component k-vector
%     W = W0 + R(9)/wlen;
    % Y = interp2(U,V,X, R(1)*(U+D(1))+R(2)*(V+D(2))+R(3)*W, R(4)*(U+D(1))+R(5)*(V+D(2))+R(6)*W, 'cubic',0) ;%.* (1 - (W==0)); % resample frequencies
    
    J = ((R(2)*R(6)-R(3)*R(5))*(U-D(1)) + (R(3)*R(4)-R(1)*R(6))*(V-D(2))) ./W + (R(1)*R(5)-R(2)*R(4)); J(isinf(J)) = 0; % jacobian   
    Y = sqrt(J).*interp2(U,V,X, R(1)*(U-D(1))+R(2)*(V-D(2))-R(3)*W, R(4)*(U-D(1))+R(5)*(V-D(2))-R(6)*W, 'spline',0) ;%.* (1 - (W==0)); % resample frequencies
    
else %% rotating from the local to the global
    
    D = enable_shift * [R(7)/wlen, R(8)/wlen];
    W = real(sqrt(wlen^(-2) -  U.^2 - V.^2));           % z-component k-vector   
%     Y = interp2(U,V,X, R(1)*U+R(2)*V+R(3)*W-D(1), R(4)*U+R(5)*V+R(6)*W-D(2), 'cubic',0) ;%.* (1 - (W==0)); % resample frequencies

    J = ((R(4)*R(8)-R(7)*R(5))*U + (R(7)*R(2)-R(1)*R(8))*V) ./W + (R(1)*R(5)-R(4)*R(2)); J(isinf(J)) = 0; % jacobian
    Y = interp2(U,V, X.*sqrt(J), R(1)*U+R(4)*V+R(7)*W-D(1), R(2)*U+R(5)*V+R(8)*W-D(2), 'spline',0) ;%.* (1 - (W==0)); % resample frequencies
    
%     fx_p = R(1)*U+R(4)*V+R(7)*W-D(1);
%     fy_p = R(2)*U+R(5)*V+R(8)*W-D(2);
%     figure,plot3(U,V,zeros(size(U)),'r.', fx_p,fy_p,zeros(size(U)),'b.' )

    
end




