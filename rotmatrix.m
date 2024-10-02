    % creates rotation matrix R = Rz*Ry*Rx given three angles A = [rx, ry, rz]
    function R = rotmatrix(A)
        Rz = [cos(A(3)) -sin(A(3)) 0; sin(A(3)) cos(A(3)) 0; 0 0 1];
        Ry = [cos(A(2)) 0 sin(A(2)); 0 1 0; -sin(A(2)) 0 cos(A(2))];
        Rx = [1 0 0; 0 cos(A(1)) -sin(A(1)); 0 sin(A(1)) cos(A(1))];
        
        R = Rz*Ry*Rx;
    end
   