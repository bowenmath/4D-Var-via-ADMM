function x = f(x, psiphys, M, N, dx, dy, dt, AH, BH)
% Compute the iteration function for prediction-correction time scheme x_{i+1} = f(x_i)

% Prediction
lap = laplacian(x,M,N,dx,dy);

xp = x - dt*arakawa(x,psiphys,M,N,dx,dy) ...
    + AH*dt*lap - BH*dt*laplacian(lap,M,N,dx,dy);

% Impose boundary condition on xp
xp = bdcondition(xp,M,N);

% Correction
lap = laplacian(xp,M,N,dx,dy);

x = x - dt*arakawa(xp,psiphys,M,N,dx,dy) ...
   + AH*dt*lap - BH*dt*laplacian(lap,M,N,dx,dy);

% Impose boundary condition on x
x = bdcondition(x,M,N);
end



