function y = f_g(x, psiphys, d, M, N, dx, dy, dt, AH, BH, optsur)
% The gradient of f(x, psiphys)
% Prediction
lap = laplacian(x,M,N,dx,dy);

xp = x - dt*arakawa(x,psiphys,M,N,dx,dy) ...
    + AH*dt*lap - BH*dt*laplacian(lap,M,N,dx,dy);

% Impose boundary condition on xp
xp = bdcondition(xp,M,N);

% Initial value for SOR(inversepoisson)
initial = zeros(M*N, 1);

[J_x, J_psi] = arakawa_g(x,psiphys,M,N,dx,dy);
[J_xp, J_psi] = arakawa_g(xp,psiphys,M,N,dx,dy);

lap = laplacian(d,M,N,dx,dy);
pd = inversepoisson(d,initial,M,N,dx,dy,optsur);
r = d - dt*J_x*pd - dt*J_psi*d...
   + AH*dt*lap - BH*dt*laplacian(lap,M,N,dx,dy);

lap = laplacian(r,M,N,dx,dy);
y = d - dt*J_xp*pd - dt*J_psi*r...
   + AH*dt*lap - BH*dt*laplacian(lap,M,N,dx,dy);

end