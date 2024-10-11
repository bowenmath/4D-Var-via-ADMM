function y = f_adj(x, psiphys, d, M, N, dx, dy, dt, AH, BH, optsur)
% Compute the transpose of the gradient of f with repect to x
% with the relation \Delta psiphys = x

% Prediction
lap = laplacian(x,M,N,dx,dy);

xp = x - dt*arakawa(x,psiphys,M,N,dx,dy) ...
    + AH*dt*laplacian(x,M,N,dx,dy) - BH*dt*laplacian(lap,M,N,dx,dy);

% Impose boundary condition on xp
xp = bdcondition(xp,M,N);

% Initial value for SOR(inversepoisson)
initial = zeros(M*N, 1);

[J_x, J_psi] = arakawa_g(x,psiphys,M,N,dx,dy);
[J_xp, J_psi] = arakawa_g(xp,psiphys,M,N,dx,dy);

J_x = transpose(J_x);
J_psi = transpose(J_psi);
J_xp = transpose(J_xp);

    function z = xpTd(subd)
        subr = inversepoisson(J_x*subd,initial,M,N,dx,dy,optsur);
        sublap = laplacian(subd,M,N,dx,dy);
        z = subd - dt*subr - dt*J_psi*subd + AH*dt*sublap - BH*dt*laplacian(sublap,M,N,dx,dy);
    end

r = inversepoisson(J_xp*d,initial,M,N,dx,dy,optsur);
lap = laplacian(d,M,N,dx,dy);
y = d - dt*r - dt*xpTd(J_psi*d) + AH*dt*xpTd(lap) - BH*dt*xpTd(laplacian(lap,M,N,dx,dy));


end