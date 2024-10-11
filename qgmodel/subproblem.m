function x1 = subproblem(x0, x1, x2, p0, p1, p2, lmda1, lmda2, xs, M, N, dx, dy, dt, rho, AH, BH, optsur, l, m, n)
% Solve the linearized subproblems of multi-block ADMM
% Use l,m,n to control which term appears in the subproblem
% x1 is the variable to update

eta = 0.1;  % The parameter \eta for linearization in ADMM
para = 1 + eta*(2*l + rho*m);

d = rho*x2 - rho*f(x1,p1,M,N,dx,dy,dt,AH,BH) - lmda2;
b = 2*l*xs + m*(rho*f(x0,p0,M,N,dx,dy,dt,AH,BH) + lmda1) + n*f_adj(x1,p1,d,M,N,dx,dy,dt,AH,BH,optsur);

x1 = (x1 + eta*b)/para;

% Impose boundary condition on x1
x1 = bdcondition(x1,M,N);

end