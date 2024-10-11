M = 50;
N = 50;

% Scale of the eddy
L = 25000;

% Resolution with respect to the scale
dx = L/5;
dy = dx;

% Diffusion and dissipation parameters
% Diffusion
%AH=0.1;
AH = 0.01;
% Biharminic diffusion
%BH=10*dx*dy;
BH = dx*dy;

% Plotting scales
omegascale = 0.00001;
psiscale = L*L*omegascale;

% Time step based on estimate of valocity scale and CFL valye
dt = 0.1*dx*dy/L/L/omegascale*10;


time_start = tic;  % Start the timer

% Estimate the optimal relaxation parameter
MM = 1 / sqrt(0.5 / (M * M) + 0.5 / (N * N));  % Calculate the parameter MM
optsur = 2 / (1 + 2*pi / MM);  % Calculate the optimal relaxation parameter

dx = pi / (M-1);  % Spacing in the x direction
dy = pi / (N-1);  % Spacing in the y direction

o = randn(M*N, 1);  % Generate a random vector represents the vorticity

initial = zeros(M*N, 1);  % Initial value for SOR(inversepoisson)

dd = randn(M*N, 1);  % Generate a random direction dd
dd = dd / norm(dd);  % Normalize dd
rr = randn(M*N, 1);  % Generate a random direction rr
rr = rr / norm(rr);  % Normalize rr

eps = 1e-7;

po = inversepoisson(o,initial,M,N,dx,dy,optsur);
po1 = inversepoisson(o+eps*dd,initial,M,N,dx,dy,optsur);
po2 = inversepoisson(o-eps*dd,initial,M,N,dx,dy,optsur);
pdd = inversepoisson(dd,initial,M,N,dx,dy,optsur);

% % The gradient of arakawa
% approx = (arakawa(o+eps*dd,po1,M,N,dx,dy)-arakawa(o-eps*dd,po2,M,N,dx,dy))/(2*eps);
% [A, B] = arakawa_g(o,po,M,N,dx,dy);
% real = A*pdd + B*dd;

% The gradient of f
approx = (f(o+eps*dd,po1,M,N,dx,dy,dt,AH,BH)-f(o-eps*dd,po1,M,N,dx,dy,dt,AH,BH))/(2*eps);
real = f_g(o,po,dd,M,N,dx,dy,dt,AH,BH,optsur);
err = norm(approx-real)/(norm(approx)+norm(real));

% % Calculate the gradient and adjoint operator
% grad = f_g(o,po,dd,M,N,dx,dy,dt,AH,BH,optsur);  % Calculate the gradient
% adj = f_adj(o,po,dd,M,N,dx,dy,dt,AH,BH,optsur);  % Calculate the adjoint operator
% % Calculate the error
% err = dot(grad, dd) - dot(dd, adj)/(norm(grad)+norm(adj));  % Calculate the inner product error

time_end = toc(time_start);  % Stop the timer and record the elapsed time

disp(['Error: ', num2str(err)]);  % Display the error
disp(['Time: ', num2str(time_end), ' seconds']);  % Display the runtime
