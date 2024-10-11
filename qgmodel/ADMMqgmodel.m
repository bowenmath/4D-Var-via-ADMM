%---------------------------------------------------------
% This program uses linearized multi-block ADMM to recover
% the 2D vorticity concentration phenomenon.
%---------------------------------------------------------
% By Bowen Li
% October 7, 2024
%---------------------------------------------------------

% M points in along channel direction. Two of them are for periodic conditions
% N points across channel, take an odd number of points (psi on center)
M = 20;
N = 20;

MN = M*N;

% Number of time steps
NT = 500;

% Scale of the eddy
% L = 25000;
L = 1;

% Resolution with respect to the scale
dx = L/5;
dy = dx;

% Diffusion and dissipation parameters
% Diffusion
%AH=0.1;
AH = 0;
% Biharminic diffusion
%BH=10*dx*dy;
BH = 0.001*dx*dy;


% Plotting scales
% omegascale = 0.00001;
omegascale = 1;
%psiscale = L*L*omegascale;

% Time step based on estimate of valocity scale and CFL valye
dt = 3*dx*dy/L/L/omegascale;
T = dt*NT;

% Initial eddy, perturbed by mode 2 structure
% for ic=1:M
%     for jc=1:N
%         xx(ic,jc)=(ic-M/2)*dx;
%         yy(ic,jc)=(jc-N/2)*dy;
%         rr=sqrt(xx(ic,jc)^2+yy(ic,jc)^2)/L;
%         th=atan2(yy(ic,jc),xx(ic,jc));
%         r=rr*(1+0.03*cos(2*th));
%         psiphys(ic,jc)=exp(-r*r*r)*L^2*omegascale;
%     end
% end

% Set the random seed
seed = 3;
rng(seed);

% Set the initial value
% psiphys = 1000*random('Normal',0,1,M*N,1);
% psiphys = psiphys - mean(psiphys);
omphys = 5*omegascale*random('Normal',0,1,M*N,1);
omphys = omphys - mean(omphys);


% Create coordinates of a grid
xx = zeros(MN, 1);
yy = zeros(MN, 1);
for i=1:M
    for j=1:N
        xx((i-1)*N+j, 1) = i*dx;
        yy((i-1)*N+j, 1) = j*dy;
    end
end

% Estimate optimal overrelaxation parameter
MM = 1/(sqrt(0.5/(M*M)+0.5/(N*N)));

optsur = 2/(1+pi/MM);

% Use customized display color
brcol = bluered(64);

xst = zeros(MN, NT+1);  % The true value of vorticity simulation
xs = zeros(MN, NT+1);  % The noised observation of vorticity
pst = zeros(MN, NT+1);  % The true value of streamfunction simulation
ps = zeros(MN, NT+1);  % The noised observation of streamfunction
initial = zeros(MN, 1);  % The initial value for inverse Poisson via SOR method

% From streamfunction, calculate vorticity
% xs(:,1) = laplacian(psiphys,M,N,dx,dy);

% From vorticity, calculate streamfunction
omphys = bdcondition(omphys,M,N);
xst(:,1) = omphys;
psiphys = inversepoisson(xst(:, 1),initial,M,N,dx,dy,optsur);

% Apply boundary conditions on fields
% xs(:,1) = bdcondition(xs(:,1),M,N);
% psiphys = bdcondition(psiphys,M,N);

psiold = psiphys;
pst(:, 1) = psiphys;

% Simulate the model
for i = 1:NT

    xst(:, i+1) = f(xst(:, i),pst(:, i),M,N,dx,dy,dt,AH,BH);

    work = inversepoisson(xst(:, i+1),2*pst(:, i+1)-psiold,M,N,dx,dy,optsur);
    psiold = pst(:, i+1);
    pst(:, i+1) = work;
    
    % Show the vorticity field every 10 time steps
    if mod(i-1,10) == 0
    omphys = reshape(xst(:, i+1),M,N);
    pcolor(real(omphys)');
    caxis([-omegascale omegascale])
    shading('interp')
    axis equal
    axis off
    title 'Vorticity'
    colormap(brcol)
    %colorbar('horiz')
    hold on

    hold off
    pause(0.000000001)
    end

end

% Add noise to each time step
for i = 1:NT+1
    xs(:, i) = xst(:, i) + 0.5*omegascale*random('Normal',0,1,M*N,1);

    psiold = ps(:, i);
    work = inversepoisson(xs(:, i),2*ps(:, i)-psiold,M,N,dx,dy,optsur);
    ps(:, i) = work;

    if mod(i-1,10) == 0
    omphys = reshape(xs(:, i),M,N);
    pcolor(real(omphys)');
    caxis([-omegascale omegascale])
    shading('interp')
    axis equal
    axis off
    title 'Vorticity'
    colormap(brcol)
    %colorbar('horiz')
    hold on

    hold off
    pause(0.000000001)
    end
end

% Recovery via ADMM

tic; % start timer


iter_steps = 5000;  % Number of iterations for ADMM
num = floor(NT / 10);  % Take observations ever num steps (Totally 10 observations)
rho = 1.5;  % The parameter \rho in ADMM
mu = 20;  % The scaling parameter \mu in ADMM
alpha = 0.1;  % The parameter for background error (Here we set the background information equal to the inital observation for simplicity)

x = zeros(M*N, NT+1);  % The voticity
p = zeros(M*N, NT+1);  % The streamfunction
pold = zeros(M*N, NT+1);  % The last step streamfunction serve as the initial value for SOR
xt = zeros(M*N, NT+1);  % Intermediate variable for x
lmda = zeros(M*N, NT+1);  % The array for \lambda (The first element is not used)
lmdat = zeros(M*N, NT+1);  % Intermediate variable for \lambda
fval = zeros(iter_steps+1, 1);  % Store the value of objective function for each iteration
cons_err = zeros(iter_steps+1, 1);  % Store the constraint error for each iteration
energy = zeros(iter_steps+1, 1);  % Store the energy for each iteration

% Calculate the value of objective function for initial value
fval(1) = dx*dy*0.5*alpha*norm(x(:,1)-xs(:,1), 2)^2;
for k = 0:floor(NT/num)
    fval(1) = fval(1)+(dx*dy*0.5/(NT/num))*norm(x(:,k*num+1)-xs(:,k*num+1))^2;
end
% Calculate the energy for initial value
energy(1) = -dx*dy*0.5*alpha*dot(p(:,1)-ps(:,1), x(:,1)-xs(:,1));
for k = 0:floor(NT/num)
    energy(1) = energy(1)-(dx*dy*0.5/(NT/num))*dot(p(:,k*num+1)-ps(:,k*num+1), x(:,k*num+1)-xs(:,k*num+1));
end
% Calculate the constraint error for initial value
for k = 1:NT
    cons_err(1) = cons_err(1) + dx*dy*norm(x(:,k+1) - f(x(:,k),p(:,k),M,N,dx,dy,dt,AH,BH))^2;
end

% Calculate the streamfunction from vorticity
for j = 1:NT+1
    p(:, j) = inversepoisson(x(:,j),initial,M,N,dx,dy,optsur);
end

pold = p;

% Multi-block ADMM with Jacobi decomposition
for i = 1:iter_steps

    xt(:, 1) = subproblem(x(:,1),x(:,1),x(:,2),pold(:,1),pold(:,1),pold(:,2),lmda(:,2),lmda(:,2),xs(:,1),M,N,dx,dy,dt,rho,AH,BH,optsur,(0.5/(NT/M))*mu+0.5*alpha*mu,0,1);
    % Calculate the streamfunction
    p(:, 1) = inversepoisson(xt(:,1), pold(:,1),M,N,dx,dy,optsur);

    for j = 2:NT
        if mod(j-1, num) == 0
            xt(:, j) = subproblem(x(:,j-1),x(:,j),x(:,j+1),pold(:,j-1),pold(:,j),pold(:,j+1),lmda(:,j),lmda(:,j+1),xs(:,j),M,N,dx,dy,dt,rho,AH,BH,optsur,(0.5/(NT/M))*mu,1,1);
        else
            xt(:, j) = subproblem(x(:,j-1),x(:,j),x(:,j+1),pold(:,j-1),pold(:,j),pold(:,j+1),lmda(:,j),lmda(:,j+1),xs(:,j),M,N,dx,dy,dt,rho,AH,BH,optsur,0,1,1);
        end
        % Calculate the streamfunction and new \lambda
        p(:, j) = inversepoisson(xt(:,j), pold(:,j),M,N,dx,dy,optsur);
        lmdat(:, j) = lmda(:,j)-rho*(xt(:,j)-f(xt(:,j-1),p(:,j-1),M,N,dx,dy,dt,AH,BH));
    end

    xt(:, NT + 1) = subproblem(x(:,NT),x(:,NT+1),x(:,NT+1),pold(:,NT),pold(:,NT+1),pold(:,NT+1),lmda(:,NT+1),lmda(:,NT+1),xs(:,NT+1),M,N,dx,dy,dt,rho,AH,BH,optsur,(0.5/(NT/M))*mu,1,0);
    % Calculate the streamfunction and new \lambda
    p(:, NT + 1) = inversepoisson(xt(:,NT+1), pold(:,NT+1),M,N,dx,dy,optsur);
    lmdat(:, NT + 1) = lmda(:,NT + 1)-rho*(xt(:,NT+1)-f(xt(:,NT),p(:,NT),M,N,dx,dy,dt,AH,BH));

    % Update the vorticity, \lambda, and streamfunction
    x = xt;
    lmda = lmdat;
    pold = p;

    % Store the value of objective function
    fval(i + 1) = dx*dy*0.5*alpha*norm(x(:,1)-xs(:,1),2)^2;
    for k = 0:floor(NT/num)
        fval(i + 1) = fval(i+1)+(dx*dy*0.5/(NT/num))*norm(x(:,k*num+1)-xs(:,k*num+1), 2)^2;
    end
    % Store the energy
    energy(i + 1) = -dx*dy*0.5*alpha*dot(p(:,1)-ps(:,1), x(:,1)-xs(:,1));
    for k = 0:floor(NT/num)
        energy(i + 1) = energy(i+1)-(dx*dy*0.5/(NT/num))*dot(p(:,k*num+1)-ps(:,k*num+1), x(:,k*num+1)-xs(:,k*num+1));
    end
    % Store the constraint error
    for k = 1:NT
        cons_err(i + 1) = cons_err(i+1) + dx*dy*norm(x(:,k+1) - f(x(:,k),p(:,k),M,N,dx,dy,dt,AH,BH))^2;
    end

    % Display the result each iteration
    disp(['Iteration: ', num2str(i)]);
    disp(['Error: ', num2str(fval(i+1))]);
    disp(['Energy Error: ', num2str(energy(i+1))]);
    disp(['Constraint Error: ', num2str(cons_err(i+1))]);

    %sectionTime = toc;
    %disp(['Time: ', num2str(sectionTime), ' Seconds']);
end

% Display the running time
elapsedTime = toc; % end timer
disp(['Time: ', num2str(elapsedTime), ' Seconds']);


% Show the recovered vorticity field every 10 time steps
for i = 1:NT+1
    if mod(i-1,10) == 0
    omphys = reshape(x(:, i),M,N);
    pcolor(real(omphys)');
    caxis([-omegascale omegascale])
    shading('interp')
    axis equal
    axis off
    title 'Recovered Vorticity'
    colormap(brcol)
    %colorbar('horiz')
    hold on

    hold off
    pause(0.000000001)
    end
end

% Plot the convergence image for objective function
iterations = 1:iter_steps+1;  % Number of iterations

figure;
plot(iterations, fval, '-o', 'LineWidth', 2, 'MarkerSize', 6); % Plot the number of iterations and the error
xlabel('Iteration', 'FontSize', 12);  % x-axis label
ylabel('fval', 'FontSize', 12);      % y-axis label
title('Function value', 'FontSize', 14); % Title of the graph
grid on; % Display the grid
set(gca, 'YScale', 'log'); % Use a logarithmic scale for clearer display


% Plot the convergence image for constraint error
iterations = 1:iter_steps+1;  % Number of iterations

figure;
plot(iterations, cons_err, '-o', 'LineWidth', 2, 'MarkerSize', 6); % Plot the number of iterations and the error
xlabel('Iteration', 'FontSize', 12);  % x-axis label
ylabel('cons_err', 'FontSize', 12);      % y-axis label
title('Constraint error', 'FontSize', 14); % Title of the graph
grid on; % Display the grid
set(gca, 'YScale', 'log'); % Use a logarithmic scale for clearer display


