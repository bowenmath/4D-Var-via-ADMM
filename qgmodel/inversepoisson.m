function psi = inversepoisson(x, psistart, M, N, dx, dy, optsur)
% Use SOR to compute the inverse Poisson operator, i.e., compute streamfunction from vorticity

psi = psistart;  % Initial value

nmax = M * N;  % Maximum iteration
tol = 1e-6;  % Tolerance of the relative error

n = 0;  % Count the iteration
relerr = 1;  % Relative error
err = zeros(M*N, 1);

% The SOR implementation
while n < nmax && relerr > tol
    n = n+1;
    errtot = 0;
    psitot = 0;
    avrpsi = 0;

    for i = 2:M-1
        for j = 2:N-1
            dx2 = (psi(i*N+j) - 2*psi((i-1)*N+j) + psi((i-2)*N+j))/(dx^2);
            dy2 = (psi((i-1)*N+j+1) - 2*psi((i-1)*N+j) + psi((i-1)*N+j-1))/(dy^2);
            err((i-1)*N+j) = (x((i-1)*N + j) - (dx2+dy2))/(2/(dx^2) + 2/(dy^2));

            psi((i-1)*N+j) = psi((i-1)*N+j) - optsur*err((i-1)*N+j);
            errtot = errtot + err((i-1)*N+j)^2;
            psitot = psitot + psi((i-1)*N+j)^2;
            avrpsi = avrpsi + psi((i-1)*N+j);
        end
    end

    errtot = errtot/(M-2)/(N-2);
    avrpsi = avrpsi/(M-2)/(N-2);

    psitot = psitot/(M-2)/(N-2) - avrpsi^2;
    % For stopping criterium
    if errtot == 0
        relerr = 0;
    else
        relerr = sqrt(errtot/psitot);
    end
    
    % Impose boundary condtion on psif
    psi = bdcondition(psi,M,N);

end
% disp(relerr);

end

