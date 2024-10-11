function jac = arakawa(x, p, M, N, dx, dy)
% Compute the Jacobian using arakawa method

jac = zeros(M*N, 1);

for i = 2:M-1
    for j = 2:N-1
        jpp = (p(i*N+j) - p((i-2)*N+j))*(x((i-1)*N+j+1) - x((i-1)*N+j-1))...
            - (x(i*N+j) - x((i-2)*N+j))*(p((i-1)*N+j+1) - p((i-1)*N+j-1));
        jpx = p((i*N+j))*(x(i*N+j+1) - x(i*N+j-1)) - p((i-2)*N+j)*(x((i-2)*N+j+1) - x((i-2)*N+j-1))...
            - p((i-1)*N+j+1)*(x(i*N+j+1) - x((i-2)*N+j+1)) + p((i-1)*N+j-1)*(x(i*N+j-1) - x((i-2)*N+j-1));
        jxp = - x(i*N+j)*(p(i*N+j+1) - p(i*N+j-1)) + x((i-2)*N+j)*(p((i-2)*N+j+1) - p((i-2)*N+j-1))...
            + x((i-1)*N+j+1)*(p(i*N+j+1) - p((i-2)*N+j+1)) - x((i-1)*N+j-1)*(p(i*N+j-1) - p((i-2)*N+j-1));
        jac((i-1)*N+j) = jpp+jpx+jxp;
    end
end
jac = jac/(12*dx*dy);

end