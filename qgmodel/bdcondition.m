function field = bdcondition(field, M, N)
% Impose boundary condition on 2d field
% Dirichlet boundary condition with value 0

for j = 1:N
    field(j) = 0;
    field((M-1)*N+j) = 0;
end

for i = 1:M
    field((i-1)*N+1) = 0;
    field((i-1)*N+N) = 0;
end

end