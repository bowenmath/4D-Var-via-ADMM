function lap = laplacian(field, M, N, dx, dy)
% Compute the laplacian of a 2d vector field

lap = zeros(M*N, 1);

for i = 2:M-1
    for j = 2:N-1
        lap((i-1)*N+j) = (field(i*N+j) - 2*field((i-1)*N+j) + field((i-2)*N+j))/(dx^2) + ...
            (field((i-1)*N+j+1) - 2*field((i-1)*N+j) + field((i-1)*N+j-1))/(dy^2);
    end
end

end