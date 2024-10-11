function [A, B] = arakawa_g(o, p, M, N, dx, dy)
% compute the gradient of arakawa method

A = spalloc(M*N, M*N, 9*M*N); % Store the gradient with respect to \psi
B = spalloc(M*N, M*N, 9*M*N); % Stor the gradient with respect to \omega

for i = 2:M-1
    for j = 2:N-1
        % Gradient of jpp
        A((i-1)*N + j, i*N + j) = A((i-1)*N + j, i*N + j) + (o((i-1)*N+j+1)-o((i-1)*N+j-1))/(12*dx*dy);
        A((i-1)*N + j, (i-2)*N + j) = A((i-1)*N + j, (i-2)*N + j) - (o((i-1)*N+j+1)-o((i-1)*N+j-1))/(12*dx*dy);
        A((i-1)*N + j, (i-1)*N + j+1) = A((i-1)*N + j, (i-1)*N + j+1) - (o(i*N+j)-o((i-2)*N+j))/(12*dx*dy);
        A((i-1)*N + j, (i-1)*N + j-1) = A((i-1)*N + j, (i-1)*N + j-1) + (o(i*N+j)-o((i-2)*N+j))/(12*dx*dy);
        
        B((i-1)*N + j, (i-1)*N + j+1) = B((i-1)*N + j, (i-1)*N + j+1) + (p(i*N+j)-p((i-2)*N+j))/(12*dx*dy);
        B((i-1)*N + j, (i-1)*N + j-1) = B((i-1)*N + j, (i-1)*N + j-1) - (p(i*N+j)-p((i-2)*N+j))/(12*dx*dy);
        B((i-1)*N + j, i*N + j) = B((i-1)*N + j, i*N + j) - (p((i-1)*N+j+1)-p((i-1)*N+j-1))/(12*dx*dy);
        B((i-1)*N + j, (i-2)*N + j) = B((i-1)*N + j, (i-2)*N + j) + (p((i-1)*N+j+1)-p((i-1)*N+j-1))/(12*dx*dy);
        
        % Gradient of jpx
        A((i-1)*N + j, i*N + j) = A((i-1)*N + j, i*N + j) + (o(i*N+j+1)-o(i*N+j-1))/(12*dx*dy);
        A((i-1)*N + j, (i-2)*N + j) = A((i-1)*N + j, (i-2)*N + j) - (o((i-2)*N+j+1)-o((i-2)*N+j-1))/(12*dx*dy);
        A((i-1)*N + j, (i-1)*N + j+1) = A((i-1)*N + j, (i-1)*N + j+1) - (o(i*N+j+1)-o((i-2)*N+j+1))/(12*dx*dy);
        A((i-1)*N + j, (i-1)*N + j-1) = A((i-1)*N + j, (i-1)*N + j-1) + (o(i*N+j-1)-o((i-2)*N+j-1))/(12*dx*dy);
        
        B((i-1)*N + j, i*N + j+1) = B((i-1)*N + j, i*N + j+1) + (p(i*N+j)-p((i-1)*N+j+1))/(12*dx*dy);
        B((i-1)*N + j, i*N + j-1) = B((i-1)*N + j, i*N + j-1) + (-p(i*N+j)+p((i-1)*N+j-1))/(12*dx*dy);
        B((i-1)*N + j, (i-2)*N + j+1) = B((i-1)*N + j, (i-2)*N + j+1) + (-p((i-2)*N+j)+p((i-1)*N+j+1))/(12*dx*dy);
        B((i-1)*N + j, (i-2)*N + j-1) = B((i-1)*N + j, (i-2)*N + j-1) + (p((i-2)*N+j)-p((i-1)*N+j-1))/(12*dx*dy);

        % Gradient of jxp
        A((i-1)*N + j, (i-2)*N + j-1) = A((i-1)*N + j, (i-2)*N + j-1) + (-o((i-2)*N+j)+o((i-1)*N+j-1))/(12*dx*dy);
        A((i-1)*N + j, (i-2)*N + j+1) = A((i-1)*N + j, (i-2)*N + j+1) + (o((i-2)*N+j)-o((i-1)*N+j+1))/(12*dx*dy);
        A((i-1)*N + j, i*N + j+1) = A((i-1)*N + j, i*N + j+1) + (-o(i*N+j)+o((i-1)*N+j+1))/(12*dx*dy);
        A((i-1)*N + j, i*N + j-1) = A((i-1)*N + j, i*N + j-1) + (o(i*N+j)-o((i-1)*N+j-1))/(12*dx*dy);
        
        B((i-1)*N + j, (i-1)*N + j+1) = B((i-1)*N + j, (i-1)*N + j+1) + (p(i*N+j+1)-p((i-2)*N+j+1))/(12*dx*dy);
        B((i-1)*N + j, (i-1)*N + j-1) = B((i-1)*N + j, (i-1)*N + j-1) - (p(i*N+j-1)-p((i-2)*N+j-1))/(12*dx*dy);
        B((i-1)*N + j, i*N + j) = B((i-1)*N + j, i*N + j) - (p(i*N+j+1)-p(i*N+j-1))/(12*dx*dy);
        B((i-1)*N + j, (i-2)*N + j) = B((i-1)*N + j, (i-2)*N + j) + (p((i-2)*N+j+1)-p((i-2)*N+j-1))/(12*dx*dy);
    end
end



end