function [x,res, cond_B] = sgmres_whiten(A,f,x1,d,k,tol)
% Input description
% A: Matrix A of size n x n (complex)
% f: Right-hand side vector f of size n (complex)
% x1: Initial guess vector x0 of size n (complex)
% d: Basis dimension for the truncated Arnoldi process
% k: Number of vectors for truncated orthogonalization
% tol: For stability check
    
% Draw subspace embedding by SRFT
s = 2*(d+1);
n = size(A,1);

% Steinhaus random variables
theta = 2*pi*rand(n,1);
omega = exp(1i * theta);

% Subset of random coordinates
idxs = randperm(n, s);

% Apply S to initial residual
r1 = f - A*x1;
f_norm = norm(f);
g = omega.*r1;
g = fft(g);
g = g(idxs);
g = g*sqrt(n/s);

% Perform k-truncated Arnoldi
B = zeros(n,d); % Store ALL basis vectors for simplicity
M = zeros(n,d); % Reduced matrix columns

b1 = r1/norm(r1);
B(:,1) = b1;
M(:,1) = A*b1;
clear r1; clear b1;

for j=2:d
    fprintf('Iteration %d \n', j)
    
    % Truncated Arnoldi: orthogonalize against last k vectors
    w = M(:,j-1);
    start_idx = max(1, j-k);
    for i = j-1:-1:start_idx
        w = w - (B(:,i)' * w) * B(:,i);
    end
    % Normalize
    B(:,j) = w/norm(w);
    M(:,j) = A*B(:,j);

    % Basis whitening
    % Sketch B
    SB = omega.*B(:,1:j);
    SB = fft(SB);
    SB = SB(idxs,:);
    SB = SB*sqrt(n/s);
    fprintf("Condition number of SB: %e\n", cond(SB));
    if cond(SB) > tol
        % whiten the basis
        fprintf("Basis whitening...\n")
        % Sketch Reduced Matrix
        C = omega.*M(:,1:j);
        C = fft(C);
        C = C(idxs,:);
        C = C*sqrt(n/s);
        % QR factorization
        [~, T] = qr(C,"econ");
        B(:,1:j) = B(:,1:j) / T;
        M(:,1:j) = M(:,1:j) / T;
    end
end
clear w;

% Sketch Reduced Matrix
C = omega.*M;
C = fft(C);
C = C(idxs,:);
C = C*sqrt(n/s);
clear omega theta idxs;

% QR factorization
fprintf("QR factorization...\n")
[U, T] = qr(C, "econ");
clear C;

if cond(T) > tol
    fprintf("Warning: condition number of T is %e\n", cond(T))
end

% Solve least squares
fprintf("Solving least squares...\n")
y = T \ (U'*g);
clear U T g;
x = x1 + B*y;
cond_B = cond(B);

res = norm(f - A*x) / f_norm;
fprintf("Relative Residual: %e \n", res);
fprintf("Finished Computation.\n")
end