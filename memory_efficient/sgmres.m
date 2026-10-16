function [x,res, cond_M] = sgmres(A,f,x1,d,k,tol)
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
B_buffer = zeros(n,k+1); % Store last k basis vectors 
C = zeros(s, d);  % Sketched reduced matrix

b1 = r1/norm(r1);
B_buffer(:,1) = b1;
Ab = A*b1;


% Sketch the first column of M(w1)
c1 = omega.*Ab;
c1= fft(c1);
C(:,1) = c1(idxs)*sqrt(n/s);
clear c1;

for j=2:d
    if mod(j,10) == 0
        fprintf('Iteration %d / %d \n', j, d)
    end

    % Truncated Arnoldi: orthogonalize against last k vectors
    b = Ab;
    start_idx = max(1, j-k);
    
    for pass = 1:2
        for i = start_idx:j-1
                buffer_idx = mod(i-1, k+1) + 1;
                b = b - (B_buffer(:,buffer_idx)' * b) * B_buffer(:,buffer_idx);
        end
    end
    % Normalize
    b_norm = norm(b);
    buffer_idx = mod(j-1, k+1) + 1;
    B_buffer(:,buffer_idx) = b / b_norm;

    % Compute A*B(:,j) for next iteration and sketch immediately
    Ab = A * B_buffer(:,buffer_idx);
    c = exp(1i * theta).*Ab;
    c = fft(c);
    C(:,j) = c(idxs)*sqrt(n/s);
    clear c;
end
clear Ab B_buffer b omega idxs; 

cond_M = cond(C); % cannot check cond M directly, this should be similar

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

%% Regenerating the basis
% Perform k-truncated Arnoldi
B_buffer = zeros(n,k+1); % Store last k basis vectors 

b1 = r1/norm(r1);
B_buffer(:,1) = b1;
Ab = A*b1;

% Regenerate and accumulate solution
x = x1 + b1* y(1);

for j=2:d
    if mod(j, 10) == 0
        fprintf('Regenerating basis vector %d / %d\n', j, d);
    end
    
    % Truncated Arnoldi: orthogonalize against last k vectors
    b = Ab;
    start_idx = max(1, j-k);
    
    for pass = 1:2
        for i = start_idx:j-1
                buffer_idx = mod(i-1, k+1) + 1;
                b = b - (B_buffer(:,buffer_idx)' * b) * B_buffer(:,buffer_idx);
        end
    end

    % Normalize
    b_norm = norm(b);
    buffer_idx = mod(j-1, k+1) + 1;
    B_buffer(:,buffer_idx) = b / b_norm;
    Ab = A*B_buffer(:,buffer_idx);

    % Accumulate solution
    x = x + B_buffer(:,buffer_idx) * y(j);
end

res = norm(f - A*x) / f_norm;
fprintf("Relative Residual: %e \n", res);
fprintf("Finished Computation.\n")
end