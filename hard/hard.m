clear all; clc;
%% 8.1.3 Symmetric linear system 
tol = 1e-20;
% Create A matrix
%A = laplacian([500,500],{'DD','DD'});
load("fs_680_1.mat");
A = Problem.A;
% Check if the matrix is sparse and how large it is
fprintf('Matrix is sparse: %d\n', issparse(A));
% Check if the matrix is single or double
if isfloat(A)
    if isa(A, 'single')
        fprintf('Matrix is of type: single\n');
    elseif isa(A, 'double')
        fprintf('Matrix is of type: double\n');
    end
else
    fprintf('Matrix is not a floating-point type.\n');
end
% Print the dimension of the problem
fprintf('Dimension of the problem: n = %d \n', size(A, 1));
% Generate traceless vector and compute the rhs
sol  = rand(size(A, 2), 1);
sol = sol - sum(sol) / size(sol, 1);
f0 = A*sol;

% Points
checkpts = [0 10 20 30 40 50 60 70 80 90 100];
%checkpts = 0:1:100;
d = checkpts(end);
%% Solve the linear system by built-in gmres
% TODO: add other frequencies and pcg
% Solve the linear system by gmres with different restart values
restarts = [10];
resvecs = zeros(5,length(checkpts)); % Initialize cell array to store residuals
times = zeros(5,length(checkpts));
for i=1:length(checkpts)
    % gmres
    tic;
    [~,~,res]  = gmres(A, f0, [], tol, checkpts(i));
    resvecs(1,i)=res;
    times(1,i) = toc;
    fprintf("Computed gmres for %d...\n",checkpts(i))
    fprintf("Elapsed time: %.0e s...\n", toc)
    % restarted gmres - GENERALIZED VERSION
    maxit = checkpts(i);
    outer_iter = ceil(maxit / restarts);
    tic;
    [~, ~, res] = gmres(A, f0, restarts, tol, outer_iter);
    resvecs(2,i)=res; % Store the residual vector for each restart value
    times(2,i) = toc;
    fprintf("Computed restarted gmres(%d) for %d iterations...\n", restarts, maxit)
    fprintf("Elapsed time: %.0e s...\n", toc)
end
%% Solve the linear system by sgmres: 2-truncated Arnoldi
x0 = zeros(size(A, 1), 1);  % Zero initial guess for sgmres
stol = 10^3;
conds2 = zeros(length(checkpts), 1);

k = 2;
for i=1:length(checkpts)
    tic;
    [~,res, cond_M] = sgmres(A,f0,x0,checkpts(i),k,stol);
    resvecs(3,i) = res;
    times(3,i) = toc;
    conds2(i) = cond_M;
    fprintf("Computed sgmres for %d...\n",checkpts(i))
    fprintf("Elapsed time: %.0e s...\n", toc)
end
%% Solve the linear system by sgmres: 10-truncated Arnoldi
k = 10;
conds10 = zeros(length(checkpts), 1);
for i=1:length(checkpts)
    tic;
    [~,res, cond_M] = sgmres(A,f0,x0,checkpts(i),k,stol);
    resvecs(4,i) = res;
    times(4,i) = toc;
    conds10(i) = cond_M;
    fprintf("Computed sgmres for %d...\n",checkpts(i))
    fprintf("Elapsed time: %.0e s...\n", toc)
end

%% Solve the linear system by sgmres: 2-truncated Arnoldi whiten
x0 = zeros(size(A, 1), 1);  % Zero initial guess for sgmres
stol = 10^3;
conds_whiten = zeros(length(checkpts), 1);

k = 2;
for i=1:length(checkpts)
    tic;
    [~,res, cond_M] = sgmres_whiten(A,f0,x0,checkpts(i),k,stol);
    resvecs(5,i) = res;
    times(5,i) = toc;
    conds_whiten(i) = cond_M;
    fprintf("Computed sgmres for %d...\n",checkpts(i))
    fprintf("Elapsed time: %.0e s...\n", toc)
end
%% Plot
% Number of sGMRES sample points

figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

%% ---- Residuals ----
nexttile;

semilogy(checkpts, resvecs(1,:), 'b-o', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES'); hold on;
hold on;
semilogy(checkpts, resvecs(2,:), 'k-s', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES(10)');
semilogy(checkpts, resvecs(3,:), 'g-d', 'LineWidth', 1.2, ...
    'DisplayName', 'sGMRES k=2');
semilogy(checkpts, resvecs(4,:), 'm-^', 'LineWidth', 1.2, ...
    'DisplayName', 'sGMRES k=10');
semilogy(checkpts, resvecs(5,:), 'c-+', 'LineWidth', 1.2, ...
    'DisplayName', 'whiten');

xlabel('Iterations','Interpreter','latex');
ylabel('Residual $\|Ax-f\|_2 / \|f\|_2$','Interpreter','latex');
title(['$n = ', num2str(size(A,1)), '$'],'Interpreter','latex');
legend('Location','best');
grid on;
xlim([checkpts(1) checkpts(end)]);
% Set x-axis ticks to match the checkpts values
xticks(checkpts);
hold off;

%% ---- Computational times ----
nexttile;
semilogy(checkpts, conds2, 'g-d', 'LineWidth', 1.2, ...
    'DisplayName', 'sGMRES k=2');
hold on;
semilogy(checkpts, conds10, 'm-^', 'LineWidth', 1.2, ...
    'DisplayName', 'sGMRES k=10');
semilogy(checkpts, conds_whiten, 'c-+', 'LineWidth', 1.2, ...
    'DisplayName', 'whiten');

xlabel('Iterations','Interpreter','latex');
ylabel('$\kappa_2(B)$','Interpreter','latex');
grid on;
legend('Location','best');
xlim([checkpts(1) checkpts(end)]);
hold off;
% Set x-axis ticks to match the checkpts values
xticks(checkpts);

saveas(gcf,'hard.pdf')