clear all; clc;
%% 8.1.3 Symmetric linear system 
tol = 1e-20;
% Create A matrix
A = laplacian([1000,1000],{'DD','DD'});
%A = laplacian([80,80],{'DD','DD'});
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
sol  = rand(size(A, 1), 1);
sol = sol - sum(sol) / size(sol, 1);
f0 = A*sol;

% Points 
checkpts = [100,200,500,1000,2000,3000];
%checkpts = [100,3000];
d = checkpts(end);
%% Solve the linear system by built-in gmres
% TODO: add other frequencies and pcg
% Solve the linear system by gmres with different restart values
restarts = [10, 30, 100];
resvecs = zeros(length(restarts)+3,length(checkpts)); % Initialize cell array to store residuals
times = zeros(length(restarts)+3,length(checkpts));
for i=1:length(checkpts)
    % pcg
    tic;
    [~,~,res]  = pcg(A, f0, tol, checkpts(i));
    resvecs(1,i)=res;
    times(1,i) = toc;
    fprintf("Computed pcg for %d...\n",checkpts(i))
    fprintf("Elapsed time: %.0e s...\n", toc)
    % gmres
    tic;
    [~,~,res]  = gmres(A, f0, [], tol, checkpts(i));
    resvecs(2,i)=res;
    times(2,i) = toc;
    fprintf("Computed gmres for %d...\n",checkpts(i))
    fprintf("Elapsed time: %.0e s...\n", toc)
    % restarted gmres - GENERALIZED VERSION
    for j = 3:length(restarts)+2
        n_restarts = restarts(j-2);
        maxit = checkpts(i);
        outer_iter = ceil(maxit / n_restarts);
        
        tic;
        [~, ~, res] = gmres(A, f0, n_restarts, tol, outer_iter);
        resvecs(j,i)=res; % Store the residual vector for each restart value
        times(j,i) = toc;
        fprintf("Computed restarted gmres(%d) for %d iterations...\n", n_restarts, maxit)
        fprintf("Elapsed time: %.0e s...\n", toc)
    end
end
%% Solve the linear system by sgmres: 2-truncated Arnoldi
k = 2;
x0 = zeros(size(A, 1), 1);  % Zero initial guess for sgmres
stol = 1e15;
conds = zeros(length(checkpts), 1);
for i=1:length(checkpts)
    tic;
    [~,res, cond_M] = sgmres(A,f0,x0,checkpts(i),k,stol);
    resvecs(6,i) = res;
    times(6,i) = toc;
    conds(i) = cond_M;
    fprintf("Computed sgmres for %d...\n",checkpts(i))
    fprintf("Elapsed time: %.0e s...\n", toc)
end

%% Plot
% Number of sGMRES sample points

figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

%% ---- Residuals ----
nexttile;
loglog(checkpts, resvecs(1,:), 'b-o', 'LineWidth', 1.2, ...
    'DisplayName', 'CG'); 
hold on;
loglog(checkpts, resvecs(2,:), 'k-s', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES');
loglog(checkpts, resvecs(3,:), 'g-d', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES(10)');
loglog(checkpts, resvecs(4,:), 'm-^', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES(30)');
loglog(checkpts, resvecs(5,:), 'c-+', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES(100)');
loglog(checkpts, resvecs(6,:), 'r-*', 'LineWidth', 1.2, ...
    'DisplayName', 'sGMRES');
% Add condition number annotations for sGMRES
x = checkpts;
y = resvecs(6,:);

for j = 1:length(checkpts)
    txt = sprintf('$\\kappa_2(AB)=%.0e$', conds(j));
    text(x(j), y(j), txt, ...
        'Interpreter','latex', ...
        'FontSize', 10, ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment','center');
end

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
loglog(checkpts, times(1,:), 'b-o', 'LineWidth', 1.2, ...
    'DisplayName', 'CG'); 
hold on;
loglog(checkpts, times(2,:), 'k-s', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES');
loglog(checkpts, times(3,:), 'g-d', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES(10)');
loglog(checkpts, times(4,:), 'm-^', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES(30)');
loglog(checkpts, times(5,:), 'c-+', 'LineWidth', 1.2, ...
    'DisplayName', 'GMRES(100)');
loglog(checkpts, times(6,:), 'r-*', 'LineWidth', 1.2, ...
    'DisplayName', 'sGMRES');

xlabel('Iterations','Interpreter','latex');
ylabel('Time (s)','Interpreter','latex');
title('Computational time','Interpreter','latex');
grid on;
legend('Location','best');
xlim([checkpts(1) checkpts(end)]);
hold off;
% Set x-axis ticks to match the checkpts values
xticks(checkpts);

saveas(gcf,'symmetric.pdf')