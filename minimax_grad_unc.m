function [x, opt, Fx, J, iter] = minimax_grad_unc(x0, lambd, max_iter, Coeff, Lip)
% solve max-min problem by minimax formulation

% initialize
lb = Coeff.lb;
ub = Coeff.ub;
x = x0;
% x_hist = zeros(n, max_iter);
Fx = zeros(max_iter, 1);
% J = zeros(max_iter, 1);
% load data
% lambd = opt.lambd;
tol = 1e-4;
% %%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on greedy FISTA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1;
z = x;
z1 = x;
y = x;
for ii=1:max_iter
    x = z1;
%     x_hist(:, ii) = x;
    % evaluate constraints
    [opt] = compute_function_derivative(x,lambd, Coeff);
    penalty_x = opt.Fx;
%     grad_Fx = opt.dFx;
    J = opt.J;
    Fx(ii) = penalty_x;
    if ii == 1
        flag = 1;
%     elseif penalty_x > penalty_y + dot(grad_Fy, x - y) + Lip/2*norm(x-y)^2
    elseif penalty_x > penalty_y  - Lip/2*norm(x-y)^2
        Lip = Lip *2;
        flag = 0;
    else
        flag = 1;
        if dot_ <  tol*norm_
            t = 1;
        else
            t = 0;
        end
    end
    if flag
        % Nesterov's optimal gradient steps
%         t1 = (1+sqrt(1+r*t^2))/2;
        y = z1 + (1-t)*(z1-z);
%         y = z1 + (t-1)/t1*(z1 - z);
        z = z1;
%         t = t1;       
        % evaluate constraints and gradients
        [opt] = compute_function_derivative(y,lambd, Coeff);
        grad_Fy = opt.dFx;
        penalty_y = opt.Fx;
    end
    % Nesterov's optimal gradient steps
    z1 = y - grad_Fy/Lip; 
    norm_ = norm(z1 - y)*norm(z1 - x);
    dot_ = dot(y - z1, x - z1);
%     dJ = opt.dJ;
%     w = opt.w;
%     dJ = dJ(w>eps);
    opt.optimal = norm(grad_Fy, 'inf');
    if ii==1
        continue;
    elseif opt.optimal<max(tol, 1e-6*Lip) %&& abs(Fx(ii)-Fx(ii-1))<tol/lambd %&& min(grad_Fx.*(x-lb<eps))>-tol && max(grad_Fx.*(ub-x<eps))<tol%
        break;
    end
end
iter = ii;
opt.Lip = Lip;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fx = Fx(1:ii);
% J = J(1:ii);
% x_hist = x_hist(:, 1:ii);
% [~, indmin] = min(Fx);
% Fx = Fx(1:indmin);
% x_hist = x_hist(:, 1:indmin);
% x = x_hist(:, end);
% figure; plot(Fx);