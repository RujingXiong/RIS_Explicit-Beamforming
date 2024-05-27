function [opt] = compute_function_derivative(omega, lambd, Coeff)
h = Coeff.h;
N = Coeff.N;
D = Coeff.D;
K = numel(h);
J = zeros(K,1);
dJ = zeros(N, K);
for k=1:K
    [J(k), dJ(:,k)] = Quadratic(h{k}, omega);
end
F = -2*lambd*J;
w = skew_proj(F, D);
if sum(w)==0
    [~,I] = max(F);
    w(I) = 1;
end
Fx = 1/(4*lambd)*(dot(2*F-w,w)+1);

dFx = -dJ*w;
opt.Fx = Fx; % smoothed function value
opt.minJ = min(J); % exact function value
opt.dFx = dFx; % gradient
opt.lambd = lambd; % 
opt.J = J; 
opt.dJ = dJ; % jacobian
opt.w = w;