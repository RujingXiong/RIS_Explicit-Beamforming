function [result, opt, Coeff] = TestMaxMin1(Coeff, omega0)
% [Coeff, omega0] = GenerateData;
lambd = 1e-4;
max_iter = 1e3;
% tol = 1e-2;
Lip = .1;
K = Coeff.K;
Coeff.D = ones(K,1);
totalIter = 0;
for ii=1:100
    [omega, opt, Fx, J, iter] = minimax_grad_unc(omega0, lambd, max_iter, Coeff, Lip);
    w = opt.w;
    Iw = w>eps;
    gap = max(J(Iw))-min(J(Iw));
    if gap>0; lambd = max(2*lambd, 1/2/gap); end
%     mu = 2*lambd*mean(J(Iw));
%     c = dot(mu-w, w);
%     dd = Coeff.D;
%     dd = ones(K,1);
%     dd(Iw) = (mu-w(Iw))/c;
%     if min(dd)>0
%         Coeff.D = dd;
%     end
    omega0 = omega;
    totalIter = totalIter + iter;
    Lip = opt.Lip;
    if gap < 1e-2
        break;
    end
end
result.omega = omega;
result.Fx = Fx;
result.J = J;
result.Iter = totalIter;
result.gap = gap;