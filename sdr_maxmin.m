function w_sdr = sdr_maxmin(N,Q,count,K)
f_tmp = 0;
w_tmp = zeros(N,1);

for k=1:count
    r = (randn(N,1)+1i*randn(N,1)).*sqrt(1/2);   % (N,1)
    %     cvx_begin
    %     variable X(N+1,N+1) symmetric semidefinite   %The variable is an (N)*(N) symmetric symmetric positive semi-definite matrix.
    %     for i = 1:K
    %     maximize ( minimize (    real(trace(Q{i}*X))      )         )
    %     subject to
    %     diag(V) == 1;
    %     end
    %     cvx_end
    cvx_begin
    variable X(N,N) symmetric semidefinite
    variable t
    maximize(t)
    subject to
    diag(X) == 1;
    for i = 1:K
        t <= real(trace(Q{i}*X));
    end
    cvx_end
    [U,Sigma] = eig(X);
 %%%%     recover 1, Gaussian randomization
    w = U*Sigma^(1/2)*r;   % (N*1)
    f = 10e10;  % f take values of inf
    for i = 1:K
        f = min(f,w'*Q{i}*w);
    end
    if f>f_tmp
        f_tmp = max(f,f_tmp);
        w_tmp = w;
    end
end
% %%%%% recover 2, the eigenvector corresponding the max eigenvalue
% [~, index] = max(diag(Sigma));
% w = U(:, index);
% f = 10e10;
% for i = 1:K
%     f = min(f,w'*Q{i}*w);
% end
% if f>f_tmp
%     f_tmp = max(f,f_tmp);
%     w_tmp = w;
% end
% end
% w_tmp = w_tmp / norm(w_tmp);

theta_opt = angle(w_tmp);
w_sdr = exp(1i*theta_opt);
end