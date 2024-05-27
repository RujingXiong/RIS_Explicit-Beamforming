function[opt,J] = QuantRand(B,h)
[K,~] = size(h);
[N,~] = size(h{1});
opt0 = exp(randi([0,2^B-1],N,1).*(1i*2*pi/(2^B)));
max_iter = 1e5;
iter = 1;
R = [];
for k = 1:K
    R = [R;h{k}'];
end
min_h = min(abs(R*opt0).^2);
while iter < max_iter
    n = randi(N);
    for i = 1:2^B
        opt0(n,1) = opt0(n,1)*exp(1i*2*pi/(2^B));
        if min(abs(R*opt0).^2) > min_h
            min_h = min(abs(R*opt0).^2);
            break
        end
    end
    iter = iter+1;
end
opt = opt0;
J = abs(R*opt0).^2;