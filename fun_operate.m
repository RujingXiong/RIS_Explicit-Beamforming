function [y] = fun_operate(R,phi,K)
y = [];
for k = 1:K
    y = [y,am(R(:,k),phi)];
end
