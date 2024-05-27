function y = skew_proj(x, d)
%SKEW_PROJ find the projection of a point onto the simplex
% y = skew_proj(x, d) project x onto the simplex
% {x: d'*x=1, x>=0} where d>0
n = length(x);
N = ones(n,1);
% I = [];
% J = N;
% y = x;
% while ~isempty(J)
%     y(J) = y(J) - (d(J)'*y(J) - 1)/norm(d(J))^2*d(J);
%     I_t = find(y<0);
%     if isempty(I_t)
%         break;
%     else
%         I = union(I, I_t);
%         J = setdiff(N, I);
%         y = max(y, 0);
%     end
% end
I = N<0;
J = N>0;
y = x;
while sum(J)
    y(J) = y(J) - (d(J)'*y(J) - 1)/norm(d(J))^2*d(J);
    I_t = y<0;
    if ~(I_t)
        break;
    else
        I = I|I_t;
        J = ~I;
        y = max(y, 0);
    end
end