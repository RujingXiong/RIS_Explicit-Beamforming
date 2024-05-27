function [energy] = am(R,phi)
w = (exp(1i*phi)).';
energy = -(abs(R'*w)^2);