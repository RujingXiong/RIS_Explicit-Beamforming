function [strength] = cal(w,RIS_row,RIS_col,TX,frequency,d1,Pt)
lambda = 3e8/frequency;
xita = -90:0.1:90;
xita_ra = xita/180*pi;
L = length(xita);
strength = zeros(1,L);
for i = 1:L
    R = lambda^2/(16*pi^2)/100*channel_ourmodel(RIS_row,RIS_col,TX,[d1*sin(xita_ra(i));0;d1*cos(xita_ra(i))],frequency);
%     strength(1,i) = 10*log10(abs(R' * w));
    strength(1,i) = Pt*abs(R' * w)^2;
end