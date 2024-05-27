function [R] = channel_ourmodel(RIS_row,RIS_col,TX,Rx_position,frequency)
numda = 3*(10^8)/frequency;
RIS_patch_position = zeros(3,RIS_row*RIS_col);
R = zeros(RIS_row*RIS_col,1);
for k = 1:RIS_row
    for j =1:RIS_col
        RIS_patch_position(1,(k-1)*RIS_col+j) = numda/2*(j-1-(RIS_col-1)/2);
        RIS_patch_position(2,(k-1)*RIS_col+j) = numda/2*(k-1-(RIS_row-1)/2);
    end
end
R0 = zeros(RIS_col*RIS_row,1);
[N,~] = size(TX.d0);
for n = 1:N
    d1 = TX.d0(n,1);
    theta = deg2rad(TX.theta_in(n,1));
    phi = deg2rad(TX.phi_in(n,1));
%    Tx_position = [d1*sin(theta)*cos(phi);d1*sin(theta)*sin(phi);d1*cos(theta)];
    u_Arr = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];

    %%%%%%%%%%%% The theta and phi of Rx %%%%%%%%%%%%
    d2 = sqrt(Rx_position(1,1)^2+Rx_position(2,1)^2+Rx_position(3,1)^2);

    u_Dep = Rx_position./d2;
    A = exp(2*pi*1i*RIS_patch_position'*u_Dep/numda);  %.'conj（A'), ' is conjugate transpose. Here A is the transpose of Bi
    B = exp(2*pi*1i*RIS_patch_position'*u_Arr/numda);
    R0 = 100/(d1*d2)*A.*B;    % For computational convenience, the numerical values are first scaled up, and then the power calculation results are scaled down after obtaining them.
    R = R + R0;
end