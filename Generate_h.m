function h = Generate_h(theta, phi, d, RIS_row, RIS_col, TX, race, frequency)
theta_out = deg2rad(theta); 
phi_out = deg2rad(phi);
Rx_position1 = [d*sin(theta_out)*cos(phi_out);d*sin(theta_out)*sin(phi_out);d*cos(theta_out)];
h = sqrt(1/race).*(channel_ourmodel(RIS_row,RIS_col,TX,Rx_position1,frequency));
end