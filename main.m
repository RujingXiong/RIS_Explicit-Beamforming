%% Author:  Rujing Xiong, Ke Yin, Jialong Lu
 % CreatTime: 2023.08.27
 % Complete:
 % Modified: 
 % E-mail: Rujing@hust.edu.cn
 % Description: Fair beam allocation through RIS
 % 
 % 

clc
clear all 
tic

TX.d0 = [5;6;8;4];
TX.theta_in = [0;40;40;60];
TX.phi_in = [0;90;180;0];


theta_out = [0;20; 30; 45; 40];
phi_out =   [0;180; 0; 0; 180];
d_out = [30;30;30;30;30];

race = [4;4;2;2;1];

K = 5;
RIS_row = 32;
RIS_col = 32;

d1 = 30;

Pt = db2pow(40);
frequency = 3.4e9;
area = (3e8/frequency/2)^2; % the area of a RIS unit 

N = RIS_row * RIS_col;
Coeff.K = K;
Coeff.D = ones(K,1);
Coeff.N = N;
lb = zeros(N, 1);
ub = 2*pi*ones(N, 1);
Coeff.lb = lb;
Coeff.ub = ub;
% omega0 = lb + (ub-lb).*rand(N, 1);
omega0 = rand(1,1)*2*pi.*ones(N, 1);
h = cell(K,1);
for k=1:K
    h{k} = Generate_h(theta_out(k), phi_out(k), d_out(k), RIS_row, RIS_col, TX, race(k), frequency);
end
Coeff.h = h;
[result, ~, Coeff] = TestMaxMin1(Coeff, omega0);
opt = exp(1i.*result.omega);
% opt = discrete(opt,1,RIS_col*RIS_row);
disp(result.J)



%% plot 3D radiation pattern
% %opt = exp(1i*2*pi.*rand(N,1));
d = 5;
theta = 0:1:90;
phi = -180:1:180;
[theta, phi] = meshgrid(theta, phi);
xMesh = d.*sin(deg2rad(theta)).*cos(deg2rad(phi));
yMesh = d.*sin(deg2rad(theta)).*sin(deg2rad(phi));
zMesh = d.*cos(deg2rad(theta));
[I,J] = size(xMesh);
Power = zeros(I,J);
for i = 1:1:I
    for j = 1:1:J
        disp((i-1)*J+j)
        R0 = channel_ourmodel(RIS_row,RIS_col,TX,[xMesh(i,j);yMesh(i,j);zMesh(i,j)],frequency);
        Power(i,j) = abs(R0' * opt)^2;
    end
end
Power = Power./max(max(Power));

figure
mesh(Power.*xMesh, Power.*yMesh, Power.*zMesh, Power)
xlabel('X')
ylabel('Y')
zlabel('Z')
axis vis3d normal equal
box on
colorbar
figure
mesh(xMesh./d1, yMesh./d1, Power)
xlabel('$$\sin(\theta)\cos(\phi)$$','Interpreter','latex','rotation',15);
ylabel('$$\sin(\theta)\sin(\phi)$$','Interpreter','latex','rotation',-18);
legend('UE\#1, 15m','UE\#2, 12m','UE\#3, 13m','UE\#4, 15m','UE\#5','Interpreter','latex')
zlabel('Power (dB)','Interpreter','latex')
axis vis3d normal
box on
colorbar
%%
strength = area*cal(opt,RIS_row,RIS_col,TX,frequency,d1,Pt);
x = -90:0.1:90;
figure
plot(x,strength,'b-','LineWidth',1)
xlabel('Angle({\circ})')
ylabel('Power (mW)')
grid on
%% poloar db TU
angle = (-90:0.1:90);
C1 = area*cal(opt,RIS_row,RIS_col,TX,frequency,d1,Pt);
C1 = C1./max(C1);
figure('NumberTitle','off','Name','Figure 11.9 in Polar','Position',[0 0 600 850]);
polardb(angle*pi/180,10*log10(abs(C1)),-30,'b');
% hold on;
% angle_of_incident = 0; % mark the incident angle
% plot([angle_of_incident, angle_of_interest], [0, 3], 'r--','LineWidth',2); % use a red dashed line to highlight the incident angle
% hold off;
%title('Beampattern');
grid on;

