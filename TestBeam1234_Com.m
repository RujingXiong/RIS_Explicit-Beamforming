%% Author:  Rujing Xiong, Ke Yin, Jialong Lu
 % CreatTime: 2023.11.27
 % Complete:
 % Modified: 
 % E-mail: Rujing@hust.edu.cn
 % Description: Beam comparison of different methods, including Fminimax,
 % QuantRand, and SDR-SDP
 % 
 % 
clc
clear all 
tic
TX.d0 = [5;6;8;4];
TX.theta_in = [0;40;40;60];
TX.phi_in = [0;90;180;0];


theta_out = [0;30;50];
phi_out = [0;0;180];
d_out = [30;30;30];
race = [1;1;1];

K = 3;
RIS_row = 1;
RIS_col = 32;

d1 = 30;


Pt = db2pow(40);
frequency = 3.4e9;
area = (3e8/frequency/2)^2; % the area of a RIS unit 

N = RIS_row * RIS_col;

R = zeros(N,K);
for k=1:K
    R(:,k) = Generate_h(theta_out(k), phi_out(k), d_out(k), RIS_row, RIS_col, TX, race(k), frequency);
end



Coeff.K = K;
Coeff.D = ones(K,1);
Coeff.N = N;
lb = zeros(N, 1);
ub = 2*pi*ones(N, 1);
Coeff.lb = lb;
Coeff.ub = ub;
omega0 = lb + (ub-lb).*rand(N, 1);
%omega0 = rand(1,1)*2*pi.*ones(N, 1);
h = cell(K,1);
for k=1:K
    h{k} = Generate_h(theta_out(k), phi_out(k), d_out(k), RIS_row, RIS_col, TX, race(k), frequency);
end
Coeff.h = h;
[result, ~, Coeff] = TestMaxMin1(Coeff, omega0);
opt = exp(1i.*result.omega);

% opt = discrete(opt,1,RIS_col*RIS_row);
disp(result.J)

%% Fminmax
fun = @(phi)fun_operate(R,phi,K);
lb = zeros(1,N);
ub = 2*pi*ones(1,N);
phi0 = rand(1,1)*2*pi.*ones(1,N);
%phi0 = 2*pi.*ones(1,N);
A = []; % No linear constraints
b = [];
Aeq = [];
beq = [];
nonlcon = [];
options = optimoptions('fminimax','UseParallel',false);
%options = optimoptions('fminimax','UseParallel',10000);
%options.SpecifyObjectiveGradient=true;
options.MaxFunctionEvaluations = 10000000;
options.MaxIterations = 500;
[phi_opt,fval] = fminimax(fun,phi0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%[phi_opt,fval] = fminimax(fun,phi0,A,b,Aeq,beq,lb,ub,nonlcon);
opt_F = exp(1i.*phi_opt).';


%% QuantRand
[opt_Q,J] = QuantRand(6,h);


%% SDR
count = 1;
q = zeros(K,N);
Q = cell(K);
for i = 1:K
q(i,:) = cell2mat(h(i));
Q{i} = q(i,:)'*q(i,:) ;
end
opt_S = sdr_maxmin(N,Q,count,K);



%% Plot
%%
strength1 = area*cal(opt,RIS_row,RIS_col,TX,frequency,d1,Pt);
strength2 = area*cal(opt_F,RIS_row,RIS_col,TX,frequency,d1,Pt);
strength3 = area*cal(opt_Q,RIS_row,RIS_col,TX,frequency,d1,Pt);
strength4 = area*cal(opt_S,RIS_row,RIS_col,TX,frequency,d1,Pt);
x = -90:0.1:90;
figure
plot(x,strength1);
hold on;
plot(x,strength2);
hold on
plot(x,strength3);
hold on 
plot(x,strength4)
xlabel('Angle({\circ})')
ylabel('Power (dBm)')
grid on
legend('Proposed MA','Fminimax','QuantRand','SDR-SDP')
box on 

toc

