% SBAIN_2TANKS COBEM
% Simulink PSO-PID twotank
% Controle via PID para o modelo SISO Não-Linear
% de um Sistema de 2 Tanques Diferencialmente Plano. 

% Nome: STORMS Group
% Data: Junho/2023

% Pré-Comandos
  close all; clc; clear all;

% Parâmetros do Sistema Não-Linear   
 
  % Tanque 01
    D1 = (5+3/8)*2.54/100;          % m
    Area1 = pi*D1^2/4;              % m^2
    H1max = 1;                      % m          30*2.54/100
    Dout = 0.05;                    % m
    Aout1 = pi*(Dout^2)/4;          % m^2
  
  % Tanque 02
    R2max = 0.176;                  % m
    H2max = 1;                      % m
    Aout2 = Aout1*3/4;              % m^2
  
  % Outros
    g = 9.81;                       % m/(s^2)
    pho = 998.23;                   % kg/(m^3)
    Fnom = 60;                      % Hz
    Qnom = 16.67;                   % kg/s  -> 60 m3/h
    
  % Agregados
    num1 = Qnom; den1 = pho*Area1*(Fnom); a1 = num1/den1;
    num2 = Aout1*sqrt(2*g); den2 = Area1; a2 = num2/den2;
    num3 = Aout1*sqrt(2*g)*H2max^2; den3 = pi*R2max^2; a3 = num3/den3;
    num4 = Aout2*sqrt(2*g)*H2max^2; den4 = pi*R2max^2; a4 = num4/den4;
    
% Condições Iniciais do Sistema
  h10 = 0.10; h20 = 0.10;
      
% Tempo de Simulação
  t0 = 0; dt = 0.01; tfinal = 25*4;
  t = t0:dt:tfinal; st = size(t,2);
  
% Planejamento de Trajetória - Saída Plana (Z)

  % Trajetória Constante 
    % Zd = 0.80 + 0*t; d1Zd = 0*t; d2Zd = 0*t;
    
  % Trajetória Constante com 2 valores (STEP)
    fP1 = 0.30; fP2 = 0.60;
    Zd = 0.30 + 0*t; Zd((round(fP1*st)):end) = 0.75; Zd((round(fP2*st)):end) = 0.50;
    d1Zd = 0*t; d2Zd = 0*t;
    
% Planejamento de Trajetória - Variáveis do Sistema
  H2d = Zd;
  
  H1d = (Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2;
     
  Fd = ((2*Zd.^(5/2)*a4 + 2*Zd.^4.*d1Zd)/(a1*a3^2)).*d2Zd + (a4^2*d1Zd + 4*Zd.^3.*d1Zd.^3 + 5*Zd.^(3/2)*a4.*d1Zd.^2 + a2*a3^2*((Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2).^(1/2))/(a1*a3^2);

% Parâmetros do Controlador (PID)
  kP = 150;     % [50:250]
  kI = 30;      % [10:50]
  kD = 5;       % [1:10]
  
  % Saturação do Atuador
    satF_sup = 90;
    satF_inf = 10;
  
% Executando a simulação e apresentando os resultados
  sim('SIMULINK_PSO.slx'); 
  plotResultados_NL;




%% Tunning of PID controller using Particle Swarm Optimization 
%
%
% Author: Wael Mansour (wael192@yahoo.com)
%
% MSc Student, Electrical Enginering Dept, 
% Faculty of Engineering Cairo University, Egypt
%


%% Initialization


% function [Melhor, Dados_Internos, Dados_Externos, Melhor_sol] = Trab_final_V1_6(n,bird_setp)
n = 30;           % Size of the swarm " no of birds "
bird_setp = 100;   % Maximum number of "birds steps"
dim = 3;          % Dimension of the problem

c2 = 2.05;        % PSO parameter C1 
c1 = 2.05;        % PSO parameter C2 
wo = 0.9;         % pso initial momentum or inertia  
wf = 0.2;         % pso final momentum or inertia  
slope = wf-wo/bird_setp; 
w = wo;
fitness=0*ones(n,bird_setp);


                                       %-----------------------------%
                                       %    initialize the parameter %
                                       %-----------------------------%
% Transfers functions                                   
s = tf('s');
Kp = 10;
Ki = 1;
Kd = 0;
N=100;


% h2d = (ones(10001           ,1))*0.75;
% h1d = (ones(10001           ,1))*0.4219;

h2 = H2';
h2d = H2d;
h1 = H1';
h1d = H1d;

%%



%ALTERAR PARA 3 GANHOS
%   Kp Ki Kd                                        
lb = [50;10;1] ;      % LOWER BOUND - LIMITE INFERIOR
ub = [250;50;10]  ;       % UPPER BOUND - LIMITE SUPERIOR
% ub = [300;100;20]  ;       % UPPER BOUND - LIMITE SUPERIOR

konstante = 200; %constante de ganho da derivada da entrada
% konstante = 150; %constante de ganho da derivada da entrada
% konstante = 175; %constante de ganho da derivada da entrada
R1 = rand(dim, n);
R2 = rand(dim, n);
current_fitness =0*ones(n,1);
                                 %------------------------------------------------%
                                 % Initializing swarm and velocities and position %
                                 %------------------------------------------------%
                                 
current_position =     lb.*ones(dim,n) +(ub - lb).*(rand(dim, n)) ;
velocity = .01*current_position;
local_best_position  = current_position ;


                                 %-------------------------------------------%
                                 %     Evaluate initial population           %           
                                 %-------------------------------------------%
                                 
                                 
                                 
                                 
                                 
                                 
                                 
                                 %-------------------------------------------%
                                 %    Initializing OPSO parameters           %           
                                 %-------------------------------------------%
                                 % OBL parameters
minFit = 0.1;
maxFNC = 100;
f_impr = 1e10;
FNC=0;
opposite=0;
cnt_obl = 0;
Ocontador = 0;

x_max = ub;
x_min = lb;
threshold = 1e-2;
L = x_max+x_min;
k = 1;  % index of iteration
Prev_FITNESS = 0;
f_ind = 1e10*ones(n,1);

Melhor_artigo = cell(30,1);
Melhor_sol = cell(30,1);

    %%
for i = 1:n
    current_fitness(i) = 1E10;%tracklsq(current_position(:,i));    
end


local_best_fitness  = current_fitness ;
[global_best_fitness,g] = min(local_best_fitness) ;
globl_best_position = zeros(dim,n);
for i=1:n
    globl_best_position(:,i) = local_best_position(:,g) ;
end
                                               %-------------------%
                                               %  VELOCITY UPDATE  %
                                               %-------------------%

velocity = w *velocity + c1*(R1.*(local_best_position-current_position)) + c2*(R2.*(globl_best_position-current_position));

                                               %------------------%
                                               %   SWARMUPDATE    %
                                               %------------------%


current_position = current_position + velocity ;

                                               %------------------------%
                                               %  evaluate anew swarm   %
                                               %------------------------%

best_fitness = zeros(bird_setp);
zeros(dim,n);
for i=1:n
    ind=find(current_position(:,i)<lb);
    current_position(ind,i)=lb(ind)+rand*0.3*lb(ind);

    ind=find(current_position(:,i)>ub);
    current_position(ind,i)=ub(ind)-rand*0.3*ub(ind);
end

%% Main Loop
iter = 0 ;        % Iterations’counter
% fprintf('iter=%d, fitness=%3f, Kp=%3f, Kd=%3f, Ki=%3f\n', iter,global_best_fitness,globl_best_position(1,1),globl_best_position(2,1),globl_best_position(3,1)); 
while  ( iter < bird_setp )
    if rem(iter,5) == 0
        fprintf('iter=%d, fitness=%3f, Kp=%3f, Ki=%3f, Kd=%3f \n', iter,global_best_fitness,globl_best_position(1,1),globl_best_position(2,1),globl_best_position(3,1)); 
    end
    iter = iter + 1;
    
    for i = 1:n
        
        ind=find(current_position(:,i)<lb);
        current_position(ind,i)=lb(ind)+rand*0.3*lb(ind);

        ind=find(current_position(:,i)>ub);
        current_position(ind,i)=ub(ind)-rand*0.3*ub(ind);
        
        kP = abs(current_position(1,i));
        kI = abs(current_position(2,i));
        kD = abs(current_position(3,i));
%           p1 = 0.75; P1 = poly(-[p1 p1]); K1 = P1(1,2); K0 = P1(1,3);
        sim('TwoTanks_PID_R2018a');
        variavel_custo_temp=0;

        for kii=1:length(t-1) %SIMULAÇÃO DEVE CONERGIR EM ATÉ 15s
            
            variavel_custo_temp = variavel_custo_temp+ ((3*(abs(h2f(kii)-h2d(kii)))+1.5*(abs(h1f(kii)-h1d(kii)))).^2)*t(kii); % func obj de acordo com artigo que oniram mandou
        end

        for ki = 1:1:(size(ff,1)-1)
            variavel_custo_temp = variavel_custo_temp + konstante*(ff(ki) -  ff(ki+1));  
        end
        current_fitness(i) =  variavel_custo_temp;
    
    end

        



    for i = 1 : n
        if current_fitness(i) < local_best_fitness(i)
            local_best_fitness(i)  = current_fitness(i);  
%             disp("melhor posição: "+ current_position(:,i))
            local_best_position(:,i) = current_position(:,i)   ;
        end   
    end
    



    [current_global_best_fitness,g] = min(local_best_fitness);

    if current_global_best_fitness < global_best_fitness
        global_best_fitness = current_global_best_fitness;
        for i=1:n
            globl_best_position(:,i) = local_best_position(:,g);
        end
    end

    best_fitness(iter)=global_best_fitness;

    if iter >1

        if (best_fitness(iter) == best_fitness(iter-1))
    %         disp(Food_fitness);
            Ocontador = Ocontador +1;
            opposite=0;
                if Ocontador == 5
                    opposite =1;
                    cnt_obl = cnt_obl+1;
                    Ocontador = 0;
                end
        else
            Prev_FITNESS = best_fitness(iter);
            Ocontador = 0;
        end

    end





    R1 = rand(dim, n);
    R2 = rand(dim, n);
    velocity = w *velocity + c1*R1.*(local_best_position-current_position) + c2*R2.*(globl_best_position-current_position);
    current_position = abs(current_position + velocity); 
    for i=1:n
        ind=find(current_position(:,i)<lb);
        current_position(ind,i)=lb(ind)+rand*0.3*lb(ind);

        ind=find(current_position(:,i)>ub);
        current_position(ind,i)=ub(ind)-rand*0.3*ub(ind);
    end


            %implementador do aleatoriazador
    if opposite == 1
        fprintf('iter=%d cnt_obl=%d \n',iter,cnt_obl);
        opposite = 0;
        [alpha, gamma] = sort(local_best_fitness); % usar o gamma
        contador = ceil(size(gamma,1)/3);
        final = gamma(end:-1:(end-contador));

        for i = 1:dim % dimensao
            iterador = 0;
            if rem(randi([1 2]),2)
                for j = 1:1:contador % particula


                    if current_position(i,gamma(j)) >L(i)/2
                        current_position(i,gamma(j)) = abs(L(i) - (x_max(i)-current_position(i,gamma(end+1-j))));
                    elseif ceil(current_position(i,gamma(j))) == ceil(L(i)/2)
                        current_position(i,gamma(j)) = abs(x_max(i) - L(i)*rand());
                    else
                        current_position(i,gamma(j)) = abs(x_max(i) - (current_position(i,gamma(end+1-j)) - x_min(i)));
                    end
                    iterador = iterador +1;

                end
            end
        end


    end            











%     ind=find(current_position<lb);
%     current_position(ind)=lb(ind);
%     ind=find(current_position>ub);
%     current_position(ind)=ub(ind);

    w = w + slope; 
%     fprintf('iter=%d, fitness=%3f, Kpi=%3f, Kii=%3f, Kdi=%3f\n', iter,global_best_fitness,globl_best_position(1,1),globl_best_position(2,1),globl_best_position(3,1)); 
end % end of while loop its mean the end of all step that the birds move it 

    %%
    %             xx=fitness(:,bird_setp);
    %             [Y,I] = min(xx);
    %             current_position(:,I)
    %saida
%

sol=globl_best_position(:,g);


% 

kP = abs(current_position(1,i));
kI = abs(current_position(2,i));
kD = abs(current_position(3,i));
sim('TwoTanks_PID_R2018a'); plotResultados_NL;

variavel_custo_temp=0;

for kii=1:length(t-1) %SIMULAÇÃO DEVE CONERGIR EM ATÉ 15s

    variavel_custo_temp = variavel_custo_temp+ ((3*(abs(h2f(kii)-h2d(kii)))+1.5*(abs(h1f(kii)-h1d(kii)))).^2)*t(kii); % func obj de acordo com artigo que oniram mandou
end

for ki = 1:1:(size(ff,1)-1)
    variavel_custo_temp = variavel_custo_temp + konstante*(ff(ki) -  ff(ki+1));  
end
% 
% plot_simu = sim(simIn);
% tempo_simu = 50;
% plot(plot_simu.yout{2}.Values.Time(1:5000),plot_simu.yout{2}.Values.Data(1:5000))
% % disp("Tempo de subida Malha Interna: "+Dados.RiseTime)
% % disp("Tempo de armotecimento Malha Interna: "+Dados.SettlingTime)
% % disp("Overshoot Interna: "+Dados.Overshoot)
% % %         
% % 
% % %
% % % disp("Tempo de subida Malha Externa: "+Dados.RiseTime)
% % % disp("Tempo de armotecimento Malha Interna: "+Dados.SettlingTime)
% % % disp("Overshoot Externa: "+Dados.Overshoot)
% % % disp("cnt_obl: " + cnt_obl);
% % disp("Terminou OPSO")
% % 
% % Plot_Conv = figure('Name','Curva de convergência');
% % Plot_Conv.Position = [1210 400 600 600];
% % semilogy(best_fitness);
% % tracklsq(globl_best_position(:,g));
% % 
% 
% save('PSO_workspace_PID_v1.0_final.mat')
% save('PSO_workspace_PID')












