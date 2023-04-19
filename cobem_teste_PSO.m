% Simulink PSO-PID twotank
clear
close all 
clc
% Parametros
h10 = 0.5;
h20 = 0.5;
pho = 998.23; % kg/(m^3)
Qe = 31.416*pho; % (m^3)/s *kg/(m^3) = kg/s

Fnom = 60; % Hz
g = 9.81; % m/(s^3)
Area1 = 3.1416; % m^2
Area2 = 0.0113; % m^2
Area4 = 0.0314; % m^2
Rmax = 1; % m
Hmax = 3; % m
%

modelname = 'cobem_simulink_alpha1_0';
simIn = Simulink.SimulationInput(modelname);

% simIn = setVariable(simIn,'a',10);

tempo_simu = 20;








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
n = 50;           % Size of the swarm " no of birds "
bird_setp = 10;   % Maximum number of "birds steps"
dim = 3;          % Dimension of the problem

c2 = 2.05;        % PSO parameter C1 
c1 = 2.05;        % PSO parameter C2 
wo = 0.5;         % pso initial momentum or inertia  
wf = 0.1;         % pso final momentum or inertia  
slope = wf-wo/bird_setp; 
w = wo;
fitness=0*ones(n,bird_setp);


                                       %-----------------------------%
                                       %    initialize the parameter %
                                       %-----------------------------%
% Transfers functions                                   
s = tf('s');
Kp1 = 1;
Ki1 = 1;
Kd1 = 0;
t_filtro = 0.0426;


% IDEAL polos padr�o Bessel
G_ideal  = 21.9/((s+4.0530 - 2.34i)*(s+4.0530 + 2.34i));
[GInum,GIden] = tfdata(G_ideal,'v');
Saida_ideal = step(G_ideal);
dado_saida = stepinfo(G_ideal);
tempo_ideal = 21; % Passo simu 1e-3 e 
nome_g_ideal = 'cobem_g_ideal'; % nome arquivo simulink para comparar
sim(nome_g_ideal); % simula��o do sistema com sa�da ideal
Vetor_obj = saida_ideal{1}.Values.Data; % sera usado na func custo
[A,B] = sort(Vetor_obj);
Overshoot = abs(Vetor_obj(B(end))-Vetor_obj(end))/Vetor_obj(B(end))*100; % calculo do overshoot em %
Temporario = ((Vetor_obj((1000:end)))-(Vetor_obj((1000-1:end-1))));
Cond_acom = mean(Temporario); % confere se acomodou

%medicao tempo de acm
passo_tempo = 50;
for i = passo_tempo+1:passo_tempo:size(Vetor_obj,1)
   if abs(mean(Vetor_obj(i-passo_tempo:i) - Vetor_obj(i-passo_tempo+1:i+1))) < 0.00005
       tempo_acm = saida_ideal{1}.Values.Time(i); % pega tempo de acomoda��o
   break;
   end
end

% medi�ao tempo de subida
Temporario = 1; %flag
for i = 1:1:size(Vetor_obj,1)
   if ((Vetor_obj(i) > Vetor_obj(end)*0.1) && (Temporario == 1))
       Temporario=0;
       tempo_rising_0 = saida_ideal{1}.Values.Time(i); % pega tempo de subida min
   end
   if Vetor_obj(i) > Vetor_obj(end)*0.9
       tempo_rising_1 = saida_ideal{1}.Values.Time(i); % pega tempo de subida max
   break;
   end
end
tempo_rising = tempo_rising_1 - tempo_rising_0; % tempo de subida
%





%   P  I  D                                       
ub=[20 10 3]'; %/*upper bounds of the parameters. */
lb=[0.1 0.0 0]';%/*lower bound of the parameters.*/ 

R1 = rand(dim, n);
R2 = rand(dim, n);
current_fitness =0*ones(n,1);
                                 %------------------------------------------------%
                                 % Initializing swarm and velocities and position %
                                 %------------------------------------------------%
                                 
current_position = ub.*(rand(dim, n));
velocity = .3*current_position;
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
%% Main Loop
iter = 0 ;        % Iterations�counter
% fprintf('iter=%d, fitness=%3f, Kp=%3f, Kd=%3f, Ki=%3f\n', iter,global_best_fitness,globl_best_position(1,1),globl_best_position(2,1),globl_best_position(3,1)); 
while  ( iter < bird_setp )
    if rem(iter,10) == 0
        fprintf('iter=%d, fitness=%3f, Kp=%3f, Ki=%3f, Kd=%3f\n', iter,global_best_fitness,globl_best_position(1,1),globl_best_position(2,1),globl_best_position(3,1)); 
    end
    iter = iter + 1;
    disp(current_position) % retirar
    for i = 1:n
        Kp = abs(current_position(1,i));
        Ki = abs(current_position(2,i));
        Kd = abs(current_position(3,i));
        simu_atual = sim(simIn);
        vetor_comparar = simu_atual.yout{2}.Values.Data(1:20000);
        
        
        custo_area_erro= sum(abs(Vetor_obj-vetor_comparar));
        
        
        
        [A,B] = sort(vetor_comparar);
        Overshoot_atual = abs(vetor_comparar(B(end))-vetor_comparar(end))/vetor_comparar(B(end)) *100; % calculo do overshoot em %
        custo_pico= 10* Overshoot_atual;
        
        
        Temporario = ((Vetor_obj((1000:end)))-(Vetor_obj((1000-1:end-1))));
        custo_acom = 0;
        Cond_acom_atual = mean(Temporario); % confere se acomodou
        if Cond_acom_atual > 0.05
            custo_acom = 1000;
        end

        current_fitness(i) = custo_area_erro + custo_pico + custo_acom;
    end

    for i = 1 : n
        if current_fitness(i) < local_best_fitness(i)
            local_best_fitness(i)  = current_fitness(i);  
            disp("melhor posi��o: ", current_position(:,i))
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
        current_position(ind,i)=lb(ind)+rand*ub(ind);
    end
    for i=1:n
        ind=find(current_position(:,i)>ub);
        current_position(ind,i)=ub(ind)-rand*ub(ind);
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



kpplot = sol(1);
kiplot = sol(2);
kdplot = sol(3);

plot_simu = sim(simIn);
tempo_simu = 50;
plot(plot_simu.yout{2}.Values.Time(1:5000),plot_simu.yout{2}.Values.Data(1:5000))
% disp("Tempo de subida Malha Interna: "+Dados.RiseTime)
% disp("Tempo de armotecimento Malha Interna: "+Dados.SettlingTime)
% disp("Overshoot Interna: "+Dados.Overshoot)
% %         
% 
% %
% % disp("Tempo de subida Malha Externa: "+Dados.RiseTime)
% % disp("Tempo de armotecimento Malha Interna: "+Dados.SettlingTime)
% % disp("Overshoot Externa: "+Dados.Overshoot)
% % disp("cnt_obl: " + cnt_obl);
% disp("Terminou OPSO")
% 
% Plot_Conv = figure('Name','Curva de converg�ncia');
% Plot_Conv.Position = [1210 400 600 600];
% semilogy(best_fitness);
% tracklsq(globl_best_position(:,g));
% 















