% SBAIN_2TANKS
% Simulink PSO-PID twotank
clear
close all 
clc
% Par�metros do Sistema N�o-Linear   
 
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
    
    
    
% Condi��es Iniciais do Sistema
  h10 = 0.25; h20 = 0.25;

% Tempo de Simula��o
  t0 = 0; dt = 0.01; tfinal = 25;
  t = t0:dt:tfinal; st = size(t,2);
  
% Planejamento de Trajet�ria - Sa�da Plana (Z)

  % Trajet�ria Constante 
    Zd = 0.75 + 0*t; d1Zd = 0*t; d2Zd = 0*t;
  
  % Trajet�ria Cosseinodal
    % ro = 0.25; w = 0.250;    
    % Zd = 0.50 + ro*cos(w*t); d1Zd = -ro*w*sin(w*t); d2Zd = -ro*(w^2)*cos(w*t);
  
% Planejamento de Trajet�ria - Vari�veis do Sistema
  H2d = Zd;
  
  H1d = (Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2;
     
  Fd = ((2*Zd.^(5/2)*a4 + 2*Zd.^4.*d1Zd)/(a1*a3^2)).*d2Zd + (a4^2*d1Zd + 4*Zd.^3.*d1Zd.^3 + 5*Zd.^(3/2)*a4.*d1Zd.^2 + a2*a3^2*((Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2).^(1/2))/(a1*a3^2);

  tempo_simu = 25;





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
dim = 2;          % Dimension of the problem

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
Kp1 = 1;
Ki1 = 1;


% IDEAL polos padr�o Bessel
% caso Gabriel
h2_G = (ones(2501           ,1))*0.75;
h1_G = (ones(2501           ,1))*0.4219;

% Caso altura H2
V_ideal = 0.75;
G_ideal  = 21.9/((s+4.0530 - 2.34i)*(s+4.0530 + 2.34i));
[GInum,GIden] = tfdata(G_ideal,'v');
Saida_ideal = step(G_ideal);
dado_saida = stepinfo(G_ideal);
tempo_ideal = tempo_simu; % Passo simu 1e-3 e 
nome_g_ideal = 'cobem_g_ideal'; % nome arquivo simulink para comparar
sim(nome_g_ideal); % simula��o do sistema com sa�da ideal
Vetor_obj = saida_ideal{1}.Values.Data; % sera usado na func custo
[A,B] = sort(Vetor_obj);
Overshoot = abs(Vetor_obj(B(end))-Vetor_obj(end))/Vetor_obj(B(end))*100; % calculo do overshoot em %
Temporario = Vetor_obj(end-1000:end)-Vetor_obj(end-1001:end-1);
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


% Caso altura H1
V_ideal = 0.4219;
saida_ideal = step(G_ideal);
dado_saida = stepinfo(G_ideal);
tempo_ideal = tempo_simu; % Passo simu 1e-2 e 
nome_g_ideal = 'cobem_g_ideal'; % nome arquivo simulink para comparar
sim(nome_g_ideal); % simula��o do sistema com sa�da ideal
Vetor_obj_H1 = saida_ideal{1}.Values.Data; % sera usado na func custo
[A,B] = sort(Vetor_obj_H1);
Overshoot_h1 = abs(Vetor_obj_H1(B(end))-Vetor_obj_H1(end))/Vetor_obj_H1(B(end))*100; % calculo do overshoot em %
Temporario = Vetor_obj_H1(end-1000:end)-Vetor_obj_H1(end-1001:end-1);
Cond_acom_h1 = mean(Temporario); % confere se acomodou

%medicao tempo de acm
passo_tempo = 50;
for i = passo_tempo+1:passo_tempo:size(Vetor_obj_H1,1)
   if abs(mean(Vetor_obj_H1(i-passo_tempo:i) - Vetor_obj_H1(i-passo_tempo+1:i+1))) < 0.00005
       tempo_acm_h1 = saida_ideal{1}.Values.Time(i); % pega tempo de acomoda��o
   break;
   end
end

% medi�ao tempo de subida
Temporario = 1; %flag
for i = 1:1:size(Vetor_obj_H1,1)
   if ((Vetor_obj_H1(i) > Vetor_obj_H1(end)*0.1) && (Temporario == 1))
       Temporario=0;
       tempo_rising_0 = saida_ideal{1}.Values.Time(i); % pega tempo de subida min
   end
   if Vetor_obj_H1(i) > Vetor_obj_H1(end)*0.9
       tempo_rising_1 = saida_ideal{1}.Values.Time(i); % pega tempo de subida max
   break;
   end
end
tempo_rising_h1 = tempo_rising_1 - tempo_rising_0; % tempo de subida


%%




%   K0 K1                                       
ub=[0.81 1.8]'; %/*upper bounds of the parameters. */
lb=[0.04 0.4]';%/*lower bound of the parameters.*/  

% %   K0 K1                                       
% ub=[10 10]'; %/*upper bounds of the parameters. */
% lb=[0.01 0.01]';%/*lower bound of the parameters.*/

R1 = rand(dim, n);
R2 = rand(dim, n);
current_fitness =0*ones(n,1);
                                 %------------------------------------------------%
                                 % Initializing swarm and velocities and position %
                                 %------------------------------------------------%
                                 
current_position = ub.*(rand(dim, n));
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
iter = 0 ;        % Iterations�counter
% fprintf('iter=%d, fitness=%3f, Kp=%3f, Kd=%3f, Ki=%3f\n', iter,global_best_fitness,globl_best_position(1,1),globl_best_position(2,1),globl_best_position(3,1)); 
while  ( iter < bird_setp )
    if rem(iter,5) == 0
        fprintf('iter=%d, fitness=%3f, Kp=%3f, Ki=%3f \n', iter,global_best_fitness,globl_best_position(1,1),globl_best_position(2,1)); 
    end
    iter = iter + 1;
    
    for i = 1:n
        
        ind=find(current_position(:,i)<lb);
        current_position(ind,i)=lb(ind)+rand*0.3*lb(ind);

        ind=find(current_position(:,i)>ub);
        current_position(ind,i)=ub(ind)-rand*0.3*ub(ind);
        
        K0 = abs(current_position(1,i));
        K1 = abs(current_position(2,i));
%           p1 = 0.75; P1 = poly(-[p1 p1]); K1 = P1(1,2); K0 = P1(1,3);
  
      % Executando a simula��o e apresentando os resultados
        sim('TwoTanks_FO_2018a'); 
          %plotResultados_NL;
        
        vetor_H1 = h1f; %simulink -> workspace 
        custo_area_erro_H1= abs(h1_G-vetor_H1);
%         [A,B] = sort(vetor_H1);
%         Overshoot_atual = abs(vetor_H1(B(end))-vetor_H1(end))/vetor_H1(B(end)) *100; % calculo do overshoot em %
%         custo_pico_H1= 10* Overshoot_atual;  
        
        
        
        
        
        
        
        
        
        
        
          
        vetor_H2 = h2f;        %simulink -> workspace 
        custo_area_erro= abs(h2_G-vetor_H2);
        
        
        
%         [A,B] = sort(vetor_H2);
%         Overshoot_atual = abs(vetor_H2(B(end))-vetor_H2(end))/vetor_H2(B(end)) *100; % calculo do overshoot em %
%         custo_pico= 10* Overshoot_atual;
        
%     temp= 0;
%     for iter_temp = 1:length(t)
%     temp = temp+ ((3*(abs(h2_G(iter_temp)-vetor_H2(iter_temp)))+1.5*(abs(h1_G(iter_temp)-vetor_H1(iter_temp))))^2)*t(iter_temp);
% 
%     end


%         current_fitness(i) = 2*custo_area_erro + 2*custo_pico + 1*custo_area_erro_H1 + 1*custo_pico_H1;
        current_fitness(i) =  (t*((3*custo_area_erro + 1.5*custo_area_erro_H1).^2) )  ;
%         current_fitness(i) = temp;
%         temp_string = ['Valor do erro do Calculo: ',num2str(temp-current_fitness(i))];
%         disp(temp_string);
    end

    for i = 1 : n
        if current_fitness(i) < local_best_fitness(i)
            local_best_fitness(i)  = current_fitness(i);  
%             disp("melhor posi��o: "+ current_position(:,i))
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
K0 = sol(1);
K1 = sol(2);
sim('TwoTanks_FO_2018a'); plotResultados_NL;
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
% % Plot_Conv = figure('Name','Curva de converg�ncia');
% % Plot_Conv.Position = [1210 400 600 600];
% % semilogy(best_fitness);
% % tracklsq(globl_best_position(:,g));
% % 
% 














