% SBAIN_2TANKS COBEM
% Simulink PSO-PID twotank
% Controle via PID para o modelo SISO Não-Linear
% de um Sistema de 2 Tanques Diferencialmente Plano. 

% Nome: STORMS Group
% Data: Junho/2023

% Pré-Comandos
  close all; 
  clc; 
  clear all;

  
  
  
  
% temp_GWO = load('GWO_runs_refnova.mat'); 
% h1f_gwo = temp_GWO.h1f;
% h2f_gwo = temp_GWO.h2f;
% ff_gwo = temp_GWO.ff;


% temp_PSO = load('PSO_workspace_PID_v1.0_agoravai_refnova.mat'); 
% h1f_pso = temp_PSO.h1f;
% h2f_pso = temp_PSO.h2f;
% ff_pso = temp_PSO.ff;
  


syms h1 h1p h2 h2p F
syms Fnom Q C1 C2

s = tf('s');

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
    num3 = Aout1*sqrt(2*g)*H2max^2; den3 = pi*(R2max^2)*(h2^2); a3 = num3/den3;
    num4 = Aout2*sqrt(2*g)*H2max^2; den4 = pi*R2max^2; a4 = num4/den4;


%     eq0 = h1p;
    eq1 = a1*F - a2*(h1^(1/2));
%     eq2 = h2p;
    eq3 = a3*(h1^(1/2)) - a4*(h2^(-3/2));
    
%     system = [eq0;eq1;eq2;eq3];
%     v =[h1 h1p h2 h2p];
    system = [eq1;eq3];
    v =[h1 h2];
        


%     P = [0.425,0,0.75,0];
%     P = [0.425,0.75];
    P = [0.3,0.5];
    J = jacobian(system,v);
    A=double(subs(J,v,P));
    
%     disp(A);

%     B = [0;a1;0;0];
%     C = [0,0,1,0];
    B = [a1;0];
    C = [0,1];
    D = 0;
    
    sysSS = ss(A,B,C,D);
    
    
    
    
    G = C*inv(s*eye(2)-A)*B+D;
    [num,den] = ss2tf(A,B,C,D);
    G2 = tf(num,den);
[num,den] = tfdata(G2,'v');
  disp(G2)
  
  
  h2p = 0.5;
  
  
  
  
  
  %PID
  K = 0.0533;
  tau1 = 1.844;
  tau2 = 4.662;
% lambda = 0.5;
  lambda = 2;
  KpIMC = (tau1+tau2)/(K*lambda); % 244.1276
  TiIMC = (tau1+tau2);
  TdIMC = (tau1*tau2)/(tau1+tau2);
  KiIMC = KpIMC/TiIMC; % 37.5235
  KdIMC = TdIMC*KpIMC; % 322.5789
  N = 100; % filtro
 
  
  
  controlador = (1/(TiIMC*s) +(1/(100*s+1))*TdIMC*s+1)*KpIMC;
  disp(controlador)
  
  
  
  
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
  h10 = 0.30; h20 = 0.10;
      
% Tempo de Simulação
  t0 = 0; dt = 0.01; tfinal = 300;
  t = t0:dt:tfinal; st = size(t,2);
  
% Planejamento de Trajetória - Saída Plana (Z)

  % Trajetória Constante com 2 valores (STEP)
    fP1 = 0.30; fP2 = 0.60;
    Zd = 0.5 + 0*t; Zd((round(fP1*st)):end) = 0.8; Zd((round(fP2*st)):end) = 0.20;
    d1Zd = 0*t; d2Zd = 0*t;
      
% Planejamento de Trajetória - Variáveis do Sistema
  H2d = Zd;
  
  H1d = (Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2;
     
  Fd = ((2*Zd.^(5/2)*a4 + 2*Zd.^4.*d1Zd)/(a1*a3^2)).*d2Zd + (a4^2*d1Zd + 4*Zd.^3.*d1Zd.^3 + 5*Zd.^(3/2)*a4.*d1Zd.^2 + a2*a3^2*((Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2).^(1/2))/(a1*a3^2);


 


  % Saturação do Atuador
    satF_sup = 90;
    satF_inf = 10;       
        
        
        %% para o PSO
        
KpIMC = 213.7679; % 244.1276
TiIMC = 5;
TdIMC = 0.04;
KiIMC = KpIMC/TiIMC; % 37.5235
KdIMC = TdIMC*KpIMC; % 322.5789
N = 100; % filtro

sim('cobem_IMC');    
        

temp_PSO.h1f= h1_saida.signals(1).values';
temp_PSO.h2f= h2_saida.signals(1).values';
temp_PSO.ff= f_saida.signals(1).values;
                %% para o GWO
        
KpIMC = 22.7335; % 244.1276
TiIMC = 2.8696;
TdIMC = 0.0472;
KiIMC = KpIMC/TiIMC; % 37.5235
KdIMC = TdIMC*KpIMC; % 322.5789
N = 100; % filtro

sim('cobem_IMC');    
        
                
temp_GWO.h1f=h1_saida.signals(1).values';
temp_GWO.h2f=h2_saida.signals(1).values';
temp_GWO.ff=f_saida.signals(1).values;
        
        
        
%%
  
  
  %PID
  K = 0.0533;
  tau1 = 1.844;
  tau2 = 4.662;
% lambda = 0.5;
  lambda = 2;
  KpIMC = (tau1+tau2)/(K*lambda); % 244.1276
  TiIMC = (tau1+tau2);
  TdIMC = (tau1*tau2)/(tau1+tau2);
  KiIMC = KpIMC/TiIMC; % 37.5235
  KdIMC = TdIMC*KpIMC; % 322.5789
  N = 100; % filtro


  
% % Executando a simulação e apresentando os resultados
%   sim('cobem_IMC'); 
%   plotResultados_NL;


  sim('cobem_IMC'); 
%   plotResultados_NL;


obj = 0;
h1 = h1_saida.signals(1).values';
h2 = h2_saida.signals(1).values';
ff = f_saida.signals(1).values;
        for dim=1:length(t) %SIMULAÇÃO DEVE CONERGIR EM ATÉ 15s
            %obj(p)=obj(p)+(abs(h2(dim)-h2d(dim)))*time(idx);
           obj = obj+ ((3*(abs(h2(dim)-H2d(dim)))+1.5*(abs(h1(dim)-H1d(dim)))).^2)*t(dim);
          %  obj(p) = obj(p) + ((abs(tC(dim)-d2zd(dim))))*time(dim); 
           % aux(p)=obj(p);
           %obj
        end
        k = 200;
        for i = 1:1:(size(ff,1)-1)
            obj = obj + k*abs(ff(i) -  ff(i+1));  
        end       
        
        
        
% fprintf("Valor da função custo:"+ string(obj))
