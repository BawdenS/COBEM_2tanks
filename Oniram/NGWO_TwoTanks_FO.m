clc; clear all;

%==============================================================================
% INICIALIZA��O DA PLANTA

% Controle por Trajet�ria para o modelo SISO N�o-Linear
% de um Sistema de 2 Tanques Diferencialmente Plano. 

% Nome: STORMS Group
% Data: Maio/2023

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
  time=t;
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

% Par�metros do Controlador VIR�O DO GREY WOLF

% GWO Parameters
Titer = 100;     % NUMERO DE INTERA��ES
Npop = 30;       % TAMANHO DA POPULA��O
var = 2;        % NUMERO DE VARIAVEIS A SEREM OTIMIZADAS

% ESPA�O DE BUSCA
lb = [0.04,0.4] ;      % LOWER BOUND - LIMITE INFERIOR
ub = [0.81,1.8]  ;       % UPPER BOUND - LIMITE SUPERIOR

% PASSOS DE OTIMIZA��O
c_cf = 0;       

% REFER�NCIA
Zd = 0.75 + 0*t; d1Zd = 0*t; d2Zd = 0*t;

H2d = Zd;

% INICIALIZA��O DO ALGORITIMO
for p=1:Npop
    
    for n=1:var
        x(p,n) = lb(n) + rand*(ub(n)-lb(n)); %POSI��O ALEAT�RIA de K0 - xa
    end
    
    % GANHOS A SEREM OTIMIZADOS NO FLATNESS
    K0 = x(p,1);
    K1 = x(p,2);
    
    % SIMULA��O DO MODELO UTILIZANDO OS GANHOS OBTIDOS
    sim('TwoTanks_FO_2018a');
    h2 = h2f';
    h2d = H2d;
    h1  = h1f';
    h1d = H1d;
    d2zd = d2Zd';
    tC = t;
    
    %FUN��O OBJETIVO - ITAE
    obj(p)=0;
    tolerance=0.02;
        idx =  find(abs((h2-0.75))>tolerance, 1, 'last');
    for dim=1:length(time); %SIMULA��O DEVE CONERGIR EM AT� 15s
        % obj(p)=obj(p)+(abs(h2(dim)-h2d(dim)))*time(dim) + time(idx) ; % ITAE+tempo de acomoda��o
        obj(p) = obj(p)+ ((3*(abs(h2(dim)-h2d(dim)))+1.5*(abs(h1(dim)-h1d(dim)))).^2)*time(dim); % func obj de acordo com artigo que oniram mandou
       % obj(p) = obj(p) + ((abs(tC(dim)-d2zd(dim))))*time(dim)
    end
    %obj
end

% ENCONTRAR A MELHOR SOLU��O
[b_sol,b_loc]=min(obj);
b_ITAE=b_sol;
best_k=x(b_loc); % RESULTADOS DOS GANHOS K0 E K1

for tt=1:Titer
    a_gwo=2-((2*tt)/Titer); % F�RMULA DE a
    % BUSCA DO ALPHA, BETA E DELTA:
    %ALPHA
    sol_cf=obj;
    sol_var=x;
    [b_sol1,loc1] = min(sol_cf);
    alph_cf=b_sol1;
    alph_var=sol_var(loc1,:);
    
    %REMOVENDO A MELHOR SOLU��O ENCONTRADA NA POPULA��O PARA BUSCAR O NOVO
    %MINIMO COMO POSI��O DO BETA
    sol_cf(loc1)=[];
    sol_var(loc1,:)=[];
    [b_sol1,loc1]=min(sol_cf);
    beta_cf=b_sol1;
    beta_var=sol_var(loc1,:);
    
    sol_cf(loc1)=[];
    sol_var(loc1,:)=[];
    [b_sol1,loc1]=min(sol_cf);
    delta_cf=b_sol1;
    delta_var=sol_var(loc1,:);
    
    %COEFICIENTES DE BUSCA DO ALGORITIMO
    c_cof=2*rand;
    a_cof=a_gwo*((2*rand)-1);
    
    for p=1:Npop
        
        for n=1:var
            
            d_a=abs(c_cof*alph_var(n)-x(p,n));
            xx1=alph_var(n)-(a_cof*d_a);
            d_b=abs(c_cof*beta_var(n)-x(p,n));
            xx2=beta_var(n)-(a_cof*d_b);
            d_d=abs(c_cof*delta_var(n)-x(p,n));
            xx3=delta_var(n)-(a_cof*d_d);
            
            x(p,n)=(xx1+xx2+xx3)/3;
            
            %CHECK NO ESPA�O DE BUSCA
            if x(p,n)>ub(n); x(p,n)=ub(n);
            elseif x(p,n)<lb(n); x(p,n)=lb(n); end
            
        end
        
        % RODAR A SIMMULA��O NOVAMENTE PARA COMPARAR SE O RESULTADO ATUAL �
        % MELHOR DO QUE O RESULTADO OBTIDO NO PASSO ANTERIOR.
    for n=1:var
        x(p,n) = lb(n) + rand*(ub(n)-lb(n)); %POSI��O ALEAT�RIA de K0 - xa
    end     
    % GANHOS A SEREM OTIMIZADOS NO FLATNESS
    K0 = x(p,1);
    K1 = x(p,2);
    
    % SIMULA��O DO MODELO UTILIZANDO OS GANHOS OBTIDOS
    sim('TwoTanks_FO_2018a'); 
    h2 = h2f';
    h2d = H2d;
    h1 = h1f';
    h1d = H1d;
    tC = tC';
    d2zd = d2Zd';
    %FUN��O OBJETIVO - ITAE
    
  obj(p)=0;
  aux(p)=0;
   tolerance = 0.02;
    idx =  find(abs((h2-0.75))>tolerance, 1, 'last');
    for dim=1:length(time); %SIMULA��O DEVE CONERGIR EM AT� 15s
        %obj(p)=obj(p)+(abs(h2(dim)-h2d(dim)))*time(idx);
       obj(p) = obj(p)+ ((3*(abs(h2(dim)-h2d(dim)))+1.5*(abs(h1(dim)-h1d(dim)))).^2)*time(dim);
      %  obj(p) = obj(p) + ((abs(tC(dim)-d2zd(dim))))*time(dim); 
       % aux(p)=obj(p);
       %obj
    end
    end
  
    %ENCONTRAR A MELHOR SOLU��O
        
    % ENCONTRAR A MELHOR SOLU��O
    [b_sol,b_loc]=min(obj);
    best2_ITAE=b_sol;
    best2_k=x(b_loc); % RESULTADOS DOS GANHOS K0 E K1
    
    if best2_ITAE<b_ITAE
        b_ITAE=best2_ITAE;
        best_k=best2_k;
    end
    
%    PLOT DA FUN��O CUSTO PARA CADA ITERA��O
    c_cf=c_cf+1;
    best_cf_gwo(c_cf)=b_ITAE;
    time1=1:c_cf;
    
    figure(2)
    plot(time1,best_cf_gwo,'k','LineWidth',2)
    xlabel('Iterations')
    ylabel('Average Fitness')
    hold on
    
    tt
end 
result (c_cf) = best_cf_gwo(c_cf) - abs(b_ITAE);
display(['Os m�nimo encontrado da fun��o objetivo foi : ', num2str(b_ITAE)]);
display(['Os melhores ganhos flatness K0 e K1 s�o : ', num2str(K0),' ', num2str(K1)]);
%POR FIM SIMULA-SE A RESPOSTA FINAL DO SISTEMA
sim('TwoTanks_FO_2018a'); 
figure(3)
plotResultados_NL
%hold on
