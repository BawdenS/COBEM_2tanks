clear
clc
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
    P = [0.425,0.75];
    J = jacobian(system,v);
    A=double(subs(J,v,P));
    
    disp(A);

%     B = [0;a1;0;0];
%     C = [0,0,1,0];
    B = [a1;0];
    C = [0,1];
    D = 0;
    
    sysSS = ss(A,B,C,D);
    
    
    
    
    G = C*inv(s*eye(2)-A)*B+D
    [num,den] = ss2tf(A,B,C,D)
    G2 = tf(num,den)