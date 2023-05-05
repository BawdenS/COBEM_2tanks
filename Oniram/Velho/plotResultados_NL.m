% Visualiza�ao do Resultados

  % Pr�-Comandos
    close all; LW = 2;

  % Vari�veis de Estado (H1 e H2) + Sinal de Control (U)
    figure(1);
    subplot(3,1,1); plot(t,H1d,'b','LineWidth', LW); hold on; plot(tf,h1f,'r--','LineWidth', LW);
    % title('Evolu��o Temporal do N�vel H1', 'fontsize', 14);
    xlabel('t [s]', 'fontsize', 18); ylabel('H1 [m]', 'fontsize', 14);
    
    subplot(3,1,2); plot(t,H2d,'b','LineWidth', LW); hold on; plot(tf,h2f,'r--','LineWidth', LW);  
    % title('Evolu��o Temporal do N�vel H2', 'fontsize', 14);
    xlabel('t [s]', 'fontsize', 18); ylabel('H2 [m]', 'fontsize', 14);

    subplot(3,1,3); plot(t,Fd,'b','LineWidth', LW); hold on; plot(tf,uf,'r--','LineWidth', LW);
    % title('Evolu��o Temporal do Sinal de Controle', 'fontsize', 14);
    xlabel('t [s]', 'fontsize', 18); ylabel('U [?]', 'fontsize', 14);