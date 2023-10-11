% Visualizaçao do Resultados

  % Pré-Comandos
    addpath('./Resultados/Export_fig/');
    close all; LW = 2; flagPrint = 0;

  % Variáveis de Estado (H1 e H2) + Sinal de Control (U)
    a = figure(1);
    subplot(3,1,1); plot(t,H1d,'b','LineWidth', LW); hold on; plot(tf,h1f,'r--','LineWidth', LW);
    % title('Evolução Temporal do Nível H1', 'fontsize', 14);
    xlabel('t [s]', 'fontsize', 18); ylabel('H1 [m]', 'fontsize', 18);
    
    subplot(3,1,2); plot(t,H2d,'b','LineWidth', LW); hold on; plot(tf,h2f,'r--','LineWidth', LW);  
    % title('Evolução Temporal do Nível H2', 'fontsize', 14);
    xlabel('t [s]', 'fontsize', 18); ylabel('H2 [m]', 'fontsize', 18);

    subplot(3,1,3); plot(t,Fd,'b','LineWidth', LW); hold on; plot(tf,ff,'r--','LineWidth', LW);
    % title('Evolução Temporal do Sinal de Controle', 'fontsize', 14);
    xlabel('t [s]', 'fontsize', 18); ylabel('F [Hz]', 'fontsize', 18);
    
    ab = 40;
    a.Position = [150 150 1080 640];
    a.Children(1).XLabel.FontSize =  ab;
    a.Children(1).YLabel.FontSize =  ab;
    a.Children(2).XLabel.FontSize =  ab;
    a.Children(2).YLabel.FontSize =  ab;
    a.Children(3).XLabel.FontSize =  ab;
    a.Children(3).YLabel.FontSize =  ab;
    
    a.Children(1).FontSize = 18;
    a.Children(2).FontSize = 18;
    a.Children(3).FontSize = 18;
    
    if(flagPrint == 1)
%       maximize(gcf); 
      export_fig('./Resultados/2Tanks', '-eps', '-transparent'); 
      eps2pdf('./Resultados/2Tanks.eps', './Resultados/2Tanks_semilogy.pdf', 1);
    end