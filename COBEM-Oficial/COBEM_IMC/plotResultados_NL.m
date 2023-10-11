% Visualizaçao do Resultados

  % Pré-Comandos
    close all; LW = 0;
% 
%   % Variáveis de Estado (H1 e H2) + Sinal de Control (U)
%     figure(1);
%     subplot(3,1,1); plot(t,H1d,'b','LineWidth', LW); hold on; plot(tf,h1f,'r--','LineWidth', LW);
%     % title('Evolução Temporal do Nível H1', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('H1 [m]', 'fontsize', 14);
%     
%     subplot(3,1,2); plot(t,H2d,'b','LineWidth', LW); hold on; plot(tf,h2f,'r--','LineWidth', LW);  
%     % title('Evolução Temporal do Nível H2', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('H2 [m]', 'fontsize', 14);
% 
%     subplot(3,1,3); plot(t,Fd,'b','LineWidth', LW); hold on; plot(tf,ff,'r--','LineWidth', LW);
%     % title('Evolução Temporal do Sinal de Controle', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('F [Hz]', 'fontsize', 14);
% Visualizaçao do Resultados

%   % Pré-Comandos
%     close all; LW = 2;
% % 
% %   % Variáveis de Estado (H1 e H2) + Sinal de Control (U)
% %     figure(1);
% %     subplot(3,1,1); plot(t,H1d,'b','LineWidth', LW); hold on; plot(tf,h1f,'r--','LineWidth', LW);
% %     % title('Evolução Temporal do Nível H1', 'fontsize', 14);
% %     xlabel('t [s]', 'fontsize', 18); ylabel('H1 [m]', 'fontsize', 14);
% %     
% %     subplot(3,1,2); plot(t,H2d,'b','LineWidth', LW); hold on; plot(tf,h2f,'r--','LineWidth', LW);  
% %     % title('Evolução Temporal do Nível H2', 'fontsize', 14);
% %     xlabel('t [s]', 'fontsize', 18); ylabel('H2 [m]', 'fontsize', 14);
% % 
% %     subplot(3,1,3); plot(t,Fd,'b','LineWidth', LW); hold on; plot(tf,ff,'r--','LineWidth', LW);
% %     % title('Evolução Temporal do Sinal de Controle', 'fontsize', 14);
% %     xlabel('t [s]', 'fontsize', 18); ylabel('F [Hz]', 'fontsize', 14);
% % Visualiza?ao do Resultados
% 
%   % Pr?-Comandos
%     addpath('./Resultado/Export_fig');
%     close all; LW = 2; flagPrint = 0;
% 
%   % Vari?veis de Estado (H1 e H2) + Sinal de Control (U)
%     a = figure(1);
%     subplot(3,1,1); plot(t,H1d,'b','LineWidth', LW); hold on; plot(tf,h1f,'r--','LineWidth', LW);
%     % title('Evolu??o Temporal do N?vel H1', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('H1 [m]', 'fontsize', 18);
%     
%     subplot(3,1,2); plot(t,H2d,'b','LineWidth', LW); hold on; plot(tf,h2f,'r--','LineWidth', LW);  
%     % title('Evolu??o Temporal do N?vel H2', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('H2 [m]', 'fontsize', 18);
% 
%     subplot(3,1,3); plot(t,Fd,'b','LineWidth', LW); hold on; plot(tf,ff,'r--','LineWidth', LW);
%     % title('Evolu??o Temporal do Sinal de Controle', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('F [Hz]', 'fontsize', 18);
%     
%     ab = 40;
%     a.Position = [150 150 1080 640];
%     a.Children(1).XLabel.FontSize =  ab;
%     a.Children(1).YLabel.FontSize =  ab;
%     a.Children(2).YTick = [0 0.2 0.5 0.8 1]
%     a.Children(2).XLabel.FontSize =  ab;
%     a.Children(2).YLabel.FontSize =  ab;
%     a.Children(3).YTick = [0 0.2 0.5 0.8 1]
%     a.Children(3).XLabel.FontSize =  ab;
%     a.Children(3).YLabel.FontSize =  ab;
%     
%     a.Children(1).FontSize = 14;
%     a.Children(2).FontSize = 14;
%     a.Children(3).FontSize = 14;
%     
%     if(flagPrint == 1)
% %       maximize(gcf);
%       export_fig('./Resultados/2Tanks', '-eps', '-transparent'); 
%       eps2pdf('./Resultados/2Tanks.eps', './Resultados/2Tanks_PSO.pdf', 1);
%     end










% Visualizaçao do Resultados

  % Pré-Comandos
%     close all; 
    LW = 2;
% 
%   % Variáveis de Estado (H1 e H2) + Sinal de Control (U)
%     figure(1);
%     subplot(3,1,1); plot(t,H1d,'b','LineWidth', LW); hold on; plot(tf,h1f,'r--','LineWidth', LW);
%     % title('Evolução Temporal do Nível H1', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('H1 [m]', 'fontsize', 14);
%     
%     subplot(3,1,2); plot(t,H2d,'b','LineWidth', LW); hold on; plot(tf,h2f,'r--','LineWidth', LW);  
%     % title('Evolução Temporal do Nível H2', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('H2 [m]', 'fontsize', 14);
% 
%     subplot(3,1,3); plot(t,Fd,'b','LineWidth', LW); hold on; plot(tf,ff,'r--','LineWidth', LW);
%     % title('Evolução Temporal do Sinal de Controle', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('F [Hz]', 'fontsize', 14);
% Visualiza?ao do Resultados

  % Pr?-Comandos
    addpath('./Resultados/Export_fig');
    close all; LW = 2; flagPrint = 1;

    str = '#F58C06';
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
  % Vari?veis de Estado (H1 e H2) + Sinal de Control (U)
    a = figure(1);
    subplot(3,1,1); 
    plot(t,H1d,'b','LineWidth', LW); 
    hold on; 
    plot(tf,h1f,'r--','LineWidth', LW);
    hold on; 
    plot(tf,temp_PSO.h1f,'Color',color,'LineStyle','--','LineWidth', LW);
    hold on; 
    plot(tf,temp_GWO.h1f,'black--','LineWidth', LW);
    legend({'Reference','Signal IMC','Signal PSO','Signal GWO'},'fontsize', 12 )


    ylabel('H1 [m]', 'fontsize', 18);
    ylim([0 1])
    yticks([0 0.2 0.5 0.8 1])
    
    
    
    subplot(3,1,2); plot(t,H2d,'b','LineWidth', LW); 
    hold on; 
    plot(tf,h2f,'r--','LineWidth', LW);  
    hold on; 
    plot(tf,temp_PSO.h2f,'Color',color,'LineStyle','--','LineWidth', LW);
    hold on; 
    plot(tf,temp_GWO.h2f,'black--','LineWidth', LW);
    legend({'Reference','Signal IMC','Signal PSO','Signal GWO'},'fontsize', 12 )
 
    ylabel('H2 [m]', 'fontsize', 18);
    ylim([0 1])
    yticks([0 0.2 0.5 0.8 1])
    subplot(3,1,3); plot(t,Fd,'b','LineWidth', LW); 
    hold on; 
    plot(tf,ff,'r--','LineWidth', LW);
        hold on; 
    plot(tf,temp_PSO.ff,'Color',color,'LineStyle','--','LineWidth', LW);
    hold on; 
    plot(tf,temp_GWO.ff,'black--','LineWidth', LW);
    legend({'Reference','Signal IMC','Signal PSO','Signal GWO'},'fontsize', 12 )
    xlabel('t [s]', 'fontsize', 18); ylabel('F [Hz]', 'fontsize', 18);
    
    ab = 40;
    a.Position = [150 150 1080 640];
%     a.Children(1).XLabel.FontSize =  ab;
%     a.Children(1).YLabel.FontSize =  ab;
%     a.Children(1).YLim = [0 100]
    
%     a.Children(2).YTick = [0 0.2 0.5 0.8 1]
%     a.Children(2).XLabel.FontSize =  ab;
%     a.Children(2).YLabel.FontSize =  ab;
%     a.Children(2).YLim = [0 1]
%     a.Children(3).YTick = [0 0.2 0.5 0.8 1]
%     a.Children(3).XLabel.FontSize =  ab;
%     a.Children(3).YLabel.FontSize =  ab;
%     a.Children(3).YLim = [0 1]
    
%     a.Children(1).FontSize = 12;
%     a.Children(2).FontSize = 12;
%     a.Children(3).FontSize = 12;
%     
    if(flagPrint == 1)
%       maximize(gcf);
      export_fig('./Resultados/2Tanks', '-eps', '-transparent'); 
      eps2pdf('./Resultados/2Tanks.eps', './Resultados/2Tanks_FINAL.pdf', 1);
    end



























