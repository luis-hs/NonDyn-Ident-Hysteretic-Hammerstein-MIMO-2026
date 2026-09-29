% Versão usada no exemplo 1 do artigo ressubmetido 2026. Esse usa o modelo
% GPI em cascata com EE, apesar disso o processo é GRDPI + EE.
%% Ensaio quase estatico que gera dados para aplicar no MPSO

clc
clear
close all

rng(1) 
nu = 2;              % numero de entradas
FR = 24;             % Fator de SNR para o ruído
R = 10;              % Qtd execuções nuvem
N_pesos = 3;        % Numero de disparos
Nv = N_pesos+5;      % Numero total de parametros
VarSize = [2 Nv];    % Size of Decision Variables Matrix
VarMin = 0;          % 0 Lower Bound of Variables
VarMax = 3;          % 1 Upper Bound of Variables
% PSO Parameters

alpha = 0.7;         % Mutation probability
Kmax = 150;          % 1000 Maximum Number of Iterations
P = 30;              % Population Size (Swarm Size)
b = ceil(0.05*P);
Ts = 0.001;

%% identificar dinâmica

num_bloco_hankel = 2;
Amp = 2.3;
f1td = 10;
f2td = 30;
f3td = 50;

f1td2 = 20;
f2td2 = 40;
f3td2 = 60;
Nd = ceil(40/min(f1td,f2td)/Ts);
t = 0:Ts:Nd*Ts - Ts;
n_transiente = ceil(1/min(f1td,f2td)/Ts);

utd = [Amp*sin(2*pi*f1td.*t)+Amp*sin(2*pi*f2td.*t)+Amp*sin(2*pi*f3td.*t);
       Amp*sin(2*pi*f1td2.*t)+Amp*sin(2*pi*f2td2.*t)+Amp*sin(2*pi*f3td2.*t)];
% utd(:,1:500) = [5*sin(2*pi*f1td.*t(1:500));5*sin(2*pi*f1td2.*t(1:500))];
ytd = Processo(length(utd),utd,Ts,2);
vtd_o = Processo_v(length(utd),utd,Ts,2);

%% Validação da dinâmica

f1d = 10;
N_val = ceil(3/f1d/Ts);

Amp = 5;
t = 0:Ts:N_val*Ts - Ts;

uvd = [Amp*sin(2*pi*f1d.*t + pi); Amp*sin(2*pi*f1d.*t + pi)];
yvd = Processo(N_val,uvd,Ts,2);
vvd = Processo_v(length(uvd),uvd,Ts,2);
% figure
% plot(uvd',vvd')
% figure
% plot(yvd')

n_transientv = 1/f1d/Ts;


%% Identificação da histerese
L = 100;
Q = 20;

ut = [];
Amp = 5;
amp = -Amp;
for k = 1:Q
    ut = [ut (Amp-2*abs(amp))*ones(2,L)];
    amp = amp + Amp/(Q/2);
end

ut = [ut(1,:) ut(1,:) zeros(1,2*size(ut(1,:),2));
      zeros(1,2*size(ut(2,:),2)) ut(2,:) ut(2,:)];
yt = Processo(length(ut),ut,Ts,1);

uta = ut(:,L:L:end);
% yta = yt(:,L:L:end);

n_transienth = 25;
uta1 = uta(1,1:2*Q);
uta2 = uta(2,length(uta)/2+1:length(uta)/2+2*Q);
u = [uta1; uta2];


t = 0:Ts:length(u)*Ts - Ts;
Ne = length(t);

%% Validação regime permanente degrau

% Construir sinal: senoidal por 2 ciclos + degrau retido
f_test = f1d;
N_sin  = ceil(1.6/f_test/Ts);
N_hold = ceil(1/f_test/Ts);    % tempo de retenção (1 ciclo)
u_hold_val = 2.67913;          % valor retido (escolher dentro da faixa operacional)

t_test = 0:Ts:(N_sin + N_hold)*Ts - Ts;
u_test_1 = [Amp*sin(2*pi*f_test*(0:Ts:N_sin*Ts-Ts) + pi), ...
             u_hold_val*ones(1, N_hold)];
u_test = [u_test_1; u_test_1];

% Processo real
y_proc_test = Processo(length(u_test), u_test,Ts, 2);

% 1000 modelos - acumular saídas
y_hold_all = zeros(R, 2, length(u_test));  % saídas durante retenção


%% PSO Parameters
phi1 = 2.05;
phi2 = 2.05;
phi = phi1+phi2;
chi = 2/abs(2-phi-sqrt(phi^2-4*phi));
c1 = phi1;       % Personal Learning Coefficient
c2 = phi2;       % Global Learning Coefficient
% Velocity Limits
VelMax = 0.05*VarMax;
%% Inicializa os parâmetros da nuvem
y1s = -inf(1,length(yvd));
y2s = y1s;
y1i = -y1s;
y2i = -y1s;
srmse = 0;
autovalor = zeros(2,R);
save_rmse = zeros(2,R);
erro = zeros(R,2,101);
save_Theta = zeros(VarSize(1),VarSize(2),R);
save_A = zeros(2,2,R);
save_B = zeros(2,2,R);
save_C = zeros(2,2,R);
save_D = zeros(2,2,R);

st = 0;
for realiza = 1:R
    tic
    yo = yt + 0*max(yt,[],2)/(1.1*FR).*randn(size(yt));
    l = L;
    FD = 1;
    for k = 1:length(uta)
        yo(:,k) = sum(yo(:,l-FD:l),2)/FD;
        l = l + L;
    end
    yta = yo(:,1:length(uta));
    yta1 = yta(1,1:2*Q);
    yta2 = yta(2,length(uta)/2+1:length(uta)/2+2*Q);
    y = [yta1; yta2];
    %% Initialization
    % Parametros de iteracao
    Theta = zeros([VarSize P]);
    V = zeros([VarSize P]);
    F = inf(VarSize(1),P);
    % Parametros globais
    Thetag = zeros(VarSize);
    Vg = zeros(VarSize);
    Fg = inf(VarSize(1),1);
    % parametros pessoais
    Thetap = zeros([VarSize P]);
    Vp = zeros([VarSize P]);
    Fp = inf(VarSize(1),P);
    
    for i = 1:P
        Theta(:,:,i) = unifrnd(0,VarMax,VarSize);
        Theta(:,N_pesos+1:end,i) = min(Theta(:,N_pesos+1:end,i),1);
        V(:,:,i) = unifrnd(-VelMax,VelMax,VarSize);
        
        % Evaluation
        F(:,i) = Custo(Ne,N_pesos,Theta(:,:,i),u,y,n_transienth,Ts);
        % Update Personal Best
        Thetap(:,:,i) = Theta(:,:,i);
        Fp(:,i) = F(:,i);
        
        % Update Global Best
        for j = 1:size(u,1)
            if Fp(j,i) < Fg(j)
                Thetag(j,:) = Thetap(j,:,i);
                Fg(j) = Fp(j,i);
            end
        end
    end
    
    
    %% PSO Main Loop de acordo com o Yang2013
    mark = 0;
    
    Mp = zeros([VarSize P]);
    
%     figure
%     hold on
    for k = 1:Kmax
        mark = mark + 1;
        Mg = zeros(VarSize);
        for l = 1:VarSize(1) % Para cada uma das entradas
            % Evaluate weighted personal and global best position
            [ord_cust, ord_ind] = sort(Fp(l,:));
            for q = 1:b
                mu = (1/ord_cust(q))/sum(1./ord_cust(1:b));
                Mg(l,:) = Mg(l,:) + mu*Thetap(l,:,ord_ind(q));
            end
            for i = 1:P
                if i ~= 1
                    Mp(l,:,ord_ind(i)) = (Fp(l,ord_ind(i))*Thetap(l,:,ord_ind(i-1))...
                        + Fp(l,ord_ind(i-1))*Thetap(l,:,ord_ind(i)))/ ...
                        (Fp(l,ord_ind(i))+Fp(l,ord_ind(i-1)));
                else
                    Mp(l,:,ord_ind(i)) = Thetap(l,:,ord_ind(i));
                end
            end
        end
        
        for i = 1:P
            % Update Velocity
            V(:,:,i) = chi*(V(:,:,i) ...
                + c1*rand(VarSize(1),1).*(Mp(:,:,i) - Theta(:,:,i)) ...
                + c2*rand(VarSize(1),1).*(Mg - Theta(:,:,i)));
            
            % Apply Velocity Limits
            V(:,:,i) = max(V(:,:,i),-VelMax);
            V(:,:,i) = min(V(:,:,i),VelMax);
            
            % Update Position
            Theta(:,:,i) = Theta(:,:,i) + V(:,:,i);
            
            for j = 1:VarSize(2)
                for l = 1:VarSize(1)
                    if (Theta(l,j,i) < VarMin || Theta(l,j,i) > VarMax)
                        % Velocity Mirror Effect if var limits are exceed
                        V(l,j,i) = -V(l,j,i);
                        % Apply Position Limits
                        Theta(l,j,i) = rand*VarMin+rand*VarMax;
                    end
                end
            end
            Theta(:,:,i) = max(Theta(:,:,i),VarMin);
            Theta(:,:,i) = min(Theta(:,:,i),VarMax);
            
            % Evaluation
            F(:,i) = Custo(Ne,N_pesos,Theta(:,:,i),u,y,n_transienth,Ts);
            
        end
        
        for i = 1:P
            for l = 1:VarSize(1)
                % Update Personal Best
%                 tempFp = Custo(Ne,N_pesos,Thetap(:,:,i),u,y,n_transienth,Ts);
                if  F(l,i) < Fp(l,i)
                    Thetap(l,:,i) = Theta(l,:,i);
                    Fp(l,i) = F(l,i);
                end
                
                % Update Global Best
                if  F(l,i) < Fg(l) 
                    if  Fg(l)-F(l,i) > 0.1
                        mark = 0;
                    end
                    Thetag(l,:) = Theta(l,:,i);
                    Fg(l) = F(l,i);
%                     disp(['Iteration ' num2str(k) ': Best Cost ' num2str(l)  ' = ' num2str(Fg(l))]);
                end
            end
        end
        rho = rand(VarSize(1),1);
        if mark > 10 % Mutation takes effect if after 10 iterations the
            if rand < alpha
                beta = randi(P);
                gamma = randi(Nv);
                for l = 1:VarSize(1)
                    if rho(l) <= 0.5
                        Theta(l,gamma,beta) = Theta(l,gamma,beta) - rho(l)*VelMax;
                    else
                        Theta(l,gamma,beta) = Theta(l,gamma,beta) + rho(l)*VelMax;
                    end
                end
            end
        end
%         plot(k,Fg(1),'ro')
    end


    %%
    save_Theta(:,:,realiza) = Thetag;
    
    vtd = Histerese(length(utd),N_pesos,Thetag,utd,Ts,2);
    ytdr = ytd + 0*max(ytd,[],2)/(1.3*FR) .* randn(size(ytd));
    vte = Histerese(length(u),N_pesos,Thetag,u,Ts,2);
%     figure
%     plot(y','k-')
%     hold on
%     plot(vte','r--')
%     rmseestat = sqrt(sum((vte(:,7:end) - y(:,7:end)).^2,2)./length(y(:,7:end)))./max(y(:,7:end),[],2)
%     figure
%     plot(vtd_o(1,500:700),'k-')
%     hold on
%     plot(vtd(1,500:700)*max(vtd_o(1,500:700))/max(vtd(1,500:700)),'r--')
%     rmsedin = sqrt(sum((vtd_o(:,500:700) - vtd(:,500:700)).^2,2)./200)./max(vtd(:,500:700),[],2)
%     plot(ytd(1,500:700)*max(vtd_o(1,500:700))/max(ytd(1,500:700)),'b--')
%     [A, B, C, D, ~] = moesp(vtd_o(:,n_transiente:end).*(max(vtd(:,n_transiente:end),...
%         [],2)./max(vtd_o(:,n_transiente:end),[],2)),ytdr(:,n_transiente:end),2,3);
      [A, B, C, D, SS] = moesp(vtd(:,n_transiente:end),ytdr(:,n_transiente:end),2,2);
%     G0 = C*((eye(2)-A)\B) + D
%     sys = ss(A, B, C, D, 1);
%     pole(sys)
%     tzero(sys)
%     vo = Processo_v(length(utd),utd,Ts,2);
%     [A, B, C, D, ~] = moesp(vo(:,n_transiente:end)*max(vtd(1,n_transiente:end))/max(vo(1,n_transiente:end)),ytdr(:,n_transiente:end),2,4);

    autovalor(:,realiza) = eig(A);
    %% Validação
    y_r = dados_nuvem(length(u_test),N_pesos,Thetag,u_test,A,B,C,D,Ts);
    y_hold_all(realiza,:,:) = y_r;
    
    y_est = dados_nuvem(length(uvd),N_pesos,Thetag,uvd,A,B,C,D,Ts);
    
%     figure
% plot(y_est')
% hold on
% plot(yvd','k--')
% figure
% plot(utd(1,1:200))
    ini = 126; %Selecione esse ponto com base no máximo ou mínimo de uma das senoides da entrada
    fim = ini+100;
    aux = y_est(1,:);
    y1s = max([aux ; y1s]);
    y1i = min([aux ; y1i]);
    aux = y_est(2,:);
    y2s = max([aux ; y2s]);
    y2i = min([aux ; y2i]);

    erro(realiza,:,:) = yvd(:, ini:fim) - y_est(:, ini:fim);
    
    rmseval = sqrt(sum((y_est - yvd).^2,2)./length(yvd))./max(yvd,[],2);
    save_rmse(:, realiza) = rmseval;
    srmse = srmse + rmseval;
    tp = toc;
    st = tp + st;
    tr = (R-realiza)*st/realiza;
    tr = seconds(tr);
    tr.Format = 'hh:mm:ss';
    fprintf('faltam aprox %s\n',tr)
end
% save exemplo1_1000_sim
%%
% load exemplo1_1000_sim
figure
subplot(1,2,1)
histogram(save_rmse(1,:), 10, 'FaceColor', [0 0 1], 'FaceAlpha', 0.5)
hold on
xline(mean(save_rmse(1,:)), 'b--', 'LineWidth', 1.5)
xlabel('NRMSE', 'FontName', 'Times', 'FontSize', 14)
xlim([0 0.05])
ylabel('count', 'FontName', 'Times', 'FontSize', 14)
box on
legend('NRMSE$(y_1)$', 'mean NRMSE$(y_1)$', 'Interpreter', 'latex')
set(gca,'FontName','Times','FontSize',14)
subplot(1,2,2)
hold on
histogram(save_rmse(2,:), 10, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.5)
xline(mean(save_rmse(2,:)), '--','color',[.2 .39 .2 1], 'LineWidth', 1.5)
xlim([0 0.05])
xlabel('NRMSE', 'FontName', 'Times', 'FontSize', 14)
ylabel('count', 'FontName', 'Times', 'FontSize', 14)
box on
legend('NRMSE$(y_2)$', 'mean NRMSE$(y_2)$', 'Interpreter', 'latex')
set(gca,'FontName','Times','FontSize',14)
set(gcf, 'Color', 'w', 'Position', [50 200 800 200])
%
disp(srmse/R)



figure
subplot(1,2,1)
hold on
for realiza = 1:R
    % O squeeze transforma y(realiza, 1, :) de 1x1xN para um vetor Nx1
    y1 = squeeze(y_hold_all(realiza, 1, :));
    y2 = squeeze(y_hold_all(realiza, 2, :));
    if realiza == 1
        h_y1 = plot((1:length(y1))*Ts,y1, 'Color', [.55 .64 .82 1]); % Azul claro
        h_y2 = plot((1:length(y2))*Ts,y2, 'Color', [.47 .75 .45 1]); % Verde claro
    else
        plot((1:length(y1))*Ts,y1, 'Color', [.55 .64 .82 1]); 
        plot((1:length(y2))*Ts,y2, 'Color', [.47 .75 .45 1]); 
    end
end
h_proc_y1 = plot((1:length(y_proc_test))*Ts,y_proc_test(1,:), 'k-');
h_proc_y2 = plot((1:length(y_proc_test))*Ts,y_proc_test(2,:), 'r-');
% Calcula a média ao longo de R (dimensão 1) e "achata" o resultado
media_y1 = squeeze(mean(y_hold_all(:, 1, :), 1)); 
media_y2 = squeeze(mean(y_hold_all(:, 2, :), 1));
% Plota as médias (ex: azul tracejado para y1, vermelho tracejado para y2)
h_med1 = plot((1:length(media_y1))*Ts,media_y1, 'b--','linewidth',1.5); 
h_med2 = plot((1:length(media_y2))*Ts,media_y2, '--','color',[.2 .39 .2 1],'linewidth',1.5);
xlabel('$kT_s$(s)', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 14)
ylabel('$y(kT_s)$', 'Interpreter', 'latex', 'FontSize', 14)
legend([h_y1, h_y2, h_proc_y1, h_proc_y2, h_med1, h_med2], ...
       {'$\hat{y}_1(kT_s)$', '$\hat{y}_2(kT_s)$', '$y_1(kT_s)$',...
       '$y_2(kT_s)$', 'mean $\hat{y}_1$', 'mean $\hat{y}_2$'}, ...
       'Interpreter', 'latex');
xlim([0 260]*Ts)
ylim([-15 15])
box on
set(gca,'FontName','Times','FontSize',14)

subplot(1,2,2)
hold on
for realiza = 1:R
    % O squeeze transforma y(realiza, 1, :) de 1x1xN para um vetor Nx1
    y1 = squeeze(y_hold_all(realiza, 1, :));
    y2 = squeeze(y_hold_all(realiza, 2, :));
    if realiza == 1
        h_y1 = plot((1:length(y1))*Ts,y1, 'Color', [.55 .64 .82 1]); % Azul claro
        h_y2 = plot((1:length(y2))*Ts,y2, 'Color', [.47 .75 .45 1]); % Verde claro
    else
        plot((1:length(y1))*Ts,y1, 'Color', [.55 .64 .82 1]); 
        plot((1:length(y2))*Ts,y2, 'Color', [.47 .75 .45 1]); 
    end
end
h_proc_y1 = plot((1:length(y_proc_test))*Ts,y_proc_test(1,:), 'k-');
h_proc_y2 = plot((1:length(y_proc_test))*Ts,y_proc_test(2,:), 'r-');
% Calcula a média ao longo de R (dimensão 1) e "achata" o resultado
media_y1 = squeeze(mean(y_hold_all(:, 1, :), 1)); 
media_y2 = squeeze(mean(y_hold_all(:, 2, :), 1));
% Plota as médias (ex: azul tracejado para y1, vermelho tracejado para y2)
h_med1 = plot((1:length(media_y1))*Ts,media_y1, 'b--','linewidth',1.5); 
h_med2 = plot((1:length(media_y2))*Ts,media_y2, '--','color',[.2 .39 .2 1],'linewidth',1.5);
xlabel('$kT_s$(s)', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 14)
ylabel('$y(kT_s)$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([160 260]*Ts)
ylim([5 15])
box on
set(gca,'FontName','Times','FontSize',14)
set(gcf, 'Color', 'w', 'Position', [50 200 800 350])


figure
zplane([],[]);
hold on
th = 0:pi/50:2*pi;
plot(cos(th), sin(th), 'k:')
axis equal
% b = 1;
xlim([-1 1])
ylim([-1 1])

cores = lines(size(autovalor, 1));  % gera N cores distintas

hold on;
for a = 1:size(autovalor, 1)
    plot(real(autovalor(a,:)), imag(autovalor(a,:)), 'x', ...
         'Color', cores(a,:), 'MarkerSize', 20);
end


% for a = 1:size(autovalor,1)
%     hold on
%     plot(real(autovalor(a,:)),imag(autovalor(a,:)),'x','color',[a-1 0 b],'markersize',20);
%     b = b-1;
% end
h = legend({'','','','','Polo 1', 'Polo 2'});
axis equal
xlim([-1 1])
ylim([-1 1])
xlabel('parte real')
ylabel('parte imaginária')
set(gca,'FontName','Times','FontSize',14)
set(gcf,'Color','w','Position',[50 200 500 350])

theta = [0.1, 0.9, 1, 0.1, 0.9, 0.3, 0.2, 0.6, 0.8, 0.4, 0.9, 0.9, 0.35, 1, 0.25;
         0.1, 1, 0.9, 0.1, 0.7, 0.1, 0.8, 0.4, 0.1, 0.5, 0.6,  1,  0.4,  1, 0.4];

for a = 1:size(save_Theta,1)
    figure
    plot(theta(a,:),'*','color',[0 0 0],'markersize',8)
    hold on
    plot(reshape(save_Theta(a,:,:),[VarSize(2),R]),'*','color',[1 0 0],'markersize',8)
    legend('processo','estimado')
    xlabel('posição da variável em \Theta')
    ylabel('valor')
    xlim([0 VarSize(2)+1])
    set(gca,'FontName','Times','FontSize',14)
    set(gcf,'Color','w','Position',[50 200 500 350])
end

%

% s - superior
% i - inferior
% d - dentro
% f - fora
qd = 1/f1d/Ts;
u1f = uvd(1,ini:ini+qd-1);
u1d = u1f;
u2f = uvd(2,ini:ini+qd-1);
u2d = u2f;


% --- CÓDIGO CORRIGIDO ---

% Índices para a segunda metade do ciclo
idx_start = ceil(qd/2) + 1;
idx_end = qd;

% Índices correspondentes nos vetores originais (y1s, y1i, etc.)
source_idx_start = ini + ceil(qd/2);
source_idx_end = ini + qd - 1;

% Construção dos laços para y1
y1d = y1s(ini:ini+qd-1);
y1d(idx_start:idx_end) = y1i(source_idx_start:source_idx_end); % Correto

y1f = y1i(ini:ini+qd-1);
y1f(idx_start:idx_end) = y1s(source_idx_start:source_idx_end); % Correto

% Construção dos laços para y2
y2d = y2s(ini:ini+qd-1);
y2d(idx_start:idx_end) = y2i(source_idx_start:source_idx_end); % Correto

y2f = y2i(ini:ini+qd-1);
y2f(idx_start:idx_end) = y2s(source_idx_start:source_idx_end); % Correto



% fecha todos caminhos
u1f = [u1f u1f(1)];
y1f = [y1f y1f(1)];
u1d = [u1d u1d(1)];
y1d = [y1d y1d(1)];
u2f = [u2f u2f(1)];
y2f = [y2f y2f(1)];
u2d = [u2d u2d(1)];
y2d = [y2d y2d(1)];

% Identifica os pontos onde tem cruzamento que ficam no
cruz1 = 0;
cruz2 = 0;
pcruz11 = 0; pcruz12 = 0; pcruz21 = 0; pcruz22 = 0;
for k = 1:ceil(qd/2)
    if y1d(k) < y1d(qd+2-k) && cruz1 == 0
        pcruz11 = k;
        cruz1 = 1;
    elseif y1d(k) > y1d(qd+2-k) && cruz1 == 1
        pcruz12 = k;
        cruz1 = 2;
    end
    if y2d(k) < y2d(qd+2-k) && cruz2 == 0
        pcruz21 = k;
        cruz2 = 1;
    elseif y2d(k) > y2d(qd+2-k) && cruz2 == 1
        pcruz22 = k;
        cruz2 = 2;
    end
end
% fecha o circulo de dentro nos pontos de cruzamento

if pcruz11 ~= 0
    u1d(1:pcruz11) = u1d(pcruz11);
    u1d(end-pcruz11+1:end) = u1d(pcruz11);
    y1d(1:pcruz11) = y1d(pcruz11);
    y1d(end-pcruz11+1:end) = y1d(pcruz11);
end
if pcruz12 ~= 0
    u1d(pcruz12:pcruz12+2*((ceil(qd/2)-pcruz12))+1) = u1d(pcruz12-1);
    y1d(pcruz12:pcruz12+2*((ceil(qd/2)-pcruz12))+1) = y1d(pcruz12-1);
end
% u2s = u2;
if pcruz21 ~= 0
    u2d(1:pcruz21) = u2d(pcruz21);
    u2d(end-pcruz21+1:end) = u2d(pcruz21);
    y2d(1:pcruz21) = y2d(pcruz21);
    y2d(end-pcruz21+1:end) = y2d(pcruz21);
end
if pcruz22 ~= 0
    u2d(pcruz22:pcruz22+2*((ceil(qd/2)-pcruz22))+1) = u2d(pcruz22-1);
    y2d(pcruz22:pcruz22+2*((ceil(qd/2)-pcruz22))+1) = y2d(pcruz22-1);
end

%

figure
patch([u1d fliplr(u1f)],[y1d fliplr(y1f)],[0, 0, 1],...
    'facealpha',0.5,'edgealpha',0,'EdgeColor',[0, 0, 1]);
hold on
patch([u2d fliplr(u2f)],[y2d fliplr(y2f)],[0 0.6 0],...
    'facealpha',0.5,'edgealpha',0,'EdgeColor',[0 0.6 0]);

plot(uvd(1,n_transientv:n_transientv+qd),yvd(1,n_transientv:n_transientv+qd),'k-')
plot(uvd(2,n_transientv:n_transientv+qd),yvd(2,n_transientv:n_transientv+qd),'r-')
yticks(-200:50:200)
xticks(-5:5)
box on
xlabel('$u(k)$','interpreter','latex')
ylabel('$y(k)$','interpreter','latex')
legend({'$\hat{y}_1(k)$','$\hat{y}_2(k)$','$y_1(k)$','$y_2(k)$'},'interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 290 500 350])
ylim([-15 15])

t = 0:Ts:length(yvd)*Ts - Ts;
figure
patch([t fliplr(t)],[y1s fliplr(y1i)],[0, 0, 1],...
    'facealpha',0.5,'edgealpha',0,'EdgeColor',[0, 0, 1]);
hold on
patch([t fliplr(t)],[y2s fliplr(y2i)],[0 0.6 0],...
    'facealpha',0.5,'edgealpha',0,'EdgeColor',[0 0.6 0]);
hold on
plot(t,yvd(1,:),'k-');
plot(t,yvd(2,:),'r-');
box on
xlim([0 t(end)])
xlabel('$kT_s(s)$','interpreter','latex')
ylabel('$y(kT_s)$','interpreter','latex')
legend({'$\hat{y}_1(kT_s)$','$\hat{y}_2(kT_s)$','$y_1(kT_s)$','$y_2(kT_s)$'},'interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
ylim([-15 15])
set(gcf,'Color','w','Position',[50 290 500 350])



erro1 = squeeze(erro(:,1,:));
erro2 = squeeze(erro(:,2,:));
u_ciclo = uvd(1, ini:fim);

% Separar ramo ascendente e descendente
idx_asc = [true; diff(u_ciclo') > 0];  % simplificado
idx_desc = ~idx_asc;
%
figure

% --- SUBPLOT 1 ---
subplot(2,1,1)
hold on

% Ramo Ascendente (Erro 1)
x_asc = u_ciclo(idx_asc);
y_asc1 = mean(erro1(:,idx_asc), 1);
% Preenchimento azul com 60% de transparência (Alpha = 0.4)
a1 = area(x_asc, y_asc1, 'FaceColor', [0 0 1], 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'LineWidth', 1);

% Ramo Descendente (Erro 1) - Ordenando X para evitar problemas no preenchimento
x_desc = u_ciclo(idx_desc);
y_desc1 = mean(erro1(:,idx_desc), 1);
[x_desc_ord, idx_sort] = sort(x_desc); % Ordena o eixo X
y_desc1_ord = y_desc1(idx_sort);       % Reorganiza o Y de acordo

% Preenchimento vermelho com 60% de transparência
a2 = area(x_desc_ord, y_desc1_ord, 'FaceColor', [1 0 0], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 1);
xticks(-5:2:5)
xlabel('$u_1(k)$', 'Interpreter', 'latex')
ylabel('$y_1(k)-\hat{y}_1(k)$', 'Interpreter', 'latex')
legend([a1, a2], {'$u_1(k)>u_1(k-1)$', '$u_1(k)<u_1(k-1)$'}, 'interpreter','latex')
set(gca, 'FontName', 'Times', 'FontSize', 14)
box on


% --- SUBPLOT 2 ---
subplot(2,1,2)
hold on

% Ramo Ascendente (Erro 2)
y_asc2 = mean(erro2(:,idx_asc), 1);
a1 = area(x_asc, y_asc2, 'FaceColor', [0 0 1], 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'LineWidth', 1);

% Ramo Descendente (Erro 2)
y_desc2 = mean(erro2(:,idx_desc), 1);
y_desc2_ord = y_desc2(idx_sort); % Usa a mesma ordenação de antes
a2 = area(x_desc_ord, y_desc2_ord, 'FaceColor', [1 0 0], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 1);
xticks(-5:2:5)
xlabel('$u_2(k)$', 'Interpreter', 'latex')
ylabel('$y_2(k)-\hat{y}_2(k)$', 'Interpreter', 'latex')
legend([a1, a2], {'$u_2(k)>u_2(k-1)$', '$u_2(k)<u_2(k-1)$'}, 'interpreter','latex')
set(gca, 'FontName', 'Times', 'FontSize', 14)
box on

% Configuração da Janela
set(gcf, 'Color', 'w', 'Position', [50 200 500 400])

% save exemplo1_1000_sim
% system('shutdown -s')
%% Funções


function v_est = Histerese(N_dados,N_pesos,pesos,ute,Ts,tipo)

theta = pesos(:,1:N_pesos+1);
theta_hat = pesos(:,N_pesos+2:end);
xi = ones(size(pesos(:,end)));
beta = 0*pesos(:,end);
v_est = zeros(size(ute,1),N_dados);
Ph = zeros(size(ute,1),N_pesos);
Bh = zeros(size(ute,1),1);
for k = 1:N_dados
    BL = max(ute,[],2).*theta_hat(:,1).*tanh(theta_hat(:,2).*ute(:,k));
%     ((theta_hat(:,1).*tanh(theta_hat(:,2).*ute(:,k)))/(theta_hat(:,3).*tanh(theta_hat(:,4).*ute(:,k))))*
    BR = max(ute,[],2).*theta_hat(:,3).*tanh(theta_hat(:,4).*ute(:,k));
    if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
        delta = zeros(size(ute,1),1);
        Ph = zeros(size(ute,1),N_pesos);
        Bh = zeros(size(ute,1),1);
    else
        delta = ute(:,k) - ute(:,k-1);
    end
    for l = 1:size(ute,1)
        if delta(l) >= 0
            Bh(l) = BR(l);
        else
            Bh(l) = BL(l);
        end
    end
    v_est(:,k) = theta(:,1).*Bh;
    for j = 1:N_pesos
        rj = xi.*vecnorm(ute,Inf,2)*(j-1)/(N_pesos);
        if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
            v_est(:,k) = theta(:,1).*Bh;
            epsilon1 = BL + rj;
            epsilon2 = BR - rj;
            for l = 1:size(ute,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v_est(:,k) = v_est(:,k) + theta(:,j+1).*Ph(:,j);
        else
            epsilon1 = BL + rj - Ph(:,j);
            epsilon2 = BR - rj - Ph(:,j);
            for l = 1:size(ute,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v_est(:,k) = v_est(:,k) + theta(:,j+1).*Ph(:,j);
        end
    end
end
end






function F = Custo(N_dados,N_pesos,pesos,ute,vte,n_transient,Ts)

% N_dados = ceil(N_dados/size(ute,1));
% ute =
theta = pesos(:,1:N_pesos+1);
theta_hat = pesos(:,N_pesos+2:end);
xi = ones(size(pesos(:,end)));
beta = 0*pesos(:,end);
v_est = zeros(size(ute,1),N_dados);
Ph = zeros(size(ute,1),N_pesos);
Bh = zeros(size(ute,1),1);
for k = 1:N_dados
    BL = max(ute,[],2).*theta_hat(:,1).*tanh(theta_hat(:,2).*ute(:,k));
    BR = max(ute,[],2).*theta_hat(:,3).*tanh(theta_hat(:,4).*ute(:,k));
    if k == 1
        delta = zeros(size(ute,1),1);
        Ph = zeros(size(ute,1),N_pesos);
        Bh = zeros(size(ute,1),1);
    else
        delta = ute(:,k) - ute(:,k-1);
    end
    for l = 1:size(ute,1)
        if delta(l) >= 0
            Bh(l) = BR(l);
        else
            Bh(l) = BL(l);
        end
    end
    v_est(:,k) = theta(:,1).*Bh;
    for j = 1:N_pesos
        rj = xi.*vecnorm(ute,Inf,2)*(j-1)/(N_pesos);
        if k == 1
            v_est(:,k) = theta(:,1).*Bh;
            epsilon1 = BL + rj;
            epsilon2 = BR - rj;
            for l = 1:size(ute,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v_est(:,k) = v_est(:,k) + theta(:,j+1).*Ph(:,j);
        else
            epsilon1 = BL + rj - Ph(:,j);
            epsilon2 = BR - rj - Ph(:,j);
            for l = 1:size(ute,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v_est(:,k) = v_est(:,k) + theta(:,j+1).*Ph(:,j);
        end
    end
end
v_est = v_est(:,n_transient+1:end-1);
vte = vte(:,n_transient+1:end-1);
F = sqrt(sum((v_est - vte).^2,2)/N_dados)./max(vte,[],2);
end







function y = Processo(N_dados,u,Ts,tipo)
% Parametros do processo de    eig      0.8010 +- 0.0115i
% MB  
% A = [0.2105   -0.1654;
%     3.5739    1.5037];
% 
% B =[0.8825    0.8825;
%     0.5107    0.5107];
% 
% C =[0.6023    0.6023;
%     0.1636    0.1636];

% Melhor
A = [0.9   0;
    0   0.7];

B = [0.4    0;
    0    1];

C = [0.5    0;
    0    0.6];


% A = [1.6   -0.4;
%     1.7    0.1];
% 
% B = 0.1*[0.9    0.9;
%     0.1    0.1];
% 
% C = [0.7    0.7;
%     0.4    0.4];


% A 0.0240    1.8873
%    -0.3311    1.6028
% 
% B 0.0185    0.0185
%     0.0858    0.0858
%    
% C 0.0717    0.0717
%     0.0539    0.0539
    
% A = [0.28    1.9
%    -0.15    1.33];
% 
% B = [0.01    0.01
%     0.08    0.08];
% 
% C = [0.07    0.07
%     0.05    0.05];
% 
% D = [0 0; 0 0];

% G0 = C*((eye(2)-A)\B) + D
% sys = ss(A, B, C, D, 1);
% pole(sys)
% tzero(sys)

% MB
% A = [1.3142   -0.1121;
%     2.7311    0.4580];
% 
% B = [0.5511    0.5511;
%     0.1465    0.1465];
% 
% C = [0.4041    0.4041;
%     0.1128    0.1128];

% ficou mais flat
% A = [0.3403   -0.0941;
%     3.0038    1.2915];
% 
% B = [0.2717    0.2717;
%     0.1150    0.1150];
% 
% C = [0.1680    0.1680;
%     0.2859    0.2859];





% A = [0.4    1.7;
%     -0.15    1.3];
% 
% B = [0.1    0.1;
%     1    1];
% 
% C = [0.1    0.1;
%     0.2    0.2]/6;
    
D = [0 0; 0 0];

theta = [0.1, 0.9, 0.4, 0.8;
         0  , 0.7, 0.2, 0.5];
theta_hat = [0.9, 0.6, 0.9, 0.6;
             0.9, 0.6, 0.9, 0.6];

xi = [1; 1];
beta = [1; 1];
N_pesos = size(theta,2)-1;
y = zeros(size(C,1),N_dados);
x = zeros(size(A,1),N_dados);
Ph = zeros(size(B,2),N_pesos);
Bh = zeros(size(B,2),1);
for k = 1:N_dados
    BL = max(u,[],2).*theta_hat(:,1).*tanh(theta_hat(:,2).*u(:,k));
    BR = max(u,[],2).*theta_hat(:,3).*tanh(theta_hat(:,4).*u(:,k));
    if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
        delta = zeros(size(u,1),1);
        Ph = zeros(size(u,1),N_pesos);
        Bh = zeros(size(u,1),1);
    else
        delta = u(:,k) - u(:,k-1);
    end
    for l = 1:size(B,2)
        if delta(l) >= 0
            Bh(l) = BR(l);
        else
            Bh(l) = BL(l);
        end
    end
    v = theta(:,1).*Bh;
    for j = 1:N_pesos
        rj = xi.*vecnorm(u,Inf,2)*(j-1)/(N_pesos);
        if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
            v = theta(:,1).*Bh;
            x(:,k) = zeros(size(x(:,k)));
            epsilon1 = BL + rj;
            epsilon2 = BR - rj;
            for l = 1:size(u,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v = v + theta(:,j+1).*Ph(:,j);
        else
            epsilon1 = BL + rj - Ph(:,j);
            epsilon2 = BR - rj - Ph(:,j);
            for l = 1:size(u,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v = v + theta(:,j+1).*Ph(:,j);
        end
    end
    x(:,k+1) = A*x(:,k) + B*v;
    y(:,k) = C*x(:,k) + D*v;
end
end


function vo = Processo_v(N_dados,u,Ts,tipo)
% Parametros do processo

% A = [0.2105   -0.1654;
%     3.5739    1.5037];
% 
% B =[0.8825    0.8825;
%     0.5107    0.5107];
% 
% C =[0.6023    0.6023;
%     0.1636    0.1636];



A = [0.9   0;
    0   0.7];

B = [0.4    0;
    0    1];

C = [0.5    0;
    0    0.6];


% A = [1.6   -0.4;
%     1.7    0.1];
% 
% B = 0.1*[0.9    0.9;
%     0.1    0.1];
% 
% C = [0.7    0.7;
%     0.4    0.4];

% A = [0.28    1.9
%    -0.15    1.33];
% 
% B = [0.01    0.01
%     0.08    0.08];
% 
% C = [0.07    0.07
%     0.05    0.05];


% A = [1.3142   -0.1121;
%     2.7311    0.4580];
% 
% B = [0.5511    0.5511;
%     0.1465    0.1465];
% 
% C = [0.4041    0.4041;
%     0.1128    0.1128];
 
% A = [0.3403   -0.0941;
%     3.0038    1.2915];
% 
% B = [0.2717    0.2717;
%     0.1150    0.1150];
% 
% C = [0.1680    0.1680;
%     0.2859    0.2859];




% A = [0.4    1.7;
%     -0.15    1.3];
% 
% B = [0.1    0.1;
%     1    1];
% 
% C = [0.1    0.1;
%     0.2    0.2]/6;



D = [0 0; 0 0];

theta = [0.1, 0.9, 0.4, 0.8;
         0  , 0.7, 0.2, 0.5];
theta_hat = [0.9, 0.6, 0.9, 0.6;
             0.9, 0.6, 0.9, 0.6];

xi = [1; 1];
beta = [1; 1];
N_pesos = size(theta,2)-1;
vo = zeros(size(C,1),N_dados);
x = zeros(size(A,1),N_dados);
Ph = zeros(size(B,2),N_pesos);
Bh = zeros(size(B,2),1);
for k = 1:N_dados
    BL = max(u,[],2).*theta_hat(:,1).*tanh(theta_hat(:,2).*u(:,k));
    BR = max(u,[],2).*theta_hat(:,3).*tanh(theta_hat(:,4).*u(:,k));
    if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
        delta = zeros(size(u,1),1);
        Ph = zeros(size(u,1),N_pesos);
        Bh = zeros(size(u,1),1);
    else
        delta = u(:,k) - u(:,k-1);
    end
    for l = 1:size(B,2)
        if delta(l) >= 0
            Bh(l) = BR(l);
        else
            Bh(l) = BL(l);
        end
    end
    v = theta(:,1).*Bh;
    for j = 1:N_pesos
        rj = xi.*vecnorm(u,Inf,2)*(j-1)/(N_pesos);
        if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
            v = theta(:,1).*Bh;
            x(:,k) = zeros(size(x(:,k)));
            epsilon1 = BL + rj;
            epsilon2 = BR - rj;
            for l = 1:size(u,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v = v + theta(:,j+1).*Ph(:,j);
        else
            epsilon1 = BL + rj - Ph(:,j);
            epsilon2 = BR - rj - Ph(:,j);
            for l = 1:size(u,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v = v + theta(:,j+1).*Ph(:,j);
        end
    end
    vo(:,k) = v;
end
end


function y_est = dados_nuvem(N_dados,N_pesos,pesos,uvd,A,B,C,D,Ts)

theta = pesos(:,1:N_pesos+1);
theta_hat = pesos(:,N_pesos+2:end);
xi = ones(size(pesos(:,end)));
beta = 0*pesos(:,end);
y_est = zeros(size(C,1),N_dados);
x = zeros(size(A,1),N_dados);
Ph = zeros(size(B,2),N_pesos);
Bh = zeros(size(B,2),1);
for k = 1:N_dados
    BL = max(uvd,[],2).*theta_hat(:,1).*tanh(theta_hat(:,2).*uvd(:,k));
    BR = max(uvd,[],2).*theta_hat(:,3).*tanh(theta_hat(:,4).*uvd(:,k));
    if k == 1
        delta = zeros(size(uvd,1),1);
        Ph = zeros(size(B,2),N_pesos);
        Bh = zeros(size(B,2),1);
    else
        delta = uvd(:,k) - uvd(:,k-1);
    end
    for l = 1:size(B,2)
        if delta(l) >= 0
            Bh(l) = BR(l);
        else
            Bh(l) = BL(l);
        end
    end
    v = theta(:,1).*Bh;
    for j = 1:N_pesos
        rj = xi.*vecnorm(uvd,Inf,2)*(j-1)/(N_pesos);
        if k == 1
            x = zeros(size(A,1),N_dados);
            v = theta(:,1).*Bh;
            epsilon1 = BL + rj;
            epsilon2 = BR - rj;
            for l = 1:size(uvd,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v = v + theta(:,j+1).*Ph(:,j);
        else
            epsilon1 = BL + rj - Ph(:,j);
            epsilon2 = BR - rj - Ph(:,j);
            for l = 1:size(uvd,1)
                if delta(l) < 0 && epsilon1(l) < 0
                    Ph(l,j) = BL(l) + rj(l);
                elseif delta(l) > 0 && epsilon2(l) > 0
                    Ph(l,j) = BR(l) - rj(l);
                end
            end
            v = v + theta(:,j+1).*Ph(:,j);
        end
    end
    x(:,k+1) = A*x(:,k) + B*v;
    y_est(:,k) = C*x(:,k) + D*v;
end
end






function [A,B,C,D,SS] = moesp(u,y,n,k)
% Sintaxe
% [A,B,C,D,SS] = moesp(u,y,n,k);
% Descrição
% O algoritmo moesp estima as matrizes de modelos na representação espaço de estados sem
% tratar ruído de processo e/ou de medição. Também pode ser utilizado para determinar a ordem
% do sistema.
% Dados de Entrada
% u -> matriz de dados medidos das entradas do sistema (linha);
% y -> matriz de dados medidos das saídas do sistema (linha);
% n -> ordem escolhida para o modelo do sistema;
% k -> número de linhas da matriz em blocos de Hankel, este número é definido como
% k = 2(maxord/nusaida). Em que, maxord é a máxima ordem que o usuário pressupõe que o
% sistema tenha e nusaida é o número de saídas que o sistema possuí. Esta relação é definida de
% forma empírica por (Van Overschee e De Moor, 1996).
% Dados de Saída
% A -> matriz dinâmica estimada para o sistema em espaço de estados;
% B -> matriz de entrada estimada para o sistema em espaço de estados;
% C -> matriz de saída estimada para o sistema em espaço de estados;
% D -> matriz de transmissão direta estimada para o sistema em espaço de estados;
% SS -> matriz de valores singulares da projeção oblíqua, utilizada para determinar a ordem
% do sistema;
% Implementado por Rodrigo Augusto Ricco: rodrigo.ricco@yahoo.com.br
% Universidade Federal de Minas Gerais - UFMG
% Última modificação 14/11/2012
% m=dim(u), l=dim(y), n=dim(x)

i = k;% i=numero de linhas;
% U=2im x j ; Y=2il x j
[l, ~] = size(y);
[m, Ndados] = size(u);
j = Ndados - 2*i;
kk = 0;
% Montando as matrizes em blocos de Hankel
for k = 1:m:2*i*m-m+1
    kk = kk+1;
    U(k:k+m-1,:) = u(:,kk:kk+j-1); % Matriz de dados U
end
kk = 0;
for k = 1:l:2*i*l-l+1
    kk = kk+1;
    Y(k:k+l-1,:) = y(:,kk:kk+j-1); % Matriz de dados Y
end
Uf = U(i*m+1:2*i*m,:); % entradas futuras
Yf = Y(i*l+1:2*i*l,:); % saídas futuras
%Up = U(1:i*m,:); % entradas passadas
%Yp = Y(1:i*l,:); % saídas passadas
%Wp = [Up; Yp]; % empilhando Up e Yp ---- Não usada nesse desenvolvimento
H = [Uf; Yf]; % matriz de dados final
% Decomposição LQ
L = triu(qr(H'))'; % ***********************
im=i*m;
il=i*l;
L11 = L(1:im,1:im);
L21 = L(im+1:im+il,1:im);
L22 = L(im+1:im+il,im+1:im+il);
Oi = L22;
% Decomposição em valores singulares
[UU,SS,~] = svd(Oi); % X = UU*SS*VV' é a forma como o sdv desmembra a matriz.
% SS tem os valores singulares em ordem decrescente
U1 = UU(:,1:n);
Ob_est_i = U1*sqrtm(SS(1:n,1:n));

C = Ob_est_i(1:l,1:n);
A = pinv(Ob_est_i(1:l*(i-1),1:n))*Ob_est_i(l+1:l*i,1:n);

U2 = UU(:,n+1:size(UU',1))';
Z = U2*L21/L11;
% A partir deste ponto todos os algortimos são iguais
mm = l*i-n;
M = zeros(mm*i,m);
LL = zeros(mm*i,l+n);
for h = 1:i
    M((h-1)*mm+1:h*mm,:) = Z(:,(h-1)*m+1:h*m);
    LL((h-1)*mm+1:h*mm,:) = [U2(:,(h-1)*l+1:h*l) U2(:,h*l+1:end)*Ob_est_i(1:end-h*l,:)];
end

DB = pinv(LL)*M;
D = DB(1:l,:);
B = DB(l+1:size(DB,1),:);

end


