% Versão usada no exemplo 1 do artigo ressubmetido 2026
%% Ensaio quase estatico que gera dados para aplicar no MPSO

clc
clear
close all

nu = 2; % numero de entradas
FR = 24;              % Fator de SNR para o ruído
R = 2;                % Qtd execuções nuvem
N_pesos = 10;         % Numero de disparos
Nv = N_pesos+7;       % Numero total de parametros
VarSize = [2 Nv];     % Size of Decision Variables Matrix
VarMin = 0;          % 0 Lower Bound of Variables
VarMax = 4;          % 1 Upper Bound of Variables
% PSO Parameters
b = 4;
alpha = 0.7;          % Mutation probability
Kmax = 200;           % 1000 Maximum Number of Iterations
P = 60;               % Population Size (Swarm Size)

Ts = 0.001;

%% identificar dinâmica

num_bloco_hankel = 2;
Amp = 1.82;
f1td = 20;
f2td = 30;
f3td = 40;
Nd = ceil(40/min(f1td,f2td)/Ts);
t = 0:Ts:Nd*Ts - Ts;
n_transiente = ceil(1/min(f1td,f2td)/Ts);

utd = [Amp*sin(2*pi*f1td.*t)+Amp*sin(2*pi*f2td.*t)+Amp*sin(2*pi*f3td.*t);
       Amp*sin(2*pi*f1td.*t)+Amp*sin(2*pi*f2td.*t)+Amp*sin(2*pi*f3td.*t)];
ytd = Processo(length(utd),utd,Ts,2);

%% Validação da dinâmica

f1d = 10;
N_val = ceil(3/f1d/Ts);

Amp = 5;
t = 0:Ts:N_val*Ts - Ts;

uvd = [Amp*sin(2*pi*f1d.*t + pi); Amp*sin(2*pi*f1d.*t + pi)];
yvd = Processo(N_val,uvd,Ts,2);
n_transientv = 1/f1d/Ts;


%% Identificação da histerese
L = 200;
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


%% PSO Parameters
phi1 = 2.05;
phi2 = 2.05;
phi = phi1+phi2;
chi = 2/abs(2-phi-sqrt(phi^2-4*phi));
c1 = phi1;       % Personal Learning Coefficient
c2 = phi2;       % Global Learning Coefficient
% Velocity Limits
VelMax = 0.2;
%% Inicializa os parâmetros da nuvem
y1s = -inf(1,length(yvd));
y2s = y1s;
y1i = -y1s;
y2i = -y1s;
srmse = 0;
autovalor = zeros(2,R);
save_pest = zeros(VarSize(1),2,R);
save_Theta = zeros(VarSize(1),VarSize(2),R);
st = 0;
for realiza = 1:R
    tic
    yo = yt + 0*max(yt,[],2)/(1.1*FR).*randn(size(yt));
    l = L;
    FD = 100;
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
    end


    %%
    save_Theta(:,:,realiza) = Thetag;

    vtd = Histerese(length(utd),N_pesos,Thetag,utd,Ts,2);
    ytdr = ytd + 0*max(ytd,[],2)/(1.3*FR).*randn(size(ytd));
    vte = Histerese(length(u),N_pesos,Thetag,u,Ts,2);
    figure
    plot(y','k-')
    hold on
    plot(vte','r--')
    [A, B, C, D, ~] = moesp(vtd(:,n_transiente:end),ytdr(:,n_transiente:end),2,3);
    G0 = C*((eye(2)-A)\B) + D
sys = ss(A, B, C, D, 1);
pole(sys)
tzero(sys)
%     vo = Processo_v(length(utd),utd,Ts,2);
%     [A, B, C, D, ~] = moesp(vo(:,n_transiente:end)*max(vtd(1,n_transiente:end))/max(vo(1,n_transiente:end)),ytdr(:,n_transiente:end),2,4);

    autovalor(:,realiza) = eig(A);
    %% Validação
    y_est = dados_nuvem(length(uvd),N_pesos,Thetag,uvd,A,B,C,D,Ts);
    
    
    ini = 126; %Selecione esse ponto com base no máximo ou mínimo de uma das senoides da entrada
    
    aux = y_est(1,:);
    y1s = max([aux ; y1s]);
    y1i = min([aux ; y1i]);
    aux = y_est(2,:);
    y2s = max([aux ; y2s]);
    y2i = min([aux ; y2i]);
    
    rmseval = sqrt(sum((y_est - yvd).^2,2)./length(yvd))./max(yvd,[],2);
    srmse = srmse + rmseval;
    tp = toc;
    st = tp + st;
    tr = (R-realiza)*st/realiza;
    tr = seconds(tr);
    tr.Format = 'hh:mm:ss';
    fprintf('faltam aprox %s\n',tr)
end
disp(srmse/R)
figure
zplane([],[]);
hold on
th = 0:pi/50:2*pi;
plot(cos(th), sin(th), 'k:')
axis equal
b = 1;
xlim([-1 1])
ylim([-1 1])
for a = 1:size(autovalor,1)
    hold on
    plot(real(autovalor(a,:)),imag(autovalor(a,:)),'x','color',[a-1 0 b],'markersize',20);
    b = b-1;
end
h = legend({'','','','','Polo 1', 'Polo 2'});
axis equal
xlim([-1 1])
ylim([-1 1])
xlabel('parte real')
ylabel('parte imaginária')
set(gca,'FontName','Times','FontSize',14)
set(gcf,'Color','w','Position',[50 200 500 350])

theta = [0.1, 0.9, 1, 0.1, 0.9, 0.3, 0.2, 0.6, 0.8, 0.4, 0.9, 0.9, 0.35, 1, 0.25, 1,1;
    0.1, 1, 0.9, 0.1, 0.7, 0.1, 0.8, 0.4, 0.1, 0.5, 0.6,  1,  0.4,  1, 0.4, 1,1];

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

%%

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
% cruz1 = 0;
% cruz2 = 0;
% pcruz11 = 0; pcruz12 = 0; pcruz21 = 0; pcruz22 = 0;
% for k = 1:ceil(qd/2)
%     if y1d(k) < y1d(qd+2-k) && cruz1 == 0
%         pcruz11 = k;
%         cruz1 = 1;
%     elseif y1d(k) > y1d(qd+2-k) && cruz1 == 1
%         pcruz12 = k;
%         cruz1 = 2;
%     end
%     if y2d(k) < y2d(qd+2-k) && cruz2 == 0
%         pcruz21 = k;
%         cruz2 = 1;
%     elseif y2d(k) > y2d(qd+2-k) && cruz2 == 1
%         pcruz22 = k;
%         cruz2 = 2;
%     end
% end
% % fecha o circulo de dentro nos pontos de cruzamento
% 
% if pcruz11 ~= 0
%     u1d(1:pcruz11) = u1d(pcruz11);
%     u1d(end-pcruz11+1:end) = u1d(pcruz11);
%     y1d(1:pcruz11) = y1d(pcruz11);
%     y1d(end-pcruz11+1:end) = y1d(pcruz11);
% end
% if pcruz12 ~= 0
%     u1d(pcruz12:pcruz12+2*((ceil(qd/2)-pcruz12))+1) = u1d(pcruz12-1);
%     y1d(pcruz12:pcruz12+2*((ceil(qd/2)-pcruz12))+1) = y1d(pcruz12-1);
% end
% % u2s = u2;
% if pcruz21 ~= 0
%     u2d(1:pcruz21) = u2d(pcruz21);
%     u2d(end-pcruz21+1:end) = u2d(pcruz21);
%     y2d(1:pcruz21) = y2d(pcruz21);
%     y2d(end-pcruz21+1:end) = y2d(pcruz21);
% end
% if pcruz22 ~= 0
%     u2d(pcruz22:pcruz22+2*((ceil(qd/2)-pcruz22))+1) = u2d(pcruz22-1);
%     y2d(pcruz22:pcruz22+2*((ceil(qd/2)-pcruz22))+1) = y2d(pcruz22-1);
% end

%%

figure
patch([u1d fliplr(u1f)],[y1d fliplr(y1f)],[0, 0, 1],...
    'facealpha',0.5,'edgealpha',0,'EdgeColor',[0, 0, 1]);
hold on
patch([u2d fliplr(u2f)],[y2d fliplr(y2f)],[0 0.6 0],...
    'facealpha',0.5,'edgealpha',0,'EdgeColor',[0 0.6 0]);

plot(uvd(1,n_transientv:n_transientv+qd),yvd(1,n_transientv:n_transientv+qd),'k-')
plot(uvd(2,n_transientv:n_transientv+qd),yvd(2,n_transientv:n_transientv+qd),'r-')
yticks(-60:20:60)
xticks(-5:5)
box on
xlabel('$u(k)$','interpreter','latex')
ylabel('$y(k)$','interpreter','latex')
legend({'$\hat{y}_1(k)$','$\hat{y}_2(k)$','$y_1(k)$','$y_2(k)$'},'interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 290 500 350])
ylim([-60 60])

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
ylabel('$y(k)$','interpreter','latex')
legend({'$\hat{y}_1(k)$','$\hat{y}_2(k)$','$y_1(k)$','$y_2(k)$'},'interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 290 500 350])

%save exemplo2_1000_sim
%system('shutdown -s')
%% Funções


function v_est = Histerese(N_dados,N_pesos,pesos,ute,Ts,tipo)

theta = pesos(:,1:N_pesos+1);
theta_hat = pesos(:,N_pesos+2:end-2);
xi = pesos(:,end-1);
beta = pesos(:,end);
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
        rj = xi.*vecnorm(ute,Inf,2)*(j-1)/(N_pesos)+ beta.*(abs(delta)./Ts)*5./(max(ute,[],2)/Ts);
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
theta_hat = pesos(:,N_pesos+2:end-2);
xi = pesos(:,end-1);
beta = pesos(:,end);
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
        rj = xi.*vecnorm(ute,Inf,2)*(j-1)/(N_pesos)+ beta.*(abs(delta)./Ts)*5./(max(ute,[],2)/Ts);
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
% A = [1.3150   -0.0757;
%     3.4918    0.2870];
% 
% B = [0.9873    0.9873;
%     0.7819    0.7819];
% 
% C = 0.01*[1.1006    1.1006;
%        .7*1.0340  .7*1.0340];
   
   
 A = [0.4    1.7;
    -0.15    1.3];

B = [0.1    0.1;
    1    1];

C = [0.1    0.1;
    0.2    0.2]/6;
    
% A = [1.2018 4.6497; -0.0480 0.2370];
% B = [1 1; 1 1];
% C = 0.01*[1 1; 2 2];
D = [0 0; 0 0];

theta = [0.1, 0.9, 1, 0.1, 0.9, 0.3, 0.2, 0.6, 0.8, 0.4, 0.9;
         0.1, 1, 0.9, 0.1, 0.7, 0.1, 0.8, 0.4, 0.1, 0.5, 0.6];
theta_hat = [0.9, 0.35, 1, 0.25;
               1,  0.4, 1, 0.4];

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
        rj = xi.*vecnorm(u,Inf,2)*(j-1)/(N_pesos) + beta.*(abs(delta)./Ts)./(max(u,[],2)/Ts);
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
% A = [1.3150   -0.0757;
%     3.4918    0.2870];
% 
% B = [0.9873    0.9873;
%     0.7819    0.7819];
% 
% C = 0.01*[1.1006    1.1006;
%     .7*1.0340    .7*1.0340];

A = [0.4    1.7;
    -0.15    1.3];

B = [0.1    0.1;
    1    1];

C = [0.1    0.1;
    0.2    0.2]/6;

% A = [1.2018 4.6497; -0.0480 0.2370];
% B = [1 1; 1 1];
% C = 0.01*[1 1; 2 2];
D = [0 0; 0 0];

theta = [0.1, 0.9, 1, 0.1, 0.9, 0.3, 0.2, 0.6, 0.8, 0.4, 0.9;
         0.1, 1, 0.9, 0.1, 0.7, 0.1, 0.8, 0.4, 0.1, 0.5, 0.6];
theta_hat = [0.9, 0.35, 1, 0.25;
               1,  0.4, 1, 0.4];

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
        rj = xi.*vecnorm(u,Inf,2)*(j-1)/(N_pesos) + beta.*(abs(delta)./Ts)./(max(u,[],2)/Ts);
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
theta_hat = pesos(:,N_pesos+2:end-2);
xi = pesos(:,end-1);
beta = pesos(:,end);
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
        rj = xi.*vecnorm(uvd,Inf,2)*(j-1)/(N_pesos) + beta.*(abs(delta)./Ts)./(max(uvd,[],2)/Ts);
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


