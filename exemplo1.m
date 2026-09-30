% Code used in example 1 of the paper
% It takes a long time to run due to R = 1000. Changing its value
% uses fewer Monte Carlo simulations, but makes execution faster.

clc
clear
close all
rng(1) 
nu = 2;              % number of inputs
FR = 24;             % SNR factor for noise
R = 1000;            % Number of swarm realizations
N_weights = 10;      % Number of weights
Nv = N_weights+5;    % Total number of parameters
VarSize = [2 Nv];    % Size of Decision Variables Matrix
VarMin = 0;          % Lower Bound of Variables
VarMax = 2;          % Upper Bound of Variables

% PSO Parameters
alpha = 0.7;         % Mutation probability
Kmax = 80;           % Maximum Number of Iterations
P = 30;              % Population Size (Swarm Size)
b = ceil(0.05*P);
Ts = 0.001;

%% Identify dynamics
num_bloco_hankel = 2;
Amp = 1.1;
f1td = 0.001*(10+ 5*rand);
f2td = 0.5*(20+ 5*rand);
f3td = 0.5*(30+ 5*rand);
f4td = 0.5*(40+ 5*rand);
f5td = 0.5*(50+ 5*rand);

f1td2 = 0.001*(15+ 5*rand);
f2td2 = 0.5*(25+ 5*rand);
f3td2 = 0.5*(35+ 5*rand);
f4td2 = 0.5*(45+ 5*rand);
f5td2 = 0.5*(55+ 5*rand);

Nd = ceil(4/min(f1td,f2td)/Ts);
t = 0:Ts:Nd*Ts - Ts;
n_transient = ceil(1/min(f1td,f2td)/Ts);

utd = [Amp*sin(2*pi*f1td.*t)+Amp*sin(2*pi*f2td.*t)+Amp*sin(2*pi*f3td.*t)+Amp*sin(2*pi*f4td.*t)+Amp*sin(2*pi*f5td.*t);
       Amp*sin(2*pi*f1td2.*t)+Amp*sin(2*pi*f2td2.*t)+Amp*sin(2*pi*f3td2.*t)+Amp*sin(2*pi*f4td2.*t)+Amp*sin(2*pi*f5td2.*t)];
Maxutd = max(abs(utd'), [], 1);
ytd = Process(length(utd),utd,Maxutd,Ts,2);

thetao = [0.1, 0.3, 0.7, 0.5, 0.1, 0.3, 0.2, 0.6, 0.4, 0.4, 0.1;
          0.1, 1, 0.8, 0.1, 0.7, 0.1, 0.8, 0.4, 0.1, 0.5, 0.6];
theta_hato = [0.9, 0.6, 1, 0.3;
               1,  0.4, 1, 0.5];
Thetao = [thetao, theta_hato];
vtd_o = Process_v(length(utd),utd,Maxutd,Ts,2,Thetao);

%% Dynamics validation
f1d = 10;
N_val = ceil(3/f1d/Ts);
Amp = 5;
t = 0:Ts:N_val*Ts - Ts;

uvd = [Amp*sin(2*pi*f1d.*t + pi); Amp*sin(2*pi*f1d.*t + pi)];
Maxuvd = max(abs(uvd'), [], 1);
yvd = Process(N_val,uvd,Maxuvd,Ts,2);
vvd = Process_v(length(uvd),uvd,Maxuvd,Ts,2,Thetao);
n_transientv = 1/f1d/Ts;

%% Hysteresis identification
L = 130;
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
Maxut = max(abs(ut'), [], 1);
yt = Process(length(ut),ut,Maxut,Ts,1);
uta = ut(:,L:L:end);
n_transienth = 25;
uta1 = uta(1,1:2*Q);
uta2 = uta(2,length(uta)/2+1:length(uta)/2+2*Q);
u = [uta1; uta2];
Maxu = max(abs(u'), [], 1);
t = 0:Ts:length(u)*Ts - Ts;
Ne = length(t);

%% Steady-state step validation
% Build signal: sinusoidal for 2 cycles + held step
f_test = f1d;
N_sin  = ceil(1.6/f_test/Ts);
N_hold = ceil(1/f_test/Ts);    % hold time (1 cycle)
u_hold_val = 2.67913;          % held value (choose within operational range)

t_test = 0:Ts:(N_sin + N_hold)*Ts - Ts;
u_test_1 = [Amp*sin(2*pi*f_test*(0:Ts:N_sin*Ts-Ts) + pi), ...
             u_hold_val*ones(1, N_hold)];
u_test = [u_test_1; u_test_1];
Maxutest = max(abs(u_test'), [], 1);

% Real process
y_proc_test = Process(length(u_test), u_test,Maxutest,Ts, 2);
y_hold_all = zeros(R, 2, length(u_test));  % outputs during hold

%% PSO Parameters
phi1 = 2.05;
phi2 = 2.05;
phi = phi1+phi2;
chi = 2/abs(2-phi-sqrt(phi^2-4*phi));
c1 = phi1;       % Personal Learning Coefficient
c2 = phi2;       % Global Learning Coefficient
VelMax = 0.05*VarMax;

%% Initialize swarm parameters
y1s = -inf(1,length(yvd));
y2s = y1s;
y1i = -y1s;
y2i = -y1s;
srmse = 0;
eigenval = zeros(3,R);
eigenval_o = zeros(3,R);
save_zero = zeros(3,R);
save_rmse = zeros(2,R);
erro = zeros(R,2,101);
save_Theta = zeros(VarSize(1),VarSize(2),R);
st = 0;

for realiza = 1:R
    tic
    yo = yt + max(yt,[],2)/(1.1*FR).*randn(size(yt));
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
    % Iteration parameters
    Theta = zeros([VarSize P]);
    V = zeros([VarSize P]);
    F = inf(VarSize(1),P);
    
    % Global parameters
    Thetag = zeros(VarSize);
    Vg = zeros(VarSize);
    Fg = inf(VarSize(1),1);
    
    % Personal parameters
    Thetap = zeros([VarSize P]);
    Vp = zeros([VarSize P]);
    Fp = inf(VarSize(1),P);
    
    for i = 1:P
        Theta(:,:,i) = unifrnd(0,VarMax,VarSize);
        Theta(:,N_weights+1:end,i) = min(Theta(:,N_weights+1:end,i),1);
        V(:,:,i) = unifrnd(-VelMax,VelMax,VarSize);
        
        % Evaluation
        F(:,i) = Cost(Ne,N_weights,Theta(:,:,i),u,Maxu,y,n_transienth,Ts);
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
    
    %% PSO Main Loop according to Yang2013
    mark = 0;
    Mp = zeros([VarSize P]);
    
    for k = 1:Kmax
        mark = mark + 1;
        Mg = zeros(VarSize);
        for l = 1:VarSize(1)
            % Evaluate weighted personal and global best position
            [sorted_cost, sorted_ind] = sort(Fp(l,:));
            for q = 1:b
                mu = (1/sorted_cost(q))/sum(1./sorted_cost(1:b));
                Mg(l,:) = Mg(l,:) + mu*Thetap(l,:,sorted_ind(q));
            end
            for i = 1:P
                if i ~= 1
                    Mp(l,:,sorted_ind(i)) = (Fp(l,sorted_ind(i))*Thetap(l,:,sorted_ind(i-1))...
                        + Fp(l,sorted_ind(i-1))*Thetap(l,:,sorted_ind(i)))/ ...
                        (Fp(l,sorted_ind(i))+Fp(l,sorted_ind(i-1)));
                else
                    Mp(l,:,sorted_ind(i)) = Thetap(l,:,sorted_ind(i));
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
            F(:,i) = Cost(Ne,N_weights,Theta(:,:,i),u,Maxu,y,n_transienth,Ts);
        end
        
        for i = 1:P
            for l = 1:VarSize(1)
                % Update Personal Best
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
                    % disp(['Iteration ' num2str(k) ': Best Cost ' num2str(l)  ' = ' num2str(Fg(l))]);
                end
            end
        end
        
        rho = rand(VarSize(1),1);
        if mark > 10 % Mutation takes effect if after 10 iterations there is no improvement
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
    
    vtd = Hysteresis(length(utd),N_weights,Thetag,utd,Maxutd,Ts,2);
    noise_p2 = max(ytd,[],2)/(1.2*FR) .* randn(size(ytd));
    ytdr = ytd + noise_p2;
    snr_p2 = 10 * log10( sum(ytd.^2, 2) ./ sum(noise_p2.^2, 2) );
    vte = Hysteresis(length(u),N_weights,Thetag,u,Maxu,Ts,2);
    
    [A2, B2, C2, D2, ~] = moesp(vtd_o(:,n_transient:end).*(max(vtd(:,n_transient:end),...
        [],2)./max(vtd_o(:,n_transient:end),[],2)),ytdr(:,n_transient:end),3,3);
    [A, B, C, D, ~] = moesp_po(vtd(:,n_transient:end),ytdr(:,n_transient:end),3,3);
    sys = ss(A, B, C, D, Ts);
    
    eigenval_o(:,realiza) = eig(A2);
    eigenval(:,realiza) = eig(A);
    save_zero(:,realiza) = tzero(sys);
    
    %% Validation
    y_r = swarm_data(length(u_test),N_weights,Thetag,u_test,Maxutest,A,B,C,D,Ts);
    y_hold_all(realiza,:,:) = y_r;
    
    y_est = swarm_data(length(uvd),N_weights,Thetag,uvd,Maxuvd,A,B,C,D,Ts);
    
    ini = 126; % Select this point based on the maximum or minimum of one of the input sine waves
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
    fprintf('approx %s remaining\n',tr)
end
save example1_1000_sim

%%
% load example1_1000_sim
yo_steadystate = yo(:,1:80);
yo(:,1:80) = yt(:,1:80) + max(yt,[],2)/(1.1*FR).*randn(1,80);

figure
plot(0:Ts:length(ut)*Ts-Ts,ut(1,:),'r-')
hold on
plot(0:Ts:length(ut)*Ts-Ts,ut(2,:),'k-')
plot(0:Ts:length(ut)*Ts-Ts,yo(1,:),'b-')
plot(0:Ts:length(ut)*Ts-Ts,yo(2,:),'g-','color',[0 0.7 0 0.7])
plot((L-1)*Ts:L*Ts:length(ut)*Ts-Ts,yo_steadystate(1,:),'dr','MarkerSize',3,'MarkerFaceColor', 'r')
plot((L-1)*Ts:L*Ts:length(ut)*Ts-Ts,yo_steadystate(2,:),'dk','MarkerSize',3,'MarkerFaceColor', 'k')
legend('$u_1$','$u_2$','$y_1$','$y_2$','$\bar{y}_1$','$\bar{y}_2$')
legend({},'interpreter','latex','orientation','horizontal')
xlabel('$kT_s$(s)','interpreter','latex')
ylabel('amplitude','interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 200 700 350])

figure
plot(0:Ts:length(ut)*Ts-Ts,1*ut(1,:),'k-')
hold on
plot(0:Ts:length(ut)*Ts-Ts,1*ut(2,:),'r-')
plot(0:Ts:length(ut)*Ts-Ts,yo(1,:),'b-')
plot(0:Ts:length(ut)*Ts-Ts,yo(2,:),'g-','color',[0 0.7 0 0.7])
xlim([2.6 5.2])
ylim([-20 20])
plot((L-1)*Ts:L*Ts:length(ut)*Ts-Ts,yo_steadystate(1,:),'*k','MarkerSize',6,'MarkerFaceColor', 'k')
plot((L-1)*Ts:L*Ts:length(ut)*Ts-Ts,yo_steadystate(2,:),'*r','MarkerSize',6,'MarkerFaceColor', 'r')
legend('$u_1$','$u_2$','$y_1$','$y_2$','$\bar{y}_1$','$\bar{y}_2$')
legend({},'interpreter','latex','orientation','horizontal','NumColumns',3)
xlabel('$kT_s$(s)','interpreter','latex')
ylabel('amplitude','interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 200 420 350])
grid off

figure
plot(0:Ts:length(ut)*Ts-Ts,1*ut(1,:),'k-')
hold on
plot(0:Ts:length(ut)*Ts-Ts,1*ut(2,:),'r-')
plot(0:Ts:length(ut)*Ts-Ts,yo(1,:),'b-')
plot(0:Ts:length(ut)*Ts-Ts,yo(2,:),'g-','color',[0 0.7 0 0.7])
xlim([7.8 10.4])
ylim([-20 20])
plot((L-1)*Ts:L*Ts:length(ut)*Ts-Ts,yo_steadystate(1,:),'*k','MarkerSize',6,'MarkerFaceColor', 'k')
plot((L-1)*Ts:L*Ts:length(ut)*Ts-Ts,yo_steadystate(2,:),'*r','MarkerSize',6,'MarkerFaceColor', 'r')
legend('$u_1$','$u_2$','$y_1$','$y_2$','$\bar{y}_1$','$\bar{y}_2$')
legend({},'interpreter','latex','orientation','horizontal','NumColumns',3)
xlabel('$kT_s$(s)','interpreter','latex')
ylabel('amplitude','interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 200 420 350])
grid off

Qf = 301;
figure
plot(0:Ts:(Qf-1)*Ts,vtd(1,1:Qf),'k--')
hold on
plot(0:Ts:(Qf-1)*Ts,vtd(2,1:Qf),'r--')
plot(0:Ts:(Qf-1)*Ts,ytdr(1,1:Qf),'b-')
plot(0:Ts:(Qf-1)*Ts,ytdr(2,1:Qf),'-','color',[0 0.7 0])
legend('$\hat{v}_1$','$\hat{v}_2$','$y_1$','$y_2$')
legend({},'interpreter','latex','orientation','horizontal')
xlabel('$kT_s$(s)','interpreter','latex')
ylabel('amplitude')
xlim([0 0.3])
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 200 500 350])

%%
figure
subplot(1,2,1)
histogram(save_rmse(1,:), 10, 'FaceColor', [0 0 1], 'FaceAlpha', 0.5)
hold on
xline(mean(save_rmse(1,:)), 'b--', 'LineWidth', 1.5)
xlabel('NRMSE', 'FontName', 'Times', 'FontSize', 16)
xlim([0 0.07])
ylim([0 400])
ylabel('count', 'FontName', 'Times', 'FontSize', 16)
box on
legend('NRMSE$(y_1)$', '$\overline{\mathrm{NRMSE}}(y_1)$', 'Interpreter', 'latex')
set(gca,'FontName','Times','FontSize',16)

subplot(1,2,2)
hold on
histogram(save_rmse(2,:), 13, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.5)
xline(mean(save_rmse(2,:)), '--','color',[.2 .39 .2 1], 'LineWidth', 1.5)
xlim([0 0.07])
ylim([0 400])
xlabel('NRMSE', 'FontName', 'Times', 'FontSize', 16)
ylabel('count', 'FontName', 'Times', 'FontSize', 16)
box on
legend('NRMSE$(y_2)$', '$\overline{\mathrm{NRMSE}}(y_2)$', 'Interpreter', 'latex')
set(gca,'FontName','Times','FontSize',16)
set(gcf, 'Color', 'w', 'Position', [50 200 800 200])

disp(srmse/R)

%%
figure
subplot(1,2,1)
hold on
for realiza = 1:R
    % squeeze transforms y(realiza, 1, :) from 1x1xN to an Nx1 vector
    y1 = squeeze(y_hold_all(realiza, 1, :));
    y2 = squeeze(y_hold_all(realiza, 2, :));
    if realiza == 1
        h_y1 = plot((1:length(y1))*Ts,y1, 'Color', [.55 .64 .82 1]); 
        h_y2 = plot((1:length(y2))*Ts,y2, 'Color', [.47 .75 .45 1]); 
    else
        plot((1:length(y1))*Ts,y1, 'Color', [.55 .64 .82 1]); 
        plot((1:length(y2))*Ts,y2, 'Color', [.47 .75 .45 1]); 
    end
end

h_proc_y1 = plot((1:length(y_proc_test))*Ts,y_proc_test(1,:), 'k-');
h_proc_y2 = plot((1:length(y_proc_test))*Ts,y_proc_test(2,:), 'r-');

mean_y1 = squeeze(mean(y_hold_all(:, 1, :), 1)); 
mean_y2 = squeeze(mean(y_hold_all(:, 2, :), 1));

h_med1 = plot((1:length(mean_y1))*Ts,mean_y1, 'b--','linewidth',1.5); 
h_med2 = plot((1:length(mean_y2))*Ts,mean_y2, '--','color',[.2 .39 .2 1],'linewidth',1.5);

xlabel('$kT_s$(s)', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 16)
ylabel('$y(kT_s)$', 'Interpreter', 'latex', 'FontSize', 16)
legend([h_y1, h_y2, h_proc_y1, h_proc_y2, h_med1, h_med2], ...
       {'$\hat{y}_1(kT_s)$', '$\hat{y}_2(kT_s)$', '$y_1(kT_s)$',...
       '$y_2(kT_s)$', 'mean $\hat{y}_1$', 'mean $\hat{y}_2$'}, ...
       'Interpreter', 'latex');
xlim([0 260]*Ts)
box on
set(gca,'FontName','Times','FontSize',16)

subplot(1,2,2)
hold on
for realiza = 1:R
    y1 = squeeze(y_hold_all(realiza, 1, :));
    y2 = squeeze(y_hold_all(realiza, 2, :));
    if realiza == 1
        h_y1 = plot((1:length(y1))*Ts,y1, 'Color', [.55 .64 .82 1]);
        h_y2 = plot((1:length(y2))*Ts,y2, 'Color', [.47 .75 .45 1]); 
    else
        plot((1:length(y1))*Ts,y1, 'Color', [.55 .64 .82 1]); 
        plot((1:length(y2))*Ts,y2, 'Color', [.47 .75 .45 1]); 
    end
end

h_proc_y1 = plot((1:length(y_proc_test))*Ts,y_proc_test(1,:), 'k-');
h_proc_y2 = plot((1:length(y_proc_test))*Ts,y_proc_test(2,:), 'r-');

mean_y1 = squeeze(mean(y_hold_all(:, 1, :), 1)); 
mean_y2 = squeeze(mean(y_hold_all(:, 2, :), 1));

h_med1 = plot((1:length(mean_y1))*Ts,mean_y1, 'b--','linewidth',1.5); 
h_med2 = plot((1:length(mean_y2))*Ts,mean_y2, '--','color',[.2 .39 .2 1],'linewidth',1.5);

xlabel('$kT_s$(s)', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 16)
ylabel('$y(kT_s)$', 'Interpreter', 'latex', 'FontSize', 16)
xlim([160 260]*Ts)
box on
set(gca,'FontName','Times','FontSize',16)
set(gcf, 'Color', 'w', 'Position', [50 200 800 350])

%%
% =========================================================================
% SCRIPT: Monte Carlo Visualization with Double Zoom (Insets)
% Goal: Analyze Bias and Variance of Estimated Poles
% =========================================================================

% 1. Definition of true process poles
true_poles = [0.7004 + 0.3918i, 0.7004 - 0.3918i, 0.5973];

th = 0:pi/100:2*pi;
circ_x = cos(th);
circ_y = sin(th);

cores = [.47 .75 .45
         .55 .64 .82
         .55 .64 .82];

figure

% =========================================================================
% MAIN AXIS: Global view of Z-plane
% =========================================================================
ax_main = axes('Position', [0.1, 0.14, 0.52, 0.74]);
hold(ax_main, 'on');

plot(ax_main, circ_x, circ_y, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Sort the estimated eigenvalues for consistent coloring
[~, idx_sort] = sort(abs(eigenval), 1);
sorted_eigenval = zeros(size(eigenval));
for col = 1:size(eigenval, 2)
    sorted_eigenval(:, col) = eigenval(idx_sort(:, col), col);
end

% Plot estimated clouds with transparency (Alpha) in the main plot
for p = 1:3
    p_est_plot(p) = scatter(ax_main, real(sorted_eigenval(p,:)), imag(sorted_eigenval(p,:)), ...
    100, 'o', 'MarkerEdgeColor', cores(p,:), 'MarkerEdgeAlpha', 0.3, ...
    'MarkerFaceColor', cores(p,:), 'MarkerFaceAlpha', 0.1, 'HandleVisibility', 'off','DisplayName', 'est. complex eig.');
end
p_est_plot = p_est_plot(1:2); % exclude one of the complex conjugate poles
p_est_plot(1).DisplayName = 'est. real eig.'; % revert to real name for the real pole
p_est_plot(1).MarkerFaceAlpha = 0.7;
p_est_plot(1).MarkerEdgeAlpha = 1;
p_est_plot(2).MarkerFaceAlpha = 0.7;
p_est_plot(2).MarkerEdgeAlpha = 1;

% Plot true process poles (High contrast)
p_real_plot = plot(ax_main, real(true_poles(1:2)), imag(true_poles(1:2)), 'k+', ...
'MarkerSize', 30, 'MarkerFaceColor', 'k', ...
'DisplayName', 'true complex eig.');

p_real_plot2 = plot(ax_main, real(true_poles(3)), imag(true_poles(3)), 'r+', ...
'MarkerSize', 30, 'MarkerFaceColor', 'r', ...
'DisplayName', 'true real eig.');

axis(ax_main, 'equal');
xlim(ax_main, [-1.1 1.1]);
ylim(ax_main, [-1.1 1.1]);
xticklabels(ax_main, [-1, -0.5, 0 , 0.5, 1]);

xl=xlabel(ax_main, 'real part', 'FontName', 'Times', 'FontSize', 16);
xl.Position = [0, -1.27, 0];
yl=ylabel(ax_main, 'imaginary part', 'FontName', 'Times', 'FontSize', 16);
yl.Position = [-1.4, 0.1, 0];

set(ax_main, 'FontName', 'Times', 'FontSize', 16);
legend(ax_main, [p_est_plot,p_real_plot,p_real_plot2], 'Location', 'southwest');
box(ax_main, 'on');

% =========================================================================
% INSET 1: Zoom on Real Pole (p3 = 0.5973)
% =========================================================================
ax_inset_real = axes('Position', [0.71, 0.14, 0.25, 0.3]);
hold(ax_inset_real, 'on');

scatter(ax_inset_real, real(sorted_eigenval(1,:)), imag(sorted_eigenval(1,:)), ...
80, 'o', 'MarkerEdgeColor', cores(1,:), 'MarkerEdgeAlpha', 0.3, ...
'MarkerFaceColor', cores(1,:), 'MarkerFaceAlpha', 0.1);

plot(ax_inset_real, real(true_poles(3)), imag(true_poles(3)), 'r+', ...
'MarkerSize', 40);

xlim(ax_inset_real, [0.57, 0.71]);
ylim(ax_inset_real, [-0.01, 0.01]);

set(ax_inset_real, 'FontName', 'Times', 'FontSize', 16);
box(ax_inset_real, 'on');

% =========================================================================
% INSET 2: Zoom on Upper Complex Pole (p1 = 0.7004 + 0.3918j)
% =========================================================================
ax_inset_complex = axes('Position', [0.71, 0.58, 0.25, 0.3]);
hold(ax_inset_complex, 'on');

% Identify which row (2 or 3) has positive imaginary part
im_pos_idx = 2;
if mean(imag(sorted_eigenval(3,:))) > 0
    im_pos_idx = 3;
end

scatter(ax_inset_complex, real(sorted_eigenval(im_pos_idx,:)), imag(sorted_eigenval(im_pos_idx,:)), ...
80, 'o', 'MarkerEdgeColor', cores(im_pos_idx,:), 'MarkerEdgeAlpha', 0.3, ...
'MarkerFaceColor', cores(im_pos_idx,:), 'MarkerFaceAlpha', 0.1);

plot(ax_inset_complex, 0.7004, 0.3918, 'k+', ...
'MarkerSize', 40, 'MarkerFaceColor', 'k');

xlim(ax_inset_complex, [0.692, 0.72]);
ylim(ax_inset_complex, [0.35, 0.405]);
yticklabels(ax_inset_complex, [0.36, 0.38, 0.4]);

set(ax_inset_complex, 'FontName', 'Times', 'FontSize', 16);
box(ax_inset_complex, 'on');
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 200 550 350])

%%
% 1. Original poles and sorting
original_poles = [0.7004 + 0.3918i, 0.7004 - 0.3918i, 0.5973];
true_magnitudes = sort(abs(original_poles));

% 2. Estimated magnitudes sorted by column
estimated_magnitudes = sort(abs(eigenval), 1);

% 3. Calculate absolute error independently
error_pole_true = abs(estimated_magnitudes(1,:) - true_magnitudes(1));
error_pole_complex = abs(estimated_magnitudes(2,:) - true_magnitudes(2));

% 4. Define intervals (edges)
val_max = max([max(error_pole_true), max(error_pole_complex)]);
edges = [0, 0.005, 0.01, 0.02, 0.05, max(0.2, val_max)];

counts_true = histcounts(error_pole_true, edges);
counts_complex = histcounts(error_pole_complex, edges);

% Assemble counts matrix for grouped plot (Size: num_bins x 2)
counts_matrix = [counts_true', counts_complex'];

% 6. Generate X-axis labels
num_bins = length(edges) - 1;
labels_x = cell(1, num_bins);
for k = 1:num_bins
    if k == num_bins
        labels_x{k} = sprintf('[%g, %g]', edges(k), edges(k+1));
    else
        labels_x{k} = sprintf('[%g, %g)', edges(k), edges(k+1));
    end
end

figure
b = bar(1:num_bins, counts_matrix, 'grouped', 'EdgeColor', 'k', 'BarWidth', 0.85);

b(1).FaceColor = [.55 .64 .82]; 
b(2).FaceColor = [.47 .75 .45]; 

xticks(1:num_bins);
xticklabels(labels_x);

xlabel('absolute eigenvalue magnitude error', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('count', 'FontName', 'Times', 'FontSize', 16);

legend({'real eigenvalue', 'complex eigenvalue'}, ...
       'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 16);
set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 200 840 350])
box on;
grid on;
ax = gca;
ax.XGrid = 'off'; 
ax.YGrid = 'on';

%%
modulos_zeros = abs(save_zero);
max_modulos = max(modulos_zeros);
val_max = ceil(max(max_modulos)); % Round up the true maximum

% 1. Define custom intervals (edges).
edges = [0, 1, 5, 10, 20, 100, max(100, val_max)];

counts = histcounts(max_modulos, edges);

num_bins = length(counts);
labels_x = cell(1, num_bins);
for k = 1:num_bins
    if k == num_bins
        labels_x{k} = sprintf('[%g, %g]', edges(k), edges(k+1));
    else
        labels_x{k} = sprintf('[%g, %g)', edges(k), edges(k+1));
    end
end

figure
b = bar(1:num_bins, counts, 'FaceColor', [.55 .64 .82], 'EdgeColor', 'k', 'BarWidth', 0.8);
hold on;

% 5. Stability/minimum phase dividing line
xline(1.5, 'k--', 'LineWidth', 2, 'Label', '|zero| = 1', ...
    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', ...
    'FontName', 'Times', 'FontSize', 16, 'Color', 'k');
xlim([0.5 6.5])

xticks(1:num_bins);
xticklabels(labels_x);

xlabel('max(|zero|)', 'FontName', 'Times', 'FontSize', 16);
ylabel('count', 'FontName', 'Times', 'FontSize', 16);

set(gca,'FontName','Times','FontSize',16)
set(gcf,'Color','w','Position',[50 200 840 350])
box on;
grid on;
ax = gca;
ax.XGrid = 'off'; 
ax.YGrid = 'on';

%%
theta = [0.1, 0.9, 1, 0.1, 0.9, 0.3, 0.2, 0.6, 0.8, 0.4, 0.9, 0.9, 0.35, 1, 0.25;
         0.1, 1, 0.9, 0.1, 0.7, 0.1, 0.8, 0.4, 0.1, 0.5, 0.6,  1,  0.4,  1, 0.4];
for a = 1:size(save_Theta,1)
    figure
    plot(theta(a,:),'*','color',[0 0 0],'markersize',8)
    hold on
    plot(reshape(save_Theta(a,:,:),[VarSize(2),R]),'*','color',[1 0 0],'markersize',8)
    legend('process','estimated')
    xlabel('variable position in \Theta')
    ylabel('value')
    xlim([0 VarSize(2)+1])
    set(gca,'FontName','Times','FontSize',14)
    set(gcf,'Color','w','Position',[50 200 500 350])
end

%%
% Bounds notation:
% s - upper (superior)
% i - lower (inferior)
% d - inside (dentro)
% f - outside (fora)

qd = 1/f1d/Ts;
u1f = uvd(1,ini:ini+qd-1);
u1d = u1f;
u2f = uvd(2,ini:ini+qd-1);
u2d = u2f;

% Indices for the second half of the cycle
idx_start = ceil(qd/2) + 1;
idx_end = qd;

% Corresponding indices in the original vectors (y1s, y1i, etc.)
source_idx_start = ini + ceil(qd/2);
source_idx_end = ini + qd - 1;

% Loop construction for y1
y1d = y1s(ini:ini+qd-1);
y1d(idx_start:idx_end) = y1i(source_idx_start:source_idx_end); 
y1f = y1i(ini:ini+qd-1);
y1f(idx_start:idx_end) = y1s(source_idx_start:source_idx_end); 

% Loop construction for y2
y2d = y2s(ini:ini+qd-1);
y2d(idx_start:idx_end) = y2i(source_idx_start:source_idx_end); 
y2f = y2i(ini:ini+qd-1);
y2f(idx_start:idx_end) = y2s(source_idx_start:source_idx_end); 

% close all paths
u1f = [u1f u1f(1)];
y1f = [y1f y1f(1)];
u1d = [u1d u1d(1)];
y1d = [y1d y1d(1)];
u2f = [u2f u2f(1)];
y2f = [y2f y2f(1)];
u2d = [u2d u2d(1)];
y2d = [y2d y2d(1)];

% Identify crossing points
cross1 = 0;
cross2 = 0;
pcross11 = 0; pcross12 = 0; pcross21 = 0; pcross22 = 0;

for k = 1:ceil(qd/2)
    if y1d(k) < y1d(qd+2-k) && cross1 == 0
        pcross11 = k;
        cross1 = 1;
    elseif y1d(k) > y1d(qd+2-k) && cross1 == 1
        pcross12 = k;
        cross1 = 2;
    end
    if y2d(k) < y2d(qd+2-k) && cross2 == 0
        pcross21 = k;
        cross2 = 1;
    elseif y2d(k) > y2d(qd+2-k) && cross2 == 1
        pcross22 = k;
        cross2 = 2;
    end
end

% close the inner circle at crossing points
if pcross11 ~= 0
    u1d(1:pcross11) = u1d(pcross11);
    u1d(end-pcross11+1:end) = u1d(pcross11);
    y1d(1:pcross11) = y1d(pcross11);
    y1d(end-pcross11+1:end) = y1d(pcross11);
end

if pcross12 ~= 0
    u1d(pcross12:pcross12+2*((ceil(qd/2)-pcross12))+1) = u1d(pcross12-1);
    y1d(pcross12:pcross12+2*((ceil(qd/2)-pcross12))+1) = y1d(pcross12-1);
end

if pcross21 ~= 0
    u2d(1:pcross21) = u2d(pcross21);
    u2d(end-pcross21+1:end) = u2d(pcross21);
    y2d(1:pcross21) = y2d(pcross21);
    y2d(end-pcross21+1:end) = y2d(pcross21);
end

if pcross22 ~= 0
    u2d(pcross22:pcross22+2*((ceil(qd/2)-pcross22))+1) = u2d(pcross22-1);
    y2d(pcross22:pcross22+2*((ceil(qd/2)-pcross22))+1) = y2d(pcross22-1);
end

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
set(gcf,'Color','w','Position',[50 290 500 350])

%%
erro1 = squeeze(erro(:,1,:));
erro1 = erro1*100/max(yvd(1, ini:fim));
erro2 = squeeze(erro(:,2,:));
erro2 = erro2*100/max(yvd(2, ini:fim));
u_ciclo = uvd(1, ini:fim);

% Separate ascending and descending branches
idx_asc = [true; diff(u_ciclo') > 0];
idx_desc = ~idx_asc;

figure
% --- SUBPLOT 1 ---
subplot(2,1,1)
hold on
x_asc = u_ciclo(idx_asc);
y_asc1 = mean(erro1(:,idx_asc), 1);
a1 = area(x_asc, y_asc1, 'FaceColor', [0 0 1], 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'LineWidth', 1);

x_desc = u_ciclo(idx_desc);
y_desc1 = mean(erro1(:,idx_desc), 1);
[x_desc_ord, idx_sort] = sort(x_desc); 
y_desc1_ord = y_desc1(idx_sort);       

a2 = area(x_desc_ord, y_desc1_ord, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.4, 'EdgeColor', [.2 .39 .2], 'LineWidth', 1);
xticks(-5:2:5)
xlabel('$u_1$', 'Interpreter', 'latex')
ylabel('$\bar{e}_1$(\%)', 'Interpreter', 'latex')
legend([a1, a2], {'$\bar{e}_1^+$', '$\bar{e}_1^-$'}, 'interpreter','latex')
set(gca, 'FontName', 'Times', 'FontSize', 16)
box on

% --- SUBPLOT 2 ---
subplot(2,1,2)
hold on
y_asc2 = mean(erro2(:,idx_asc), 1);
a1 = area(x_asc, y_asc2, 'FaceColor', [0 0 1], 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'LineWidth', 1);

y_desc2 = mean(erro2(:,idx_desc), 1);
y_desc2_ord = y_desc2(idx_sort); 

a2 = area(x_desc_ord, y_desc2_ord, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.4, 'EdgeColor', [.2 .39 .2], 'LineWidth', 1);
xticks(-5:2:5)
xlabel('$u_2$', 'Interpreter', 'latex')
ylabel('$\bar{e}_2$(\%)', 'Interpreter', 'latex')
legend([a1, a2], {'$\bar{e}_2^+$', '$\bar{e}_2^-$'}, 'interpreter','latex')
set(gca, 'FontName', 'Times', 'FontSize', 16)
box on
set(gcf, 'Color', 'w', 'Position', [50 200 500 350])

%% Functions
function v_est = Hysteresis(N_dados,N_weights,pesos,ute,Maxute,Ts,tipo)
theta = pesos(:,1:N_weights+1);
theta_hat = pesos(:,N_weights+2:end);
xi = ones(size(pesos(:,end)));
beta = 0*pesos(:,end);
v_est = zeros(size(ute,1),N_dados);
Ph = zeros(size(ute,1),N_weights);
Bh = zeros(size(ute,1),1);
for k = 1:N_dados
    BL = Maxute.*theta_hat(:,1).*tanh(theta_hat(:,2).*ute(:,k));
    BR = Maxute.*theta_hat(:,3).*tanh(theta_hat(:,4).*ute(:,k));
    if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
        delta = zeros(size(ute,1),1);
        Ph = zeros(size(ute,1),N_weights);
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
    for j = 1:N_weights
        rj = xi.*Maxute*(j-1)/(N_weights);
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

function F = Cost(N_dados,N_weights,pesos,ute,Maxute,vte,n_transient,Ts)
theta = pesos(:,1:N_weights+1);
theta_hat = pesos(:,N_weights+2:end);
xi = ones(size(pesos(:,end)));
beta = 0*pesos(:,end);
v_est = zeros(size(ute,1),N_dados);
Ph = zeros(size(ute,1),N_weights);
Bh = zeros(size(ute,1),1);
for k = 1:N_dados
    BL = Maxute.*theta_hat(:,1).*tanh(theta_hat(:,2).*ute(:,k));
    BR = Maxute.*theta_hat(:,3).*tanh(theta_hat(:,4).*ute(:,k));
    if k == 1
        delta = zeros(size(ute,1),1);
        Ph = zeros(size(ute,1),N_weights);
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
    for j = 1:N_weights
        rj = xi.*Maxute*(j-1)/(N_weights);
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

function y = Process(N_dados,u,Maxu,Ts,tipo)
% Process parameters
A = [-0.4349  -0.9854  -1.169
      0.936    1.335   0.7665
      0.3947   0.4053  1.098];
B = [-0.147  -0.1228
     -0.1222  -0.6355
     -0.3618   0.1866];
C = [0.2091   0.009051  -0.2766
     0.1895    -0.2875  0.0009427];
    
D = [0 0; 0 0];
theta = [0.1, 0.3, 0.7, 0.5, 0.1, 0.3, 0.2, 0.6, 0.4, 0.4, 0.1;
         0.1, 1, 0.8, 0.1, 0.7, 0.1, 0.8, 0.4, 0.1, 0.5, 0.6];
theta_hat = [0.9, 0.6, 1, 0.3;
              1,  0.4, 1, 0.5];
           
xi = [1; 1];
beta = [1; 1];
N_weights = size(theta,2)-1;
y = zeros(size(C,1),N_dados);
x = zeros(size(A,1),N_dados);
Ph = zeros(size(B,2),N_weights);
Bh = zeros(size(B,2),1);
for k = 1:N_dados
    BL = Maxu.*theta_hat(:,1).*tanh(theta_hat(:,2).*u(:,k));
    BR = Maxu.*theta_hat(:,3).*tanh(theta_hat(:,4).*u(:,k));
    if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
        delta = zeros(size(u,1),1);
        Ph = zeros(size(u,1),N_weights);
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
    for j = 1:N_weights
        rj = xi.*Maxu*(j-1)/(N_weights);
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

function vo = Process_v(N_dados,u,Maxu,Ts,tipo, Theta)
% Process parameters
A = [-0.4349  -0.9854  -1.169
      0.936    1.335   0.7665
      0.3947   0.4053  1.098];
B = [-0.147  -0.1228
     -0.1222  -0.6355
     -0.3618   0.1866];
C = [0.2091   0.009051  -0.2766
     0.1895    -0.2875  0.0009427];
    
D = [0 0; 0 0];
theta = Theta(:,1:11);
theta_hat = Theta(:,12:15);
xi = [1; 1];
beta = [1; 1];
N_weights = size(theta,2)-1;
vo = zeros(size(C,1),N_dados);
x = zeros(size(A,1),N_dados);
Ph = zeros(size(B,2),N_weights);
Bh = zeros(size(B,2),1);
for k = 1:N_dados
    BL = Maxu.*theta_hat(:,1).*tanh(theta_hat(:,2).*u(:,k));
    BR = Maxu.*theta_hat(:,3).*tanh(theta_hat(:,4).*u(:,k));
    if k == 1 || (k == ceil(N_dados/2)+1 && tipo == 1)
        delta = zeros(size(u,1),1);
        Ph = zeros(size(u,1),N_weights);
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
    for j = 1:N_weights
        rj = xi.*Maxu*(j-1)/(N_weights);
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

function y_est = swarm_data(N_dados,N_weights,pesos,uvd,Maxuvd,A,B,C,D,Ts)
theta = pesos(:,1:N_weights+1);
theta_hat = pesos(:,N_weights+2:end);
xi = ones(size(pesos(:,end)));
beta = 0*pesos(:,end);
y_est = zeros(size(C,1),N_dados);
x = zeros(size(A,1),N_dados);
Ph = zeros(size(B,2),N_weights);
Bh = zeros(size(B,2),1);
for k = 1:N_dados
    BL = Maxuvd.*theta_hat(:,1).*tanh(theta_hat(:,2).*uvd(:,k));
    BR = Maxuvd.*theta_hat(:,3).*tanh(theta_hat(:,4).*uvd(:,k));
    if k == 1
        delta = zeros(size(uvd,1),1);
        Ph = zeros(size(B,2),N_weights);
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
    for j = 1:N_weights
        rj = xi.*Maxuvd*(j-1)/(N_weights);
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
% Syntax
% [A,B,C,D,SS] = moesp(u,y,n,k);
%
% Description
% The moesp algorithm estimates the model matrices in state-space representation without
% treating process and/or measurement noise. It can also be used to determine the system order.
%
% Input Data
% u -> measured input data matrix of the system (row);
% y -> measured output data matrix of the system (row);
% n -> chosen order for the system model;
% k -> number of rows of the Hankel block matrix, defined as
% k = 2(maxord/nusaida). Where maxord is the maximum order the user assumes the
% system has and nusaida is the number of outputs the system has. This relationship is defined
% empirically by (Van Overschee and De Moor, 1996).
%
% Output Data
% A -> estimated dynamic matrix for the system in state space;
% B -> estimated input matrix for the system in state space;
% C -> estimated output matrix for the system in state space;
% D -> estimated direct transmission matrix for the system in state space;
% SS -> singular values matrix of the oblique projection, used to determine the system order;
%
% Implemented by Rodrigo Augusto Ricco: rodrigo.ricco@yahoo.com.br
% Universidade Federal de Minas Gerais - UFMG
% Last modification 11/14/2012

i = k;

[l, ~] = size(y);
[m, Ndados] = size(u);
j = Ndados - 2*i;
kk = 0;

for k = 1:m:2*i*m-m+1
    kk = kk+1;
    U(k:k+m-1,:) = u(:,kk:kk+j-1); 
end
kk = 0;
for k = 1:l:2*i*l-l+1
    kk = kk+1;
    Y(k:k+l-1,:) = y(:,kk:kk+j-1); 
end
Uf = U(i*m+1:2*i*m,:); 
Yf = Y(i*l+1:2*i*l,:); 

H = [Uf; Yf]; 

L = triu(qr(H'))'; 
im=i*m;
il=i*l;
L11 = L(1:im,1:im);
L21 = L(im+1:im+il,1:im);
L22 = L(im+1:im+il,im+1:im+il);
Oi = L22;

[UU,SS,~] = svd(Oi); 

U1 = UU(:,1:n);
Ob_est_i = U1*sqrtm(SS(1:n,1:n));
C = Ob_est_i(1:l,1:n);
A = pinv(Ob_est_i(1:l*(i-1),1:n))*Ob_est_i(l+1:l*i,1:n);
U2 = UU(:,n+1:size(UU',1))';
Z = U2*L21/L11;

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

function [A,B,C,D,SS]=moesp_po(u,y,n,k)
% Syntax
% [A,B,C,D,SS]=moesp_po(u,y,n,k);
%
% Description
% The moesp_po algorithm estimates the model matrices in state-space 
% representation for systems with colored process and measurement noise. 
% It can also be used to determine the system order.
%
% Input Data
% u -> measured input data matrix of the system;
% y -> measured output data matrix of the system;
% n -> chosen order for the system model;
% k -> number of rows of the Hankel block matrix, defined as
% k = 2(maxord/nusaida). 
% Where maxord is the maximum order the user assumes the
% system has and nusaida is the number of outputs the system has. 
% This relationship is defined empirically by (Van Overschee and De Moor, 1996).
%
% Output Data
% A -> estimated dynamic matrix for the system in state space;
% B -> estimated input matrix for the system in state space;
% C -> estimated output matrix for the system in state space;
% D -> estimated direct transmission matrix for the system in state space;
% SS -> singular values matrix of the oblique projection, used to determine the
% system order;
%
% Implemented by Rodrigo Augusto Ricco: rodrigo.ricco@yahoo.com.br
% Universidade Federal de Minas Gerais - UFMG
% Last modification 11/14/2012 

i=k; 
[l,Ndados] = size(y); [m,Ndados] = size(u); j = Ndados-2*i;
kk = 0;
for k = 1:m:2*i*m-m+1
    kk = kk+1; U(k:k+m-1,:) = u(:,kk:kk+j-1); 
end
kk = 0;
for k = 1:l:2*i*l-l+1
    kk = kk+1;
    Y(k:k+l-1,:) = y(:,kk:kk+j-1); 
end
Uf = U(i*m+1:2*i*m,:); 
Yf = Y(i*l+1:2*i*l,:); 
Up = U(1:i*m,:);       
Yp = Y(1:i*l,:);       
Wp = [Up; Yp];         
H = [Uf; Wp; Yf];      

L = triu(qr([H]'))';
im=i*m;
il=i*l;
L11 = L(1:im,1:im);
L21 = L(im+1:2*im,1:im);
L22 = L(im+1:2*im,im+1:2*im);
L31 = L(2*im+1:2*im+il,1:im);
L32 = L(2*im+1:2*im+il,im+1:2*im);
L33 = L(2*im+1:2*im+il,2*im+1:2*im+il);
L41 = L(2*im+il+1:2*im+2*il,1:im);
L42 = L(2*im+il+1:2*im+2*il,im+1:2*im);
L43 = L(2*im+il+1:2*im+2*il,2*im+1:2*im+il);
L44 = L(2*im+il+1:2*im+2*il,2*im+il+1:2*im+2*il);
R11 = L11;
R21 = [L21;L31];
R22 = [L22 zeros(im,il); L32 L33];
R31 = L41;
R32 = [L42 L43];
R33 = L44;

[UU,SS,VV]=svd([R32]);
U1 = UU(:,1:n);
Oi = U1*sqrtm(SS(1:n,1:n)); 

C = Oi(1:l,1:n);
A = pinv(Oi(1:l*(i-1),1:n))*Oi(l+1:i*l,1:n);

U2 = UU(:,n+1:size(UU',1))';
Z = U2*[R31]/[R11];

mm = l*i-n;
M = zeros(mm*i,m);
LL = zeros(mm*i,l+n);
for h = 1:i
    M((h-1)*mm+1:h*mm,:)=Z(:,(h-1)*m+1:h*m);
    LL((h-1)*mm+1:h*mm,:)=[U2(:,(h-1)*l+1:h*l) U2(:,h*l+1:end)*Oi(1:end-h*l,:)];
end
DB = pinv(LL)*M;
D = DB(1:l,:); 
B = DB(l+1:size(DB,1),:);
end
