% Code for example 2 of the paper
% It takes a long time to run due to R = 334. Changing its value
% uses fewer Monte Carlo simulations (3xR), but makes execution faster.
clc
clear
close all
rng(2)

% Default case
Q0 = 10;
Nop0 = 5;
SNR0 = 20;
R = 334;

% Sweep vectors
Q_vals = [4 6 8 10 12 14 16 20 24];
Nop_vals = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 20 24 28 32];
SNR_vals = [1 2.5 5 7.5 10 12.5 15 20 25 30 35 40];

rmse_Q_M1 = zeros(length(Q_vals), R, 2);
rmse_Q_M2 = zeros(length(Q_vals), R, 2);
rmse_Q_M3 = zeros(length(Q_vals), R, 2);
rmse_Q_M = zeros(length(Q_vals), 3*R, 2);
counter = 0;
st = 0;

for iq = 1:length(Q_vals)
    tic
    rmse_Q_M1(iq,:,:) = run_case(Q_vals(iq), Nop0, SNR0,R);
    counter = counter + 1;
    tp = toc;
    st = tp + st;
    tr = (3*length(Q_vals)-counter)*st/counter;
    tr = seconds(tr);
    tr.Format = 'hh:mm:ss';
    fprintf('approx %s remaining\n',tr)
    
    tic
    rmse_Q_M2(iq,:,:) = run_case(Q_vals(iq), Nop0, SNR0,R);
    counter = counter + 1;
    tp = toc;
    st = tp + st;
    tr = (3*length(Q_vals)-counter)*st/counter;
    tr = seconds(tr);
    tr.Format = 'hh:mm:ss';
    fprintf('approx %s remaining\n',tr)
    
    tic
    rmse_Q_M3(iq,:,:) = run_case(Q_vals(iq), Nop0, SNR0,R);
    counter = counter + 1;
    tp = toc;
    st = tp + st;
    tr = (3*length(Q_vals)-counter)*st/counter;
    tr = seconds(tr);
    tr.Format = 'hh:mm:ss';
    fprintf('approx %s remaining\n',tr)
    
    rmse_Q_M(iq,:,:) = [rmse_Q_M1(iq,:,:),rmse_Q_M2(iq,:,:),rmse_Q_M3(iq,:,:)];
    fprintf('Q=%d, N_weights=%d, SNR=%d -> RMSE1=%.4f, RMSE2=%.4f\n', Q_vals(iq), Nop0, SNR0, mean(rmse_Q_M(iq,:,1)), mean(rmse_Q_M(iq,:,2)));   
end

%% Calculate medians for each SNR
median1 = median(rmse_Q_M(:,:,1), 2);
median2 = median(rmse_Q_M(:,:,2), 2);

% Figure with boxchart and Costm median (manual overwrite)
figure('Color','w','Position',[50 290 900 350])

subplot(1,2,1);
hold on;
x1 = repelem(1:length(Q_vals), 3*R);
y1 = rmse_Q_M(:,:,1).';
y1 = y1(:);
b1 = boxchart(x1, y1, 'BoxWidth', 0.9);
b1.MarkerStyle = '.';
b1.MarkerColorMode = 'manual';
b1.MarkerColor = [0.4 0.5 1];
b1.MarkerSize = 10;

% Overwrite median: plot horizontal lines
for i = 1:length(Q_vals)
    x_line = [i - 0.45, i + 0.45];
    y_line = [median1(i), median1(i)];
    plot(x_line, y_line, 'Color', 'k', 'LineWidth', 2.5);
end

set(gca, 'XTick', 1:length(Q_vals));
xticklabels(Q_vals);  
xlabel('$N_N$','interpreter','latex');
ylabel('NRMSE($y_1$)','interpreter','latex');
grid on;
box on;
set(gca, 'FontName', 'Times', 'FontSize', 16);

% Means as red dots
mean_val1 = mean(rmse_Q_M(:, :, 1), 2);
plot(1:length(Q_vals), mean_val1, 'rs-', 'MarkerSize',5, 'LineWidth', 1, 'MarkerFaceColor', 'r');
% legend('Mean RMSE', 'Median', 'Location', 'best');

subplot(1,2,2);
hold on;
x2 = repelem(1:length(Q_vals), 3*R);
y2 = rmse_Q_M(:,:,2).';
y2 = y2(:);
b2 = boxchart(x2, y2, 'BoxWidth', 0.9);
b2.MarkerStyle = '.';
b2.MarkerColorMode = 'manual';
b2.MarkerColor = [0.4 0.5 1];
b2.MarkerSize = 10;

for i = 1:length(Q_vals)
    x_line = [i - 0.45, i + 0.45];
    y_line = [median2(i), median2(i)];
    plot(x_line, y_line, 'Color', 'k', 'LineWidth', 2.5);
end

ylim([0 0.3])
set(gca, 'XTick', 1:length(Q_vals));
xticklabels(Q_vals);
xlabel('$N_N$','interpreter','latex');
ylabel('NRMSE($y_2$)','interpreter','latex');
grid on;
box on;
set(gca, 'FontName', 'Times', 'FontSize', 16);

mean_val2 = mean(rmse_Q_M(:, :, 2), 2);
plot(1:length(Q_vals), mean_val2, 'rs-', 'MarkerSize',5, 'LineWidth', 1, 'MarkerFaceColor', 'r');
% legend('Mean RMSE', 'Median', 'Location', 'best');

set(gcf, 'Color', 'w', 'Position', [50 290 900 350]);

%%
rmse_Nop_M1 = zeros(length(Nop_vals), R, 2);
rmse_Nop_M2 = zeros(length(Nop_vals), R, 2);
rmse_Nop_M3 = zeros(length(Nop_vals), R, 2);
rmse_Nop_M = zeros(length(Nop_vals), 3*R, 2);

for inop = 1:length(Nop_vals)
    rmse_Nop_M1(inop,:,:) = run_case(Q0, Nop_vals(inop), SNR0,R);
    rmse_Nop_M2(inop,:,:) = run_case(Q0, Nop_vals(inop), SNR0,R);
    rmse_Nop_M3(inop,:,:) = run_case(Q0, Nop_vals(inop), SNR0,R);
    rmse_Nop_M(inop,:,:) = [rmse_Nop_M1(inop,:,:),rmse_Nop_M2(inop,:,:),rmse_Nop_M3(inop,:,:)];
    fprintf('Q=%d, N_weights=%d, SNR=%d -> RMSE1=%.4f, RMSE2=%.4f\n', Q0, Nop_vals(inop), SNR0, mean(rmse_Nop_M(inop,:,1)), mean(rmse_Nop_M(inop,:,2)));   
end

%% Calculate medians for each SNR
median1 = median(rmse_Nop_M(:,:,1), 2);
median2 = median(rmse_Nop_M(:,:,2), 2);

% Figure with boxchart and Costm median (manual overwrite)
figure

subplot(1,2,1);
hold on;
x1 = repelem(1:length(Nop_vals), 3*R);
y1 = rmse_Nop_M(:,:,1).';
y1 = y1(:);
b1 = boxchart(x1, y1, 'BoxWidth', 0.9);
b1.MarkerStyle = '.';
b1.MarkerColorMode = 'manual';
b1.MarkerColor = [0.4 0.5 1];
b1.MarkerSize = 10;

% Overwrite median: plot horizontal lines
for i = 1:length(Nop_vals)
    x_line = [i - 0.45, i + 0.45];
    y_line = [median1(i), median1(i)];
    plot(x_line, y_line, 'Color', 'k', 'LineWidth', 2.5);
end

ylim([0 0.3123])
xlim([0 21])
set(gca, 'XTick', 1:length(Nop_vals));
xticklabels(Nop_vals);
xlabel('$N_P$','interpreter','latex');
ylabel('NRMSE($y_1$)','interpreter','latex');
grid on;
box on;
set(gca, 'FontName', 'Times', 'FontSize', 16);

% Means as red dots
mean_val1 = mean(rmse_Nop_M(:, :, 1), 2);
plot(1:length(Nop_vals), mean_val1, 'rs-', 'MarkerSize',5, 'LineWidth', 1, 'MarkerFaceColor', 'r');
% legend('Mean RMSE', 'Median', 'Location', 'best');

subplot(1,2,2);
hold on;
x2 = repelem(1:length(Nop_vals), 3*R);
y2 = rmse_Nop_M(:,:,2).';
y2 = y2(:);
b2 = boxchart(x2, y2, 'BoxWidth', 0.9);
b2.MarkerStyle = '.';
b2.MarkerColorMode = 'manual';
b2.MarkerColor = [0.4 0.5 1];
b2.MarkerSize = 10;

for i = 1:length(Nop_vals)
    x_line = [i - 0.45, i + 0.45];
    y_line = [median2(i), median2(i)];
    plot(x_line, y_line, 'Color', 'k', 'LineWidth', 2.5);
end

xlim([0 21])
set(gca, 'XTick', 1:length(Nop_vals));
xticklabels(Nop_vals);
xlabel('$N_P$','interpreter','latex');
ylabel('NRMSE($y_2$)','interpreter','latex');
grid on;
box on;
set(gca, 'FontName', 'Times', 'FontSize', 16);

mean_val2 = mean(rmse_Nop_M(:, :, 2), 2);
plot(1:length(Nop_vals), mean_val2, 'rs-', 'MarkerSize',5, 'LineWidth', 1, 'MarkerFaceColor', 'r');
% legend('Mean RMSE', 'Median', 'Location', 'best');

set(gcf, 'Color', 'w', 'Position', [50 290 1200 350]);

%%
rmse_SNR_M1 = zeros(length(SNR_vals), R, 2);
rmse_SNR_M2 = zeros(length(SNR_vals), R, 2);
rmse_SNR_M3 = zeros(length(SNR_vals), R, 2);
rmse_SNR_M = zeros(length(SNR_vals), 3*R, 2);

for isnr = 1:length(SNR_vals)
    rmse_SNR_M1(isnr,:,:) = run_case(Q0, Nop0, SNR_vals(isnr),R);
    rmse_SNR_M2(isnr,:,:) = run_case(Q0, Nop0, SNR_vals(isnr),R);
    rmse_SNR_M3(isnr,:,:) = run_case(Q0, Nop0, SNR_vals(isnr),R);
    rmse_SNR_M(isnr,:,:) = [rmse_SNR_M1(isnr,:,:),rmse_SNR_M2(isnr,:,:),rmse_SNR_M3(isnr,:,:)];
    fprintf('Q=%d, N_weights=%d, SNR=%d -> RMSE1=%.4f, RMSE2=%.4f\n', Q0, Nop0, SNR_vals(isnr), mean(rmse_SNR_M(isnr,:,1)), mean(rmse_SNR_M(isnr,:,2)));   
end

%% Calculate medians for each SNR
median1 = median(rmse_SNR_M(:,:,1), 2);
median2 = median(rmse_SNR_M(:,:,2), 2);

% Figure with boxchart and Costm median (manual overwrite)
figure('Color','w','Position',[50 290 900 350])

subplot(1,2,1);
hold on;
x1 = repelem(1:length(SNR_vals), 3*R);
y1 = rmse_SNR_M(:,:,1).';
y1 = y1(:);
b1 = boxchart(x1, y1, 'BoxWidth', 0.9);
b1.MarkerStyle = '.';
b1.MarkerColorMode = 'manual';
b1.MarkerColor = [0.4 0.5 1];
b1.MarkerSize = 10;

% Overwrite median: plot horizontal lines
for i = 1:length(SNR_vals)
    x_line = [i - 0.45, i + 0.45];
    y_line = [median1(i), median1(i)];
    plot(x_line, y_line, 'Color', 'k', 'LineWidth', 2.5);
end

xlim([0 13])
set(gca, 'XTick', 1:length(SNR_vals));
xticklabels(SNR_vals);
xlabel('SNR');
ylabel('NRMSE($y_1$)','interpreter','latex');
grid on;
box on;
set(gca, 'FontName', 'Times', 'FontSize', 16);

% Means as red dots
mean_val1 = mean(rmse_SNR_M(:, :, 1), 2);
plot(1:length(SNR_vals), mean_val1, 'rs-', 'MarkerSize',5, 'LineWidth', 1, 'MarkerFaceColor', 'r');
% legend('Mean RMSE', 'Median', 'Location', 'best');

subplot(1,2,2);
hold on;
x2 = repelem(1:length(SNR_vals), 3*R);
y2 = rmse_SNR_M(:,:,2).';
y2 = y2(:);
b2 = boxchart(x2, y2, 'BoxWidth', 0.9);
b2.MarkerStyle = '.';
b2.MarkerColorMode = 'manual';
b2.MarkerColor = [0.4 0.5 1];
b2.MarkerSize = 10;

for i = 1:length(SNR_vals)
    x_line = [i - 0.45, i + 0.45];
    y_line = [median2(i), median2(i)];
    plot(x_line, y_line, 'Color', 'k', 'LineWidth', 2.5);
end

xlim([0 13])
set(gca, 'XTick', 1:length(SNR_vals));
xticklabels(SNR_vals);
xlabel('SNR');
ylabel('NRMSE($y_2$)','interpreter','latex');
grid on;
box on;
set(gca, 'FontName', 'Times', 'FontSize', 16);

mean_val2 = mean(rmse_SNR_M(:, :, 2), 2);
plot(1:length(SNR_vals), mean_val2, 'rs-', 'MarkerSize',5, 'LineWidth', 1, 'MarkerFaceColor', 'r');
% legend('Mean RMSE', 'Median', 'Location', 'best');

set(gcf, 'Color', 'w', 'Position', [50 290 900 350]);
set(gcf, 'PaperPositionMode', 'auto');

%%
save data_example1

%% Functions 
function rmse_final = run_case(Q, N_weights, SNR,R)
Ts = 0.001;

%% Fixed process parameters and dynamic input (kept from original)
Amp_dyn = 2.5;
f1td = pi; f2td = 6*pi;
f1td2 = 4; f2td2 = 21;
Nd = ceil(40/min(f1td,f2td)/Ts);
t = 0:Ts:Nd*Ts - Ts;
n_transient = ceil(1/min(f1td,f2td)/Ts);

utd = [Amp_dyn*sin(2*pi*f1td.*t)+Amp_dyn*sin(2*pi*f2td.*t);
       Amp_dyn*sin(2*pi*f1td2.*t)+Amp_dyn*sin(2*pi*f2td2.*t)];

% Generates ytd and vtd_o only for reference if needed, focus is the sweep
Maxutd = max(abs(utd'), [], 1);
ytd = Process(length(utd),utd,Maxutd,Ts,2);

%% PSO Parameters (fixed)
alpha = 0.7;
Kmax = 150;
P = 60;
b = ceil(0.05*P);
phi1 = 2.05;
phi2 = 2.05;
phi = phi1+phi2;
chi = 2/abs(2-phi-sqrt(phi^2-4*phi));
c1 = phi1;
c2 = phi2;

%% SWEEP DEFINITION
% Variable parameters for PSO
VarMin = 0;
VarMax = 4;

%% Dynamics validation
f1d = 10;
N_val = ceil(3/f1d/Ts);
Amp = 5;
t = 0:Ts:N_val*Ts - Ts;

uvd = [Amp*sin(2*pi*f1d.*t + pi); Amp*sin(2*pi*f1d.*t + pi)];
Maxuvd = max(abs(uvd'), [], 1);
yvd = Process(N_val,uvd,Maxuvd,Ts,2);
 
Nv = N_weights + 5;
VarSize = [2 Nv];
VelMax = 0.05*VarMax;

L = 200;
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
n_transienth = Q;
uta1 = uta(1,1:2*Q);
uta2 = uta(2,length(uta)/2+1:length(uta)/2+2*Q);
u = [uta1; uta2];
Maxu = max(abs(u'), [], 1);
t = 0:Ts:length(u)*Ts - Ts;
Ne = length(t);

%% 2. PSO execution for hysteresis identification
rmse_final = zeros(R,2);
for realiza = 1:R
    
    % Add noise to the static part according to original logic
    yo = add_noise_snr_per_channel(yt, SNR);
    l = L;
    FD = 100;
    for k = 1:length(uta)
        yo(:,k) = sum(yo(:,l-FD:l),2)/FD;
        l = l + L;
    end
    yta = yo(:,1:length(uta));
    yta1 = yta(1,1:2*Q);
    yta2 = yta(2,length(uta)/2+1:length(uta)/2+2*Q);
    y_hyster = [yta1; yta2];
    
    %% PSO Initialization
    % Iteration parameters
    Theta = zeros([VarSize P]);
    V = zeros([VarSize P]);
    F = inf(VarSize(1),P);
    
    % Global parameters
    Thetag = zeros(VarSize);
    Fg = inf(VarSize(1),1);
    
    % Personal parameters
    Thetap = zeros([VarSize P]);
    Fp = inf(VarSize(1),P);
    
    for i = 1:P
        Theta(:,:,i) = unifrnd(0,VarMax,VarSize);
        Theta(:,N_weights+1:end,i) = min(Theta(:,N_weights+1:end,i),1);
        V(:,:,i) = unifrnd(-VelMax,VelMax,VarSize);
        
        % Evaluation
        F(:,i) = Cost(Ne,N_weights,Theta(:,:,i),u,Maxu,y_hyster,n_transienth,Ts);
        
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
            F(:,i) = Cost(Ne,N_weights,Theta(:,:,i),u,Maxu,y_hyster,n_transienth,Ts);
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
%                     disp(['Iteration ' num2str(k) ': Best Cost ' num2str(l)  ' = ' num2str(Fg(l))]);
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
    
    %% 3. Dynamics identification and Validation
    Thetag_final = Thetag;
    vtd = Hysteresis(length(utd),N_weights,Thetag_final,utd,Maxutd,Ts,2);
    ytdr = add_noise_snr_per_channel(ytd, 25);
    [A, B, C, D, ~] = moesp(vtd(:,n_transient:end),ytdr(:,n_transient:end),2,2);
    
    y_est = swarm_data(length(uvd),N_weights,Thetag_final,uvd,Maxuvd,A,B,C,D,Ts);
    rmseval = sqrt(sum((y_est - yvd).^2,2)./length(yvd))./max(yvd,[],2);
    rmse_final(realiza,:,:) = rmseval;
end
end

function yo = add_noise_snr_per_channel(yt, val_SNR)
    yo = zeros(size(yt));
    for i = 1:size(yt,1)
        Ps = mean(yt(i,:).^2);
        Pn = Ps / (10^(val_SNR/10));
        yo(i,:) = yt(i,:) + sqrt(Pn) * randn(size(yt,2),1).';
    end
end

function v_est = Hysteresis(N_dados,N_weights,pesos,ute,Maxute,Ts,type_sim)
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
    if k == 1 || (k == ceil(N_dados/2)+1 && type_sim == 1)
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
        if k == 1 || (k == ceil(N_dados/2)+1 && type_sim == 1)
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

function y = Process(N_dados,u,Maxu,Ts,type_sim)
A = [0.8 0.1; 
    -0.1 0.8];
B = [0.5    0.1;
     0.2    0.5];
C = eye(2)-A;
    
D = [0 0; 0 0];
theta = [0.1, 0.9, 1, 0.1, 0.9, 0.3, 0.2, 0.6, 0.8, 0.4, 0.9;
         0.1, 1, 0.9, 0.1, 0.7, 0.1, 0.8, 0.4, 0.1, 0.5, 0.6];
theta_hat = [0.9, 0.35, 1, 0.25;
               1,  0.4, 1, 0.4];
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
    if k == 1 || (k == ceil(N_dados/2)+1 && type_sim == 1)
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
        if k == 1 || (k == ceil(N_dados/2)+1 && type_sim == 1)
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

i = k; % i=number of rows;

[l, ~] = size(y);
[m, N_data] = size(u);
j = N_data - 2*i;
kk = 0;

% Assembling Hankel block matrices
for k = 1:m:2*i*m-m+1
    kk = kk+1;
    U(k:k+m-1,:) = u(:,kk:kk+j-1); % U data matrix
end
kk = 0;
for k = 1:l:2*i*l-l+1
    kk = kk+1;
    Y(k:k+l-1,:) = y(:,kk:kk+j-1); % Y data matrix
end

Uf = U(i*m+1:2*i*m,:); % future inputs
Yf = Y(i*l+1:2*i*l,:); % future outputs

H = [Uf; Yf]; % final data matrix

% LQ Decomposition
L = triu(qr(H'))'; 
im=i*m;
il=i*l;
L11 = L(1:im,1:im);
L21 = L(im+1:im+il,1:im);
L22 = L(im+1:im+il,im+1:im+il);
Oi = L22;

% Singular value decomposition
[UU,SS,~] = svd(Oi); % X = UU*SS*VV' is how svd decomposes the matrix.
% SS contains the singular values in descending order

U1 = UU(:,1:n);
Ob_est_i = U1*sqrtm(SS(1:n,1:n));
C = Ob_est_i(1:l,1:n);
A = pinv(Ob_est_i(1:l*(i-1),1:n))*Ob_est_i(l+1:l*i,1:n);
U2 = UU(:,n+1:size(UU',1))';
Z = U2*L21/L11;

% From this point on, all algorithms are the same
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
