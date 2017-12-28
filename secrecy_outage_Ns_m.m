% Laplace transform of Interference
clc;
clear all;
fig_num = 11;
%% parameters for simulations
sim = 10^3; % number of simulations
lambda_bs = 0.01;
lambda_e = 0.001;
L = 1e3; % range of PPP
Ns = 3; % rnumber of antennas at BS
m = 1; % nakagami-m facotr
x_u = 0; %location of typical UE;
y_u = 0; % location of typical UE;
%% System model
alpha_l = 3; %path loss exponent LoS
alpha_nl = 4; % path loss exponent NLoS
beta = 0.14; % LoS prob coefficient
varphi = 11.95; % LoS prob coefficient
%% Approximation parameter
% h = 10; 
% c = 0.001969;
% h = 15;
% c = 0.0008752;
h = 20;
c = 0.0004923;
d_typical = zeros(sim,1);
d_typical_a = zeros(sim,1);
%s = 1;
%% Simulations
counter = 1;
threshold = linspace(1e-3,1e-1,10);
x_threshold = length(threshold);
outage_sim = zeros(1,x_threshold);
x_sim = zeros(1,x_threshold);
for n = 1 : x_threshold
    display(['progress_sim: ',num2str(floor(n/x_threshold*100)),'%']);
    r_th = threshold(n);
    outage_number = 0;
    outage_number_a = 0;
    for isim = 1 : sim
        %isim
        %display(['progress_inner: ',num2str(floor(isim/sim*100)),'%']);
        %generate PPP of BS
        BSnum = poissrnd(lambda_bs*L^2);
        BSXary = unifrnd(-L/2,L/2,1,BSnum); % X-axis of BS
        BSYary = unifrnd(-L/2,L/2,1,BSnum); % Y-axis of BS
        % generate ppp of Eve
        Enum = poissrnd(lambda_e*L^2);
        EXary = unifrnd(-L/2,L/2,1,Enum); % X-axis of Eavesdroppers
        EYary = unifrnd(-L/2,L/2,1,Enum); % X-axis of Eavesdroppers
        % find the typical BS
        l_si = sqrt((BSXary).^2 + (BSYary).^2); %horizontal distance of BS-origin
        d_si = sqrt(l_si.^2 + h^2); % 3-D distance
        d_si_re  = reshape(d_si,1,1,BSnum);
        angle_si = atan(h./l_si); % elevation angle
        Los_p = 1./(1+varphi.*(exp(-beta.*(180/pi.*angle_si-varphi)))); % LOS possibility
        Los_index = rand(1,BSnum);
        Los_f = Los_index<Los_p; % determination of los propagation;
        Los_f_re = reshape(Los_f,1,1,BSnum);
        h_si_nl = sqrt(0.5)*(randn(1,Ns,BSnum) + 1i*randn(1,Ns,BSnum));% rayleigh fading;
        gain_l = gamrnd(m,1/m,[1,Ns,BSnum]); %nakagami-m gain
        phase_l = unifrnd(-pi,pi,[1,Ns,BSnum]);% nakagami-m phase
        h_si_l = sqrt(gain_l).*exp(1i*phase_l); %nakagami-m channel
        h_si = h_si_l.*Los_f_re + h_si_nl.*(1-Los_f_re);
        alpha_si = alpha_l.*Los_f_re + alpha_nl.*(1-Los_f_re);% actual path loss after los determination
        avg_power = d_si_re.^(-alpha_si);
        [x_a,x_b] = max(avg_power);
        h_si_typical = h_si(:,:,x_b);
        w = h_si_typical/norm(h_si_typical,'fro');
        BSXary_typical  = BSXary(x_b);
        BSXary(x_b)=[];
        BSXary_rest = BSXary;
        BSXary_re = [BSXary_typical BSXary_rest];
        BSYary_typical  = BSYary(x_b);
        BSYary(x_b)=[];
        BSYary_rest = BSYary;
        BSYary_re = [BSYary_typical BSYary_rest];
        gamma_se = zeros(1,Enum);
        for n_eve = 1: Enum
            EXary_t = EXary(n_eve);
            EYary_t = EYary(n_eve);
            l_se = sqrt((BSXary_re-EXary_t).^2 + (BSYary_re-EYary_t).^2); %2-D distance from BSs to Eve.
            d_se = sqrt(l_se.^2 + h^2); % 3-D distance from BSs to Eve.
            d_se_re  = reshape(d_se,1,1,BSnum);
            angle_se = atan(h./l_se); % elevation angle
            Los_p_e = 1./(1+varphi.*(exp(-beta.*(180/pi.*angle_se-varphi)))); % LOS possibility
            Los_index_e = rand(1,BSnum);
            Los_f_e = Los_index_e<Los_p_e; % determination of los propagation;
            Los_f_e_re = reshape(Los_f,1,1,BSnum);
            h_se_nl = sqrt(0.5)*(randn(1,Ns,BSnum) + 1i*randn(1,Ns,BSnum));
            gain_e_l = gamrnd(m,1/m,[1,Ns,BSnum]); %nakagami-m gain
            phase_e_l = unifrnd(-pi,pi,[1,Ns,BSnum]);% nakagami-m phase
            h_se_l = sqrt(gain_e_l).*exp(1i*phase_e_l); %nakagami-m channel
            h_se = h_se_l.*Los_f_e_re + h_se_nl.*(1-Los_f_e_re);
            alpha_se = alpha_l.*Los_f_e_re + alpha_nl.*(1-Los_f_e_re);
            avg_e_power = d_se_re.^(-alpha_se);
            gamma_se_typical = abs(h_se(:,:,1)*w')^2*avg_e_power(:,:,1);
            gamma_se_interf = sum(abs(sum(h_se(:,:,2:BSnum).*conj(w))).^2.*avg_e_power(:,:,2:BSnum));
            gamma_se(n_eve) = gamma_se_typical/gamma_se_interf;
        end
        gamma_se_real=max(gamma_se);
        if gamma_se_real>r_th
            outage_number = outage_number + 1;
        else
            outage_number = outage_number;
        end
    end
    outage_sim(n) = outage_number/sim;
    x_sim(n) = r_th;
end

figure(fig_num);

plot(x_sim,outage_sim,'ro');
hold on;
%semilogy(x_sim,outage_sim_a,'k-');
%% Analysis
% threshold = linspace(1e-1,1,20);
% x_threshold = length(threshold);
% outage_ana = zeros(1,x_threshold);
% outage_ana_app = zeros(1,x_threshold);
% x_ana = zeros(1,x_threshold);
% for n = 1 : x_threshold
%     display(['progress_ana: ',num2str(floor(n/x_threshold*100)),'%']);
%     r_th = threshold(n);
%     outage_ana(n) = MyLaplace(r_th,Ns,m,c,h,lambda_bs,alpha_l,alpha_nl,varphi,beta);
%     outage_ana_app(n) = MyLaplace_app(r_th,Ns,m,c,h,lambda_bs,alpha_l,alpha_nl,varphi,beta);
%     x_ana(n) = r_th;
% end
% 
% plot(x_ana,outage_ana,'k-');
% plot(x_ana,outage_ana_app,'k--');
% hold on;