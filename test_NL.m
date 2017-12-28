% Laplace transform of Interference
clc;
clear all;
fig_num = 11;
%% parameters for simulations
sim = 10^4; % number of simulations
lambda_bs = 0.01;
%lambda_e = 0.01;
L = 2e3; % range of PPP
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
h = 10; 
c = 0.001969;
% h = 15;
% c = 0.0008752;
% h = 20;
% c = 0.0004923;
d_typical = zeros(sim,1);
d_typical_a = zeros(sim,1);
%s = 1;
%% Simulations
counter = 1;
threshold = linspace(1e-1,1,20);
x_threshold = length(threshold);
outage_sim = zeros(1,x_threshold);
x_sim = zeros(1,x_threshold);
for n = 1 : x_threshold
    display(['progress_sim: ',num2str(floor(n/x_threshold*100)),'%']);
    r_th = threshold(n);
    outage_number = 0;
    outage_number_a = 0;
    BSnum = 39928;
%     BSXary = unifrnd(-L/2,L/2,1,BSnum); % X-axis of BS
%     BSYary = unifrnd(-L/2,L/2,1,BSnum); % Y-axis of BS
    % %         Enum = poissrnd(lambda_e*L^2);
    % %         EXary = unifrnd(-L/2,L/2,1,Enum); % X-axis of Eavesdroppers
    % %         EYary = unifrnd(-L/2,L/2,1,Enum); % X-axis of Eavesdroppers
%     l_si = sqrt((BSXary).^2 + (BSYary).^2); %horizontal distance of BS-origin
%     d_si = sqrt(l_si.^2 + h^2); % 3-D distance
%     angle_si = atan(h./l_si); % elevation angle
    load l_si.mat
    load d_si.mat
    angle_si = atan(h./l_si); % elevation angle
    d_si_re  = reshape(d_si,1,1,BSnum);
    [x_a,x_b]=min(d_si_re)
    for isim = 1 : sim
        %generate PPP of BS
        Los_p = 1./(1+varphi.*(exp(-beta.*(180/pi.*angle_si-varphi)))); % LOS possibility
        Los_p_a = exp(-c*l_si.^2);
        Los_index = rand(1,BSnum);
        Los_f = Los_index<Los_p; % determination of los propagation;
        Los_f_a = Los_index<Los_p_a;
        Los_f_re = reshape(Los_f,1,1,BSnum);
        Los_f_a_re = reshape(Los_f_a,1,1,BSnum);
        h_si_nl = sqrt(0.5)*(randn(1,Ns,BSnum) + 1i*randn(1,Ns,BSnum));% rayleigh fading;
        gain_l = gamrnd(m,1/m,[1,Ns,BSnum]); %nakagami-m gain
        phase_l = unifrnd(-pi,pi,[1,Ns,BSnum]);% nakagami-m phase
        h_si_l = sqrt(gain_l).*exp(1i*phase_l); %nakagami-m channel
        h_si = h_si_l.*Los_f_re + h_si_nl.*(1-Los_f_re);
        h_si_a = h_si_l.*Los_f_a_re + h_si_nl.*(1-Los_f_a_re);
%         for n_bs = 1:BSnum
%             h_si(:,:,n_bs) = h_si_l(:,:,n_bs).*Los_f(n_bs) + h_si_nl(:,:,n_bs).*(1-Los_f(n_bs)); %actual channel after los determination
%             h_si_a(:,:,n_bs) = h_si_l(:,:,n_bs).*Los_f_a(n_bs) + h_si_nl(:,:,n_bs).*(1-Los_f_a(n_bs)); %actual channel after los determination
%         end
        alpha_si = alpha_l.*Los_f_re + alpha_nl.*(1-Los_f_re);% actual path loss after los determination
        alpha_si_a = alpha_l.*Los_f_a_re + alpha_nl.*(1-Los_f_a_re); 
        avg_power = d_si_re.^(-alpha_si);
        avg_power_a = d_si_re.^(-alpha_si_a);
        [x_a,x_b] = max(avg_power);
        [x_a_a,x_b_a] = max(avg_power_a);
        h_si_typical = h_si(:,:,x_b);
        h_si_typical_a = h_si_a(:,:,x_b_a);
        avg_power_typical = avg_power(:,:,x_b); 
        avg_power_typical_a = avg_power_a(:,:,x_b_a);
        h_si(:,:,x_b)=[];
        h_si_interf = h_si;
        h_si_a(:,:,x_b_a)=[];
        h_si_interf_a = h_si_a;
        avg_power(:,:,x_b)=[];
        avg_power_interf = avg_power;
        avg_power_a(:,:,x_b_a)=[];
        avg_power_interf_a = avg_power_a;
        Typical_signal = norm(h_si_typical,'fro').^2*avg_power_typical;
        Typical_signal_a = norm(h_si_typical_a,'fro')^2*avg_power_typical_a;
        
%         for n_bs = 1: BSnum-1
%             Typical_interf(n_bs) = norm(h_si_interf(:,:,n_bs)*h_si_typical'/norm(h_si_typical,'fro'))^2*avg_power_interf(:,:,n_bs);
%             Typical_interf_a(n_bs) = norm(h_si_interf_a(:,:,n_bs)*h_si_typical_a'/norm(h_si_typical_a,'fro'))^2*avg_power_interf_a(:,:,n_bs);
%         end
        if Ns==1
           Typical_interf=sum(abs((h_si_interf.*conj(h_si_typical))./norm(h_si_typical,'fro')).^2.*avg_power_interf);
          % Typical_interf_a = sum(abs((h_si_interf_a.*conj(h_si_typical_a))./norm(h_si_typical_a,'fro')).^2.*avg_power_interf_a);
        else
            Typical_interf=sum(abs(sum(h_si_interf.*conj(h_si_typical))./norm(h_si_typical,'fro')).^2.*avg_power_interf);
           % Typical_interf_a = sum(abs(sum(h_si_interf_a.*conj(h_si_typical_a))./norm(h_si_typical_a,'fro')).^2.*avg_power_interf_a);
        end
%     sum = 0;
%     for p = 1:BSnum-1
%         sumterm = norm(h_si_interf(:,:,p)*h_si_typical'/norm(h_si_typical,'fro'),'fro')^2*avg_power(:,:,p);
%         sum = sum + sumterm;
%     end
%     Typical_interf = sum;
        
        gamma_typical = Typical_signal/Typical_interf;
       % gamma_typical_a = Typical_signal_a/Typical_interf_a;
        if gamma_typical>r_th
            outage_number = outage_number + 1;
        else
            outage_number = outage_number;
        end
        
%         if gamma_typical_a>r_th
%             outage_number_a = outage_number_a + 1;
%         else
%             outage_number_a = outage_number_a;
%         end
    end
    outage_sim(n) = 1-outage_number/sim;
   % outage_sim_a(n) =1- outage_number_a/sim;
    x_sim(n) = r_th;
end

figure(fig_num);

semilogy(x_sim,outage_sim,'ro');
hold on;
%semilogy(x_sim,outage_sim_a,'k-');
%% Analysis
threshold = linspace(1e-1,1,20);
x_threshold = length(threshold);
outage_ana = zeros(1,x_threshold);
outage_ana_app = zeros(1,x_threshold);
x_ana = zeros(1,x_threshold);
for n = 1 : x_threshold
    display(['progress_ana: ',num2str(floor(n/x_threshold*100)),'%']);
    r_th = threshold(n);
    outage_ana(n) = MyLaplace(r_th,Ns,m,c,h,lambda_bs,alpha_l,alpha_nl,varphi,beta);
    outage_ana_app(n) = MyLaplace_app(r_th,Ns,m,c,h,lambda_bs,alpha_l,alpha_nl,varphi,beta);
    x_ana(n) = r_th;
end

semilogy(x_ana,outage_ana,'kp');
semilogy(x_ana,outage_ana_app,'k--');
hold on;