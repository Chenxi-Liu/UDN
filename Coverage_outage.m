clc;
clear all;

%Initialization
alpha_nl = 4; %path loss exponent non-line-of-sight
alpha_l = 3; %path loss exponent line-of-sight
L = 1e3; % range of tangular area of PPP
%lambda_bs = 1e-4; % Base station density;
x_u = 0; %location of typical UE;
y_u = 0; % location of typical UE;
N = 1e5; % number of ppp realizations
m = 10; % nakagami-m facotr
beta = 0.14; % LOS prob coefficient
varphi = 11.95; % LOS prob coefficient
% h = 10; 
%% c = 0.001969;
%h = 15;
%%c = 0.0008752;
h = 20;
c = 0.0004923;
lambda_bs = 1e-4;
% simulation
counter = 1;
for r_th = 1e-4 : (1e-3-1e-4)/50 : 1e-3
    r_th
    outage_number = 0;
    outage_number_2 = 0;
    avg_I = 0;
    avg_I_2 = 0;
    for n = 1 : N
        BSnum = poissrnd(lambda_bs*L^2);
        BSXary = unifrnd(-L,L,1,BSnum); % X-axis of BS
        BSYary = unifrnd(-L,L,1,BSnum); % Y-axis of BS
        l_si = sqrt((BSXary).^2 + (BSYary).^2); %horizontal distance of BS-origin
        [a,b] = min(l_si); % b denotes the index of the BS with minimum distance;
        d_si = sqrt(l_si.^2 + h^2);
        angle_si = atan(h./l_si);
        Los_p = 1./(1+varphi.*(exp(-beta.*(180/pi.*angle_si-varphi))));
        Los_p_a = 1.*exp(-c*l_si.^2);
        Ns = BSnum;
        Los_index = rand(1,BSnum);
        Los_f = Los_index<Los_p; % determination of los propagation;
        Los_f_a = Los_index<Los_p_a; % determination of los propagation_approximation;
        h_si_nl = sqrt(0.5)*(randn(1,Ns) + 1i*randn(1,Ns));% rayleigh fading;
        gain_l = gamrnd(m,1/m,[1,Ns]); %nakagami-m gain
        phase_l = unifrnd(-pi,pi,[1,Ns]);% nakagami-m phase
        h_si_l = sqrt(gain_l).*exp(1i*phase_l); %nakagami-m channel
        h_si = h_si_l.*Los_f + h_si_nl.*(1-Los_f); %actual channel after los determination
        h_si_a = h_si_l.*Los_f_a + h_si_nl.*(1-Los_f_a); %actual channel after los determination
        alpha_si = alpha_l.*Los_f + alpha_nl.*(1-Los_f);
        alpha_si_a = alpha_l.*Los_f_a + alpha_nl.*(1-Los_f_a); % actual path loss after los determination
        h_si_typical = h_si(b);
        h_si_typical_a = h_si_a(b);
        alpha_si_typical = alpha_si(b); 
        alpha_si_typical_a = alpha_si_a(b);
        d_si_typical = d_si(b);
        h_si(b)=[];
        h_si_a(b)=[];
        h_si_interf = h_si;
        h_si_interf_a = h_si_a;
        alpha_si(b) = [];
        alpha_si_a(b) = [];
        alpha_si_interf = alpha_si;
        alpha_si_interf_a = alpha_si_a;
        d_si(b)=[];
        d_si_interf = d_si;
        Typical_signal = abs(h_si_typical)^2*d_si_typical^(-alpha_si_typical);
        Typical_signal_a = abs(h_si_typical_a)^2*d_si_typical^(-alpha_si_typical_a);
        Typical_interf = abs(h_si_interf).^2.*d_si_interf.^(-alpha_si_interf);
        Typical_interf_a = abs(h_si_interf_a).^2.*d_si_interf.^(-alpha_si_interf_a);
        gamma_typical = sum(Typical_interf);
        %Typical_signal/
        gamma_typical_a = sum(Typical_interf_a);
        avg_I = avg_I + gamma_typical_a;
        avg_I_2 = avg_I_2 + gamma_typical_a^2;
        %Typical_signal/
        if gamma_typical<r_th
            outage_number = outage_number + 1;
        else
            outage_number = outage_number;
        end
        if gamma_typical_a<r_th
            outage_number_2 = outage_number_2 + 1;
        else
            outage_number_2 = outage_number_2;
        end
    end
    mean_I(counter) = avg_I/N;
    var_I(counter) = avg_I_2/N - mean_I(counter)^2;
    outage_sim(counter) = outage_number/N;
    outage_sim_a(counter) = outage_number_2/N;
    x_sim(counter) = r_th;
    counter = counter + 1;
end
hold on;
semilogy(x_sim,outage_sim,'r-');
semilogy(x_sim,outage_sim_a,'k-');

