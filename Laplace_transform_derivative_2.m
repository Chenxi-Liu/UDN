% Laplace transform of Interference
clc;
clear all;
fig_num = 11;
%% parameters for simulations
sim = 10^4; % number of simulations
lambda_bs = 0.01;
L = 1e3; % range of PPP
%Ns = 1; % number of antennas at BS
m = 2; % nakagami-m facotr
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
s = 1;
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
    for isim = 1 : sim
        %generate PPP of BS
        BSnum = poissrnd(lambda_bs*L^2);
        BSXary = unifrnd(-L/2,L/2,1,BSnum); % X-axis of BS
        BSYary = unifrnd(-L/2,L/2,1,BSnum); % Y-axis of BS
        l_si = sqrt((BSXary).^2 + (BSYary).^2); %horizontal distance of BS-origin
        d_si = sqrt(l_si.^2 + h^2); % 3-D distance
        angle_si = atan(h./l_si); % elevation angle
        Los_p = 1./(1+varphi.*(exp(-beta.*(180/pi.*angle_si-varphi)))); % LOS possibility
        Los_p_a = 1.*exp(-c*l_si.^2);
        Los_index = rand(1,BSnum);
        Los_f = Los_index<Los_p; % determination of los propagation;
        Los_f_a = Los_index<Los_p_a;
        h_si_nl = sqrt(0.5)*(randn(1,BSnum) + 1i*randn(1,BSnum));% rayleigh fading;
        gain_l = gamrnd(m,1/m,[1,BSnum]); %nakagami-m gain
        phase_l = unifrnd(-pi,pi,[1,BSnum]);% nakagami-m phase
        h_si_l = sqrt(gain_l).*exp(1i*phase_l); %nakagami-m channel
        h_si = h_si_l.*Los_f + h_si_nl.*(1-Los_f); %actual channel after los determination
        h_si_a = h_si_l.*Los_f_a + h_si_nl.*(1-Los_f_a);
        alpha_si = alpha_l.*Los_f + alpha_nl.*(1-Los_f);% actual path loss after los determination
        alpha_si_a = alpha_l.*Los_f_a + alpha_nl.*(1-Los_f_a); 
        avg_power = d_si.^(-alpha_si);
        avg_power_a = d_si.^(-alpha_si_a);
        [x_a,x_b] = max(avg_power);
        [x_a_a,x_b_a] = max(avg_power_a);
        h_si_typical = h_si(x_b);
        h_si_typical_a = h_si(x_b_a);
        avg_power_typical = avg_power(x_b); 
        avg_power_typical_a = avg_power(x_b_a);
        h_si(x_b)=[];
        h_si_interf = h_si;
        h_si_a(x_b_a)=[];
        h_si_interf_a = h_si_a;
        avg_power(x_b)=[];
        avg_power_interf = avg_power;
        avg_power_a(x_b_a)=[];
        avg_power_interf_a = avg_power_a;
        Typical_signal = abs(h_si_typical)^2*avg_power_typical;
        Typical_interf = sum(abs(h_si_interf).^2.*avg_power_interf);
        Typical_signal_a = abs(h_si_typical_a)^2*avg_power_typical_a;
        Typical_interf_a = sum(abs(h_si_interf_a).^2.*avg_power_interf_a);
        gamma_typical = Typical_signal/Typical_interf;
        gamma_typical_a = Typical_signal_a/Typical_interf_a;
        if gamma_typical<r_th
            outage_number = outage_number + 1;
        else
            outage_number = outage_number;
        end
        
        if gamma_typical_a<r_th
            outage_number_a = outage_number_a + 1;
        else
            outage_number_a = outage_number_a;
        end
    end
    outage_sim(n) = outage_number/sim;
    outage_sim_a(n) = outage_number_a/sim;
    x_sim(n) = r_th;
end

figure(fig_num);

semilogy(x_sim,outage_sim,'r-');
hold on;
semilogy(x_sim,outage_sim_a,'k-');
%% Analysis
outage_ana = zeros(1,x_threshold);
x_ana = zeros(1,x_threshold);
for n = 1 : x_threshold
    display(['progress_ana: ',num2str(floor(n/x_threshold*100)),'%']);
    r_th = threshold(n);
    part_1 = 1-integral(@(x)(1-exp(-c.*x.^2)).*2*pi*lambda_bs.*exp(-pi*lambda_bs.*x.^(2*alpha_l/alpha_nl))...
        .*exp(pi*lambda_bs/c.*(exp(-c.*x.^2)-exp(-c.*x.^(2*alpha_l/alpha_nl))))...
        .*(alpha_l/alpha_nl.*x.^(2*alpha_l/alpha_nl-1)-alpha_l/alpha_nl.*x.^(2*alpha_l/alpha_nl-1).*exp(-c.*x.^(2*alpha_l/alpha_nl))+x.*exp(-c.*x.^2))...
        .*exp(-2*pi*lambda_bs.*arrayfun(@(xx)integral(@(y)(1-1./((1+r_th/m.*(xx^2+h^2)^(alpha_nl/2).*(y.^2+h^2).^(-alpha_l/2))).^m).*exp(-c.*y.^2).*y,xx,inf),x))...
        .*exp(-2*pi*lambda_bs.*arrayfun(@(xx)integral(@(z) (1-1./(1+r_th*(xx^2+h^2)^(alpha_nl/2).*(z.^2+h^2).^(-alpha_nl/2))).*(1-exp(-c.*z.^2)).*z,xx^(alpha_l/alpha_nl),inf),x)),0,inf)...
        - integral(@(x)(exp(-c.*x.^2)).*2*pi*lambda_bs.*exp(-pi*lambda_bs.*x.^(2*alpha_l/alpha_nl))...
        .*exp(pi*lambda_bs/c.*(exp(-c.*x.^2)-exp(-c.*x.^(2*alpha_l/alpha_nl))))...
        .*(alpha_l/alpha_nl.*x.^(2*alpha_l/alpha_nl-1)-alpha_l/alpha_nl.*x.^(2*alpha_l/alpha_nl-1).*exp(-c.*x.^(2*alpha_l/alpha_nl))+x.*exp(-c.*x.^2))...
        .*exp(-2*pi*lambda_bs.*arrayfun(@(xx)integral(@(y)(1-1./((1+r_th.*(xx^2+h^2)^(alpha_l/2).*(y.^2+h^2).^(-alpha_l/2))).^m).*exp(-c.*y.^2).*y,xx,inf),x))...
        .*exp(-2*pi*lambda_bs.*arrayfun(@(xx)integral(@(z) (1-1./(1+m*r_th*(xx^2+h^2)^(alpha_l/2).*(z.^2+h^2).^(-alpha_nl/2))).*(1-exp(-c.*z.^2)).*z,xx^(alpha_l/alpha_nl),inf),x))...
        .*(1+(-1)*m*r_th.*(x.^2+h^2).^(alpha_l/2).*(pi*lambda_bs*exp(c*h^2)*(-m).*arrayfun(@(xx) integral(@(y) exp(-c.*y).*(y.^(-alpha_l/2)/m)...
        .*(1+r_th*(xx^2+h^2)^(alpha_l/2).*y.^(-alpha_l/2)).^(-m-1),xx^2+h^2,inf),x)...
        + pi*lambda_bs*exp(c*h^2)*(-1).*arrayfun(@(xx)integral(@(z) (1-exp(-c.*z)).*(z.^(-alpha_nl/2)).*(1+m*r_th*(xx^2+h^2)^(alpha_l/2).*z.^(-alpha_nl/2)).^(-2),xx^(2*alpha_l/alpha_nl)+h^2,inf),x))),0,inf);


    outage_ana(n) = part_1;
    x_ana(n) = r_th;
%     for k = 1 : m-1;
%     end
end




    
    
% (1-exp(-c*x.^2))*
% .*exp(-pi*lambda_bs.*exp(-c*h^2).*arrayfun(@(xx) integral(@(y) (1-1./(1+r_th/m.*(xx^2+h^2)^(alpha_nl/2).*y.^(-alpha_l/2)).^m).*(exp(-c.*xx^2)),xx^2+h^2,1000),x))...
%                 .*exp(-pi*lambda_bs.*exp(-c*h^2).*arrayfun(@(xx) integral(@(z)(1-1./(1+r_th.*(xx.^2+h^2).^(alpha_nl/2).*z.^(-alpha_nl/2))).*(1-exp(-c*xx.^2)),xx.^(2*alpha_l/alpha_nl+h^2),1000),x))
semilogy(x_ana,outage_ana,'ro');
hold on;