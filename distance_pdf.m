%% pdf of distance

clc;
clear all;
fig_num = 11;
%% parameters for simulations
sim = 10^4; % number of simulations
lambda_bs = 0.01;
L = 1e3; % range of PPP
%% System model
alpha_l = 3; %path loss exponent line-of-sight
alpha_nl = 4; %path loss exponent non-line-of-sight
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
%% Simulations
for isim = 1 : sim
    display(['progress: ',num2str(floor(isim/sim*100)),'%']);
    %generate PPP of BS
    BSnum = poissrnd(lambda_bs*L^2);
%     BSLengr=L.*sqrt(rand(1,BSnum));
%     BSangle = 2*pi.*rand(1,BSnum);
    BSXary = unifrnd(-L/2,L/2,1,BSnum); % X-axis of BS
    %BSLengr.*cos(BSangle);
    BSYary = unifrnd(-L/2,L/2,1,BSnum); % Y-axis of BS
    %BSLengr.*sin(BSangle);
    l_si = sqrt((BSXary).^2 + (BSYary).^2); %horizontal distance of BS-origin
    d_si = sqrt(l_si.^2 + h^2);
    angle_si = atan(h./l_si);
    Los_p = 1./(1+varphi.*(exp(-beta.*(180/pi.*angle_si-varphi))));
    Los_p_a = 1.*exp(-c*l_si.^2);
    Los_index = rand(1,BSnum);
    Los_f = Los_index<Los_p; % determination of los propagation;
    Los_f_a = Los_index<Los_p_a; % determination of los propagation_approximation;
    alpha_si = alpha_l.*Los_f + alpha_nl.*(1-Los_f);
    alpha_si_a = alpha_l.*Los_f_a + alpha_nl.*(1-Los_f_a); % actual path loss after los determination
    avg_power = d_si.^(alpha_si);
    avg_power_a = d_si.^(alpha_si_a);
    [x_a,x_b] = min(avg_power); % associated to the BS with largest avg power
    [x_a_a,x_b_a]= min(avg_power_a); % 
    d_typical(isim,1) = l_si(x_b); % distance to the associated BS
    d_typical_a(isim,1) = l_si(x_b_a);
end
%% PDF
spacing_d = max(d_typical)/100;
X_d = [0:spacing_d:max(d_typical)];
% figure(fig_num),hist(d_typical,X_d),hold on;
% [rr_freq x_centers] = hist(d_typical_a,X_d);
% plot(x_centers,rr_freq,'r','Linewidth',2);

%% CDF
cdf_d_sim = f_cdf(d_typical,X_d);
% spacing_d_a = max(d_typical)/20;
% X_d_a = [0:spacing_d_a:max(d_typical)];
cdf_d_sim_2 = f_cdf(d_typical_a,X_d);
cdf_d_a = 1-exp(pi*lambda_bs/c*exp(-c.*X_d.^2)).*exp(-pi*lambda_bs/c*exp(-c*X_d.^(2*alpha_l/alpha_nl)))...
    .*exp(-pi*lambda_bs*X_d.^(2*alpha_l/alpha_nl));
%cdf_d_a_2 = 1-exp(-pi*lambda_bs*X_d.^2); 
% Num  = length(X_d);
% const=zeros(1,Num);
% for n = 1 : Num
% fun = @(x) 1./(1+varphi.*(exp(-beta.*(180/pi.*atan(h./x)-varphi)))).*x;
% const(n) = integral(fun,X_d(n).^(alpha_l/alpha_nl),X_d(n));
% end
% cdf_d_r = 1-exp(-pi*lambda_bs*X_d.^(2*alpha_l/alpha_nl)-2*pi*lambda_bs.*const);
figure(fig_num+1),
plot(X_d(round(linspace(1,100,20))),cdf_d_sim(round(linspace(1,100,20))),'ro','Linewidth',1.5), hold on;
plot(X_d,cdf_d_a,'-b','Linewidth',1.5);
%plot(X_d,cdf_d_a_2,'-k','Linewidth',1.5);
%plot(X_d,cdf_d_r,'-r','Linewidth',1.5);
