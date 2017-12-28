clc;
clear all;

L = 1e3; % range of tangular area of PPP
N = 1e4; % number of ppp realizations
lambda_bs = 1e-2;
eta = 4;
avg_I = 0;
    avg_I_2 = 0;
    for n = 1 : N
        display(['progress: ',num2str(floor(n/N*100)),'%']);
        BSnum = poissrnd(lambda_bs*L^2);
        BSXary = unifrnd(-L,L,1,BSnum); % X-axis of BS
        BSYary = unifrnd(-L,L,1,BSnum); % Y-axis of BS
        l_si = sqrt((BSXary).^2 + (BSYary).^2); %horizontal distance of BS-origin
        [a,b] = min(l_si);
        l_si(b)=[];
        l_si_interf = l_si;
        Typical_interf(n) = sum(l_si_interf.^(-eta));
        Typical_interf_a(n) = sum(l_si_interf.^(-2*eta));   
    end