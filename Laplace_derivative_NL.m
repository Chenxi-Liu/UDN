function [r] = Laplace_derivative_NL(r_th,Ns,m,c,h,lambda_bs,alpha_l,alpha_nl,x,varphi,beta)

% r_th: threshold;
% m: Nakagami-m factor;
% c: approximation parameter for LoS probability;
% h: height;
% lambda_bs: density for BS, range should be larger than 0.0019 (c for 10);
% alpha_l: LoS path loss exponent
% alpha_nl: NLoS path loss exponent
n_total = Ns;
sum = 0;
for k = 0 : n_total-1
    if k==0
        sumterm=1;
    else
       [a,b] = Get_parameter_result(k);
       sum_2 = 0;
       for q = 1:b
           m_l = a(:,:,q);
           prodd = 1;
           for l = 1 : k
               prodd_2 = 1;
               prodd_3 = 1;
               for t = 0 : l-1
                   prodd_term_2 = -m-t;
                   prodd_term_3 = -1-t;
                   prodd_2 = prodd_2*prodd_term_2;
                   prodd_3 = prodd_3*prodd_term_3;
               end
               prodd_term = 1/(factorial(m_l(l))*(factorial(l))^(m_l(l)))...
                   .*(2*pi*lambda_bs*prodd_2.*arrayfun(@(xx)integral(@(y)y.*exp(-c.*y.^2).*((y.^2+h^2).^(-alpha_l/2)/m).^l.*(1+r_th/m*(xx.^2+h^2).^(alpha_nl/2).*(y.^2+h^2).^(-alpha_l/2)).^(-m-l),xx,inf),x)...
                   +2*pi*lambda_bs*prodd_3.*arrayfun(@(xx)integral(@(z)z.*(1-exp(-c.*z.^2)).*((z.^2+h^2).^(-alpha_nl/2)).^l.*(1+r_th*(xx.^2+h^2).^(alpha_nl/2).*(z.^2+h^2).^(-alpha_nl/2)).^(-1-l),xx.^(alpha_l/alpha_nl),inf),x)).^(m_l(l));
               prodd = prodd.*prodd_term;
           end
           sum_2 = sum_2 + prodd;
       end
       sumterm = (-1)^k.*(r_th*(x.^2+h^2).^(alpha_nl/2)).^k.*sum_2;
    end
    sum = sum + sumterm;
end
%r_1 = sumterm
r = sum;