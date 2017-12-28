%curve fitting
%h = 20; %antenna height
varphi = 11.95;
beta = 0.14;
counter = 1;
r = linspace(0,pi/2,20);
x = length(r);
for x_r = 1 : x
    p(x_r) = 1./(1+varphi*exp(-beta*(180/pi*atan(r(x_r))-varphi)));
    %p_f(counter) = 1.*exp(-0.001969*r^2);
    %x(x_r) = r;
    %counter = counter + 1;
end
plot(rad2deg(r),p,'k-');
hold on;
%plot(x,p_f,'r--');