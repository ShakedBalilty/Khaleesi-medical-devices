function y = analytic_potential_func(z,phi) %function of the potential according to the analytic solution
    a=(2+1)/3;
    L1=((3+6+3)/5)*a;
    alpha=(3+2)*pi/20;
    arg=pi/(a*(pi-alpha));
    y=(-(1+cosh(arg*L1))/sinh(arg*L1)*sinh(z.*arg)+cosh(z.*arg))*cos(pi*phi./(pi-alpha));
end
