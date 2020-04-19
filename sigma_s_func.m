function s = sigma_s_func(z,phi) %function of the conductivity as given at question 3
	a=(2+1)/3;
    L1=((3+6+3)/5)*a;
    alpha=(3+2)*pi/20;
    c=-(2+6)/12;
    sigma_0=1;
    d_sigma=-exp(c);
    z_c=L1/2;
    phi_c=0;
    w_phi=(pi-alpha)/(2+2/5);
    w_z=L1/(3+6/5);

   s=sigma_0+0.25*d_sigma*((1+cos(2*pi*(phi-phi_c)./w_phi))*(1+cos((2*pi*(z-z_c)./w_z))));
 
end
