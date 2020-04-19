clc;
clear all;

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

l_lim_phi=phi_c-0.5*w_phi;
r_lim_phi=phi_c+0.5*w_phi;
l_lim_z=z_c-0.5*w_z;
r_lim_z=z_c+0.5*w_z;

N=3000;
phi_divided=sqrt(N*(2*pi-2*alpha)*a/L1);       %the devision of phi into elements in order to get a square
Z_divided=N/phi_divided;       %the devision of L1 into elements in order to get a square
phi_divided=round(phi_divided);   
Z_divided=round(Z_divided);
delta_z=L1/Z_divided;        % length of element in Z axis
delta_phi=(2*pi-2*alpha)*a/phi_divided;      % length of element in phi axis
total_num_dots=(phi_divided+1)*(Z_divided+1);    %the total num of dots we need to look for
phi_jump=(2*pi-2*alpha)/phi_divided;
potential_mat=zeros(total_num_dots,total_num_dots);      %potential matrix restarted to zeros
b_column_vector=zeros(total_num_dots,1);         % column vector b in the equation potential_mat*x=b

for p=0:phi_divided      %loop over phi axis
    b_column_vector(p+1)=0; %updating b_column_vector for potential in z=0
end

for p=0:phi_divided      %loop over phi axis
    b_column_vector(Z_divided*(phi_divided+1)+p+1)=1; %updating b_column_vector for potential in z=L1
end

conductivity_mat=ones(Z_divided+1,phi_divided+1);
for z=1:Z_divided+1       %loop over rows
    for p=1:phi_divided+1       %loop over the columns
        if (l_lim_z<((z-1)*delta_z)&&((z-1)*delta_z)<r_lim_z) % Unevenly conductivity
            if(l_lim_phi<((p-1)*phi_jump-pi+alpha)&&((p-1)*phi_jump-pi+alpha)<r_lim_phi) 
                conductivity_mat(z,p)=sigma_s_func((z-1)*delta_z,(p-1)*phi_jump-pi+alpha);
            end
        end
     end
end
        

for z=1:Z_divided+1       %loop over z axis
    for p=1:phi_divided+1       %loop over phi axis
        if (z==1)                   %z=0
             potential_mat(p,p)=1; %dot itself  
        elseif (z==Z_divided+1)         %z=L1    
             potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p)=1; %dot itself  
        elseif (p==1)                 %phi=0
             potential_mat((z-1)*(phi_divided+1)+1,(z-1)*(phi_divided+1)+2)=1; %right dot
             potential_mat((z-1)*(phi_divided+1)+1,(z-1)*(phi_divided+1)+1)=-1; %dot itself
        elseif (p==phi_divided+1)      %phi=2*pi-2*alpha
             potential_mat(z*(phi_divided+1),z*(phi_divided+1))=1; %dot itself
             potential_mat(z*(phi_divided+1),z*(phi_divided+1)-1)=-1; %left dot 
        else                         % dot in the middle
             sigma_s_dot=conductivity_mat(z,p);
             sigma_s_dot_back=conductivity_mat(z,p-1);
             sigma_s_dot_forward=conductivity_mat(z,p+1);
             sigma_s_dot_up=conductivity_mat(z+1,p);
             sigma_s_dot_down=conductivity_mat(z-1,p);
             
             potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p-1)=(sigma_s_dot/((a*delta_phi)^2))-(sigma_s_dot_forward-sigma_s_dot_back)/((2*a*delta_phi)^2); %back dot
             potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p+1)=(sigma_s_dot/((a*delta_phi)^2))+(sigma_s_dot_forward-sigma_s_dot_back)/((2*a*delta_phi)^2); %forward dot
             potential_mat((z-1)*(phi_divided+1)+p,(z-2)*(phi_divided+1)+p)=(sigma_s_dot/delta_z^2)-(sigma_s_dot_up-sigma_s_dot_down)/((2*delta_z)^2); %down dot
             potential_mat((z-1)*(phi_divided+1)+p,(z)*(phi_divided+1)+p)=(sigma_s_dot/delta_z^2)+(sigma_s_dot_up-sigma_s_dot_down)/((2*delta_z)^2); %up dot
             potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p)=((-2/(a*delta_phi)^2-2/(delta_z)^2))*sigma_s_dot; %dot itself  
        end
    end
end
x_potential_vector=potential_mat\b_column_vector;

potential_matrix=zeros(Z_divided+1,phi_divided+1);%converting the vector into matrix
for z=1:Z_divided+1       %loop over z axis
    for p=1:phi_divided+1       %loop over phi axis
        potential_matrix(z,p)=x_potential_vector(p+(z-1)*(phi_divided+1));
    end
end

k_z_vec=zeros(Z_divided+1 ,1);  %z_value at each point
k_phi_vec=zeros(1 ,phi_divided+1);  %phi_value at each point
for z=1:Z_divided+1       %loop over rows
        k_z_vec(z,1)= (z-1)*delta_z;
end
for p=1:phi_divided+1       %loop over the columns
    k_phi_vec(1,p)=(p-1)*phi_jump-pi+alpha;
end
figure(1);
contourf(k_phi_vec,k_z_vec,potential_matrix),title('map of the potential(phi,z)'), xlabel('phi'),ylabel('z'); % potential plot

z_potential_dermat=zeros(Z_divided+1 ,phi_divided+1);      %derivative potential matrix in z axis
phi_potential_dermat=zeros(Z_divided+1 ,phi_divided+1);     %derivative potential matrix in phi axis
%updating  z potential derivative matrix
for z=1:Z_divided+1       %loop over rows
    for p=1:phi_divided+1       %loop over the columns
            if (z==1)                   %z=0
                 z_potential_dermat(z,p)=(potential_matrix(2,p)-potential_matrix(1,p))/delta_z; %forward derivative
            elseif (z==Z_divided+1)         %z=L1  
                 z_potential_dermat(z,p)=(potential_matrix(Z_divided+1,p)-potential_matrix(Z_divided,p))/delta_z; %back derivative
            else                         % dot in the middle
                 z_potential_dermat(z,p)=(potential_matrix(z+1,p)-potential_matrix(z-1,p))/(2*delta_z); %central derivative
            end
    end
end

%updating  phi potential derivative matrix
for z=1:Z_divided+1       %loop over rows
    for p=1:phi_divided+1       %loop over the columns
            if (p==1)                   %p=-pi+alpha
                 phi_potential_dermat(z,p)=(potential_matrix(z,2)-potential_matrix(z,1))/(a*delta_phi); %forward derivative
            elseif (p==phi_divided+1)         %p=pi-alpha  
                 phi_potential_dermat(z,p)=(potential_matrix(z,phi_divided+1)-potential_matrix(z,phi_divided))/(a*delta_z); %back derivative
            else                         % dot in the middle
                 phi_potential_dermat(z,p)=(potential_matrix(z,phi_divided+1)-potential_matrix(z,phi_divided-1))/(2*a*delta_z); %central derivative
            end
    end
end

z_points_mat=zeros(Z_divided+1 ,phi_divided+1);    %z_value at each point
phi_points_mat=zeros(Z_divided+1 ,phi_divided+1);  %phi_value at each point
k_z_mat=zeros(Z_divided+1 ,phi_divided+1);  %z_value at each point
k_phi_mat=zeros(Z_divided+1 ,phi_divided+1);  %phi_value at each point
for z=1:Z_divided+1       %loop over rows
    for p=1:phi_divided+1       %loop over the columns
        z_points_mat(z,p)=(z-1)*delta_z;
        phi_points_mat(z,p)=(p-1)*phi_jump-pi+alpha; 
    end
end

z_potential_dermat=-z_potential_dermat; % E at the z direction
phi_potential_dermat=-phi_potential_dermat;% E at the phi direction
figure(2);
quiver(phi_points_mat, z_points_mat,phi_potential_dermat, z_potential_dermat),title('map of the electric field(phi,z)'), xlabel('phi'),ylabel('z');%,title('map of the electric field(phi,z)'), xlabel('phi'),ylabel('z'); % Vector field

k_z_mat= conductivity_mat.*z_potential_dermat;
k_phi_mat= conductivity_mat.*phi_potential_dermat;
figure(3);
quiver(phi_points_mat, z_points_mat,k_phi_mat,k_z_mat),title('map of the surface current(phi,z)'), xlabel('phi'),ylabel('z'); % surface current map
p_d_z_matrix= z_potential_dermat.*k_z_mat;% p_d_ at the z direction
p_d_phi_matrix=phi_potential_dermat.*k_phi_mat;% p_d_ at the phi direction
p_d_matrix=p_d_z_matrix+p_d_phi_matrix;
figure(4);
contourf(k_phi_vec,k_z_vec,p_d_matrix),title('map of p_d(phi,z)'), xlabel('phi'),ylabel('z'); % p_d plot

%integration over p_d
electric_power=0;
for z=1:Z_divided+1       %loop over rows
    for p=1:phi_divided+1       %loop over the columns
        electric_power=electric_power+p_d_matrix(z,p);
    end
end
electric_power=electric_power*a*delta_phi*delta_z;
disp(electric_power);



