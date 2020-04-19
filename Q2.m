clc;
clear all;

a=(2+1)/3;
L1=((3+6+3)/5)*a;
alpha=(3+2)*pi/20;
discritArr=[200,500,1000,2000,3000];
density_matrix=zeros(5,1);%density for every N
% a part

for k=1:length(discritArr)       %iteration on the number of discritizations  
    N=discritArr(k);    
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
                 potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p-1)=(1/((a*delta_phi)^2))*(1+delta_z*(z-1)/L1); %back dot
                 potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p+1)=(1/((a*delta_phi)^2))*(1+delta_z*(z-1)/L1); %forward dot
                 potential_mat((z-1)*(phi_divided+1)+p,(z-2)*(phi_divided+1)+p)=(1/delta_z^2)*(1+delta_z*(z-1)/L1)-1/(2*L1*delta_z); %down dot
                 potential_mat((z-1)*(phi_divided+1)+p,(z)*(phi_divided+1)+p)=(1/delta_z^2)*(1+delta_z*(z-1)/L1)+1/(2*L1*delta_z); %up dot
                 potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p)=(-2/(a*delta_phi)^2-2/(delta_z)^2)*(1+delta_z*(z-1)/L1); %dot itself  
            end
        end
    end
    x_potential_vector=potential_mat\b_column_vector;
end
% transform vector potential into matrix potential
potential_matrix=zeros(Z_divided+1 ,phi_divided+1);     % potential matrix - for every point
for z=1:Z_divided+1       %loop over z axis
        for p=1:phi_divided+1     %loop over phi axis
            potential_matrix(z,p)=x_potential_vector((z-1)*(phi_divided+1)+p);
        end
end
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
z_potential_dermat=-z_potential_dermat;
phi_potential_dermat=-phi_potential_dermat;
multi_z_potential_dermat=z_potential_dermat.*z_potential_dermat;  %E-Z^2
multi_phi_potential_dermat=phi_potential_dermat.*phi_potential_dermat;    %E-PHI^2
multi_z_potential_dermat=multi_z_potential_dermat+multi_phi_potential_dermat;  %E-Z^2+E-PHI^2
%calculate p_d
p_d_mat=zeros(Z_divided+1 ,phi_divided+1);
for z=1:Z_divided+1       %loop over rows
    for p=1:phi_divided+1       %loop over the columns
        p_d_mat(z,p)=multi_z_potential_dermat(z,p)*(1+(z-1)*delta_z/L1);
    end
end
%integration over p_d
electric_power=0;
for z=1:Z_divided+1       %loop over rows
    for p=1:phi_divided+1       %loop over the columns
        electric_power=electric_power+p_d_mat(z,p);
    end
end
electric_power=electric_power*a*delta_phi*delta_z;
disp(electric_power);
res=1/electric_power;
disp(res);