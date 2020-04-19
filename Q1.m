clc;
clear all;

a=(2+1)/3;
L1=((3+6+3)/5)*a;
alpha=(3+2)*pi/20;
discritArr=[200,500,1000,2000,3000];
error_matrix=zeros(5,2);%help matrix for calculating the error for every N

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
        b_column_vector(p+1)=cos(pi*(-pi+alpha+p*phi_jump)/(pi-alpha)); %updating b_column_vector for potential in z=0
    end
    
    for p=0:phi_divided      %loop over phi axis
        b_column_vector(Z_divided*(phi_divided+1)+p+1)=-cos(pi*(-pi+alpha+p*phi_jump)/(pi-alpha)); %updating b_column_vector for potential in z=L1
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
                 potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p-1)=1/(a*delta_phi)^2; %back dot
                 potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p+1)=1/(a*delta_phi)^2; %forward dot
                 potential_mat((z-1)*(phi_divided+1)+p,(z-2)*(phi_divided+1)+p)=1/delta_z^2; %down dot
                 potential_mat((z-1)*(phi_divided+1)+p,(z)*(phi_divided+1)+p)=1/delta_z^2; %up dot
                 potential_mat((z-1)*(phi_divided+1)+p,(z-1)*(phi_divided+1)+p)=-2/(a*delta_phi)^2-2/delta_z^2; %dot itself    
            end 
        end
    end
    x_potential_vector=potential_mat\b_column_vector; 
    %calculating the error - c part
    for j=1:size(x_potential_vector)
        if(mod(j,(phi_divided+1))==0)  %if phi=pi-alpha
            error_matrix(k,1)=error_matrix(k,1)+abs(x_potential_vector(j)-analytic_potential_func(floor(j/(phi_divided+1))*delta_z,+pi-alpha));
            error_matrix(k,2)=error_matrix(k,2)+abs(analytic_potential_func(floor(j/(phi_divided+1))*delta_z,+pi-alpha));
        else   %if phi!=pi-alpha
            error_matrix(k,1)=error_matrix(k,1)+abs(x_potential_vector(j)-analytic_potential_func(floor(j/(phi_divided+1))*delta_z,(mod(j,(phi_divided+1))-1)*phi_jump-pi+alpha));
            error_matrix(k,2)=error_matrix(k,2)+abs(analytic_potential_func(floor(j/(phi_divided+1))*delta_z,(mod(j,(phi_divided+1))-1)*phi_jump-pi+alpha));    
        end
    end  
end

% b part

% z=L1/4
z_index=round(L1/4/delta_z)+1;    %calculating the z'th index for L/4
numeric_vector=zeros(1,phi_divided+1);
analytic_vector=zeros(1,phi_divided+1);
phi_vector=linspace(-pi+alpha,pi-alpha,phi_divided+1);
for p=1:phi_divided+1       %loop over phi axis
    numeric_vector(p)=x_potential_vector((z_index-1)*(phi_divided+1)+p);  
    analytic_vector(p)=analytic_potential_func(L1/4,phi_vector(p));
end
figure(1);
plot(phi_vector,numeric_vector,'color','r'); hold on;
plot(phi_vector,analytic_vector,'color','b'),title('numeric and analytic solution for Z=L/4'), xlabel('-pi+alpha<phi< pi-alpha'),ylabel('potential');
% z=3L1/4 
z_index=round(3*L1/4/delta_z)+1;    %calculating the z'th index for 3L/4
numeric_vector=zeros(1,phi_divided+1);
analytic_vector=zeros(1,phi_divided+1);
phi_vector=linspace(-pi+alpha,pi-alpha,phi_divided+1);
for p=1:phi_divided+1       %loop over phi axis
    numeric_vector(p)=x_potential_vector((z_index-1)*(phi_divided+1)+p);  
    analytic_vector(p)=analytic_potential_func(3*L1/4,phi_vector(p));
end
figure(2);
plot(phi_vector,numeric_vector,'color','r'); hold on;
plot(phi_vector,analytic_vector,'color','b'),title('numeric and analytic solution for Z=3L/4'), xlabel('-pi+alpha<phi< pi-alpha'),ylabel('potential');
% phi=(pi-alpha)/3
phi_index=round((((pi-alpha)/3)-(-pi+alpha))/phi_jump)+1;    %calculating the phi'th index for (pi-alpha)/3
numeric_vector=zeros(1,Z_divided+1);
analytic_vector=zeros(1,Z_divided+1);
z_vector=linspace(0,L1,Z_divided+1);
for z=1:Z_divided+1       %loop over Z axis
    numeric_vector(z)=x_potential_vector(phi_index+(z-1)*(phi_divided+1));  
    analytic_vector(z)=analytic_potential_func((z-1)*delta_z,(pi-alpha)/3);
end
figure(3);
plot(z_vector,numeric_vector,'color','r'); hold on;
plot(z_vector,analytic_vector,'color','b'), title('numeric and analytic solution for phi=(pi-alpha)/3'), xlabel('0<z<L'),ylabel('potential');


% plot of part c after calculated in the loop at the beginning
error_vector=zeros(1,5);
for k=1:length(discritArr)       %iteration on the number of discritizations 
    error_vector(k)=error_matrix(k,1)/error_matrix(k,2);  %final error for every N
end
figure(4);
plot(discritArr,error_vector), title('E-relative graph(N)'), xlabel('N'),ylabel('E-relative');         % plot (error for N)
disp(error_vector);
