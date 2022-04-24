function [phi,coef] = Gram_schmidt(m,T)

[n,M]=size(m);
           
Tb=T/n;
m2= m.^2;
m2_sum=sum(m2);
energy=sum(m2)*Tb;

phi=zeros(n,M);

phi(:,1)=m(:,1)/sqrt(energy(1));
coef= zeros(M,M);
coef(1,1)=sqrt(energy(1));

for j= 2 : M                 %%%iterates on col(phi loop)
    
    g = m(:,j);
    
    for k = 1:j-1             %%%%%% filling coef above s(j,j)              
        coef(k,j) = (m(:,j)' * phi(:,k))*Tb;
        g = g - coef(k,j)*phi(:,k);
    end
    
    energy_g = sum(g.^2)*Tb;
    if (energy_g == 0 || energy_g < 1e-5) 
        continue;
    end
    phi(:,j) = g/sqrt(energy_g);
   
    
end

for z = 2 : M  %calculate the diagonal coefficients 
    coef(z,z)= (m(:,z)' * phi(:,z))*Tb;
end
phi=phi(:,any(phi));                %remove zeros in phi
coef2=coef';
coef=coef2(:,any(coef2));           %remove zeros in phi
              
end