% Mean squared displacement for a particular dynein track. 

function [msd,E,StdE] = MSD_Standalone(x,y,pixelfactor)

m = length(y);

k = 30;
msd = zeros(1,k);
E = zeros(1,k);
StdE = zeros(1,k);
c = cell(1,k);

for i = 1:k;
     j = 1;
     q = 0;
         
     while i+j <=  m
        q = q+1;
        c{i}(q) = ((pixelfactor)^2)*((x(i+j)-x(j))^2 + (y(j+i)-y(j))^2);
        j = j+i;              
     end
    
    E(i) = std(c{i});
    msd(i) = mean(c{i});
    StdE(i) = E(i)/sqrt(q);
        
 end