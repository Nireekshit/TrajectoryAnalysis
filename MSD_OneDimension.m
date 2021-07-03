% One dimensional MSD. 

% The input data x should be the distance vector. 

function [msd,E,StdE] = MSD_OneDimension(x)

m = length(x);

k = floor(length(x)/2);

msd = zeros(1,k);
E = zeros(1,k);
StdE = zeros(1,k);
c = cell(1,k);

for i = 1:k;
     j = 1;
     q = 0;
         
     while i+j <=  m
        q = q+1;
        c{i}(q) = ((x(i+j)-x(j))^2);
        j = j+i;              
     end
    
    E(i) = std(c{i});
    msd(i) = mean(c{i});
    StdE(i) = E(i)/sqrt(q);
        
 end