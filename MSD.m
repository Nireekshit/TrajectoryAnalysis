% function to calculate mean squared displacement

% Trackinfo should be a structure file which has the data about the tracks
% in the field data. 

function [msd,E,StdE] = MSD(Trackinfo)

n = length(Trackinfo);

for s = 1:n

A{s} = 0.16*Trackinfo(s).data;

k=0;
c={};

[m(s),~] = size(A{s});

k = 10;
msd{s} = zeros(1,k);
E{s} = zeros(1,k);
StdE{s} = zeros(1,k);
c = cell(1,k);

 for i = 1:k;
     j = 1;
     q = 0;
         
     while i+j <=  m(s)
        q = q+1;
        c{i}(q) = ((A{s}(i+j,3)-A{s}(j,3))^2 + (A{s}(j+i,4)-A{s}(j,4))^2);
        j = j+i;              
     end
    
    E{s}(i) = std(c{i});
    msd{s}(i) = mean(c{i});
    StdE{s}(i) = E{s}(i)/sqrt(q);
    
end
 
end

end
