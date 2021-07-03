% One dimensional MSD. 

% The input data x should be a matrix with the distance data of many
% particles. The columns should represent total number of particles. All 
% particle tracks should be of same length, the total number of rows. 

function [EnsembleMovement] = MSD_OneDimension_Ensemble(x)

[m,n] = size(x);

Maxspace = floor(m/2);

         C = cell(1,m);
         
         for i = 1:n
             C{i} = zeros(m-1,Maxspace);
         end

for a = 1:n   
    
    A = x(:,a);

for  i = 1:Maxspace

     j = 1;

     q = 0;
     
     while i+j <=  m
        q = q+1;

        C{a}(j,i) = (A(i+j)-A(j))^2 ;

        j = j+i;              
     end
end

% The loop calculates the squared displacements for different time spaces.

end

i = 1;

while i <= Maxspace
    
     MSD_All{i} = 0;

     for j = 1:n
         MSD_All{i} = [MSD_All{i},(C{j}(:,i))'];
     end
     
     i = i+1;
     
end

for i = 1:length(MSD_All)
    
    index = 0;
    
    index = find(MSD_All{i} ~= 0);
    
    MSD_All{i} = MSD_All{i}(index);
    
    EnsembleMovement(i).e = std(MSD_All{i});
    
    EnsembleMovement(i).MSD = mean(MSD_All{i});
    
    EnsembleMovement(i).stde = EnsembleMovement(i).e/sqrt(length(MSD_All{i}));
end

