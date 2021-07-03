% This function will calculate the ensemble average mean squared 
% displacement of the given set of data. 

% The data should be in the form of a struct with the field data having the
% x-y location in columns. Mention column number of x and y in input. 

% Also input pixel factor. 

% Don't forget to remove the centers from it. 


function [EnsembleMovement] = Ensemble(ParticleInfo,x,y,pixelfactor,k)
               
         d = length(ParticleInfo);
         
          
         for i = 1:d
             [p,~] = size(ParticleInfo(i).data);
             Tracklength(i) = p;
         end
         
         Maxspace = k;
         C = cell(1,length(Tracklength));
         
         for i = 1:d
             C{i} = zeros(max(Tracklength)-1,Maxspace);
         end

for a = 1:d
    
    [m,n] = size(ParticleInfo(a).data);

for  i = 1:k;
     j = 1;
     q = 0;
     
     while i+j <=  m
        q = q+1;
        C{a}(j,i) = (pixelfactor^2)*(((ParticleInfo(a).data(i+j,x)-ParticleInfo(a).data(j,x))^2 + (ParticleInfo(a).data(j+i,y)-ParticleInfo(a).data(j,y))^2));
        j = j+i;              
     end
end

% The loop calculates the squared displacements for different time spaces.

end

i = 1;

while i <= k
    
     MSD_All{i} = 0;

     for j = 1:d
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













    

