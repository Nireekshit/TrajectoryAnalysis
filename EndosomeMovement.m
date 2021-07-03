% This function can be used to understand the dynamics of endosomes
% moving inside a cell. 

function [Endosome] = EndosomeMovement(location,fps,pixelfactor,TH)

% location - consists of multiple folders each corresponding to a cell
%            in which many endosomes are tracked. The tracking is done using the 
%            Low Light Tracking Tool (Krull et al., 2014) and each of the 'trcT'
%            files contains x,y information with respect to every frame. Other data 
%            like intensity, flux of photons is part of the trcT file.
%            Enter with a '/' in the end. Helps with strcat later. 
        
% fps      - is the image acquisition rate. For example 10 frames per
%            second. 

% pixelfactor - is the distance in micrometers corresponding to each pixel
%               movement. For example, 0.16 um/pixel. [For 100x objective
%               and 16 um pixel size]

% TH          - is the threshold displacements below which displacements can be 
%               ignored. For example, 50nm = 0.05. 

d  = dir(location); % Extracts files in the location


d(strncmp({d.name}, '.', 1)) = []; % removes some extra files from the imported directory d
d(strncmp({d.name}, '.tif', 1)) = []; % same as above. 

% Next we store all the image names (subfolders) seperately. Each of the images along
% with trcT files is stored in a folder. 

imagenames = cell(length(d),1);

for i = 1:length(d)
    imagenames{i} = d(i).name;
end


% Each of the files contains multiple endosome tracks, the tif
% image and also the location of a center. We save all the trcT files to a 
% directory and save the center seperately as centrosome. 

for q = 1:length(imagenames)
    
    files = struct();
    files = dir(strcat(location,imagenames{q})); % importing files present within
                                                 % each folder.  

    files(strncmp({files.name}, '.', 1)) = [];
    files(strncmp({files.name},'errors',1)) = []; % removing a file with the name 'errors.txt'
 
           String = 'centre'; % Will be used to find a file with the name center. 

           Name  = cell(length(files),1);
           found = cell(length(files),1);
            
           for i = 1:length(files)
                Name{i} = files(i).name;
                found{i} = strfind(Name{i},String);
           end


           % Centrosome is that file in which found == 1, from the previous loop
           
           Centrosome =  find(~cellfun(@isempty,found)); % Finding the index of the file corresponding
                                                         % to the center. 

           Centrosomedetails = load(files(Centrosome).name);

           Centrosomeposition = [Centrosomedetails(3),Centrosomedetails(4)];
           
           Endosomeinfo = struct();

           files(Centrosome) = []; % remove the centrosome file from the list of files. 

           
          % These are the files contained within one folder(i.e representing one cell) 
          
          te = 0; 
          
          for i = 1:length(files) 
              
              [path,filename,ext] = fileparts(files(i).name);     
              
              % Thw switch will add only those files with trcT extension
              % into the 'Endosome' variable. 
                  
              switch ext
                  
                  case '.trcT'
                     
                   te = te + 1;
                     
                   Endosomeinfo(te).name        = files(i).name;
                   Endosomeinfo(te).data        = load(Endosomeinfo(te).name);
                   [m,~]                        = size(Endosomeinfo(te).data);
                   Endosomeinfo(te).length      = m;
              end
             
          end
          
          % The indices of the tracks b/w files and Endosomeinfo have
          % changed. Have to be careful when writing the rest of the code. 
          
         
        % Next we calculate various parameters to understand the dynamics of the endosome. 
        
counter = 0;
        
for i = 1:length(Endosomeinfo)   
    
        counter = counter + 1;
    
        Endosome(counter).name         = Endosomeinfo(i).name;
        Endosome(counter).intensity    = Endosomeinfo(i).data(1,11);
        Endosome(counter).allintensity = Endosomeinfo(i).data(:,11);
        Endosome(counter).data         = Endosomeinfo(i).data;
        Endosome(counter).length       = Endosomeinfo(i).length;
        Endosome(counter).centrosome   = Centrosomeposition;
        
        % The next loop calculates the distance of the endosome from the
        % cell center at various time points
        
        for j = 1:Endosome(counter).length
            Endosome(counter).distance(j) = pixelfactor*(sqrt((Endosome(counter).data(j,3)-Centrosomeposition(1))^2 + (Endosome(counter).data(j,4)-Centrosomeposition(2))^2));
        end
        
       % In the next loop we calculate the displacement in each frame and
       % the direction of movement. If the displacement is less than 30nm
       % it is within the average error with the given intensity of
       % particles and we assign it a zero value. 
              
            
        for j = 1:Endosome(counter).length-1
            
            % In the field direction index j indicates the direction
            % between j and j+1. So naturally, there is no direction in the 
            % index (length-1).
            
            Endosome(counter).Direction(j) = sign(Endosome(counter).distance(j+1)-Endosome(counter).distance(j));
            
            % Movement is displacement between each frame and therefore
            % there is no movement in the last frame. 
            
            Endosome(counter).Movement(j) = (Endosome(counter).Direction(j))* pixelfactor * (sqrt((Endosome(counter).data(j+1,3)-Endosome(counter).data(j,3))^2 + (Endosome(counter).data(j+1,4)-Endosome(counter).data(j,4))^2));
            
            if (Endosome(counter).Movement(j) > -TH) && (Endosome(counter).Movement(j) < TH)
                
               Endosome(counter).Movement(j) = 0; 
                 
            end 
                      
      
          
            % Displacement is in comparison to the first frame
            
            Endosome(counter).NetDisplacement(j) =  pixelfactor * (sqrt((Endosome(counter).data(j+1,3)-Endosome(counter).data(1,3))^2 + (Endosome(counter).data(j+1,4)-Endosome(counter).data(1,4))^2));
        end
        
          
        % Cumulative movement is the summation of the displacements in
        % every frame. 
        
        Endosome(counter).Cumulative(1) = 0;
        
        
        for j = 1:length(Endosome(counter).Movement)
            Endosome(counter).Cumulative(j+1) = Endosome(counter).Cumulative(j) + Endosome(counter).Movement(j);
        end
        
        Endosome(counter).CumulativeVelocity = Endosome(counter).Cumulative(length(Endosome(counter).Cumulative))*fps/(Endosome(counter).length-1);
        
        
        % We have to find runs in the track of an endosome. Since we
        % reduced movements less than 50nm to zero we do not need to add
        % any thresholds. 
        
        Endosome(counter).Run       = zeros(1,length(Endosome(counter).Movement));
        Endosome(counter).Runlength = zeros(1,length(Endosome(counter).Movement));
        Endosome(counter).Velocity  = zeros(1,length(Endosome(counter).Movement));

        j = 1;
        
while j <= length(Endosome(counter).Direction)
    
        
      if sign(Endosome(counter).Movement(j)) == -1 
          
         Endosome(counter).Run(j) = -1; 
         
         c = 1; k = j+1;
         
         while c == 1 && k <= length(Endosome(counter).Movement)
               c =  (sign(Endosome(counter).Movement(k)) == -1);
               k = k + c;
               Endosome(counter).Run(j) = Endosome(counter).Run(j) - c;
         end
        
         Endosome(counter).Runlength(j) = sum(Endosome(counter).Movement(j:k-1));
         
         Endosome(counter).Velocity(j)  = -(Endosome(counter).Runlength(j)/(Endosome(counter).Run(j)*0.1));
         
         j = k;
         
      elseif sign(Endosome(counter).Movement(j)) == 1

         Endosome(counter).Run(j) = 1; 
         
         c = 1; k = j+1;

         while c == 1 && k <= length(Endosome(counter).Movement) 
               c =  (sign(Endosome(counter).Movement(k)) == 1);
               k = k + c;
               Endosome(counter).Run(j) = Endosome(counter).Run(j) + c;
         end
         
         Endosome(counter).Runlength(j) = sum(Endosome(counter).Movement(j:k-1));
       
         Endosome(counter).Velocity(j)  = (Endosome(counter).Runlength(j))/(Endosome(counter).Run(j)*0.1);
         
         j = k;
         
      else 
          
         j = j + 1;

      end
      
end
      Endosome(counter).Minusrun            = Endosome(counter).Runlength(find(Endosome(counter).Runlength < 0 ));
      Endosome(counter).Minustime           = Endosome(counter).Run(find(Endosome(counter).Run < 0 )); 
      Endosome(counter).Minusvelocity       = Endosome(counter).Velocity(find(Endosome(counter).Velocity < 0));
      Endosome(counter).Totalminus          = sum(Endosome(counter).Run(Endosome(counter).Run < 0));
      Endosome(counter).TotalMinusmovement  = sum(Endosome(counter).Movement(Endosome(counter).Movement < 0));
      Endosome(counter).Minusfraction       = Endosome(counter).Totalminus/(Endosome(counter).length-1);
   
      
      Endosome(counter).Plusrun             = Endosome(counter).Runlength(find(Endosome(counter).Runlength > 0 ));
      Endosome(counter).Plustime            = Endosome(counter).Run(find(Endosome(counter).Run > 0 )); 
      Endosome(counter).Plusvelocity        = Endosome(counter).Velocity(find(Endosome(counter).Velocity > 0));
      Endosome(counter).Totalplus           = sum(Endosome(counter).Run(Endosome(counter).Run > 0));
      Endosome(counter).TotalPlusmovement   = sum(Endosome(counter).Movement(Endosome(counter).Movement > 0));
      Endosome(counter).Plusfraction        = Endosome(counter).Totalplus/(Endosome(counter).length-1);

      Endosome(counter).Totalmovement       = abs(Endosome(counter).TotalMinusmovement) + Endosome(counter).TotalPlusmovement;
      Endosome(counter).Moving              = abs(sum(Endosome(counter).Run(Endosome(counter).Run < 0))) + sum(Endosome(counter).Run(Endosome(counter).Run > 0));
      Endosome(counter).Movingfraction      = (Endosome(counter).Moving)/(Endosome(counter).length-1);
        
        
        
        
end

end


    for i = 1:length(Endosome)
        Endosome(i).Runloc  = find(Endosome(i).Run );
        Endosome(i).NonZeroRun = Endosome(i).Run(find(Endosome(i).Run));
    end 
    
    
    for i = 1:length(Endosome)
        Endosome(i).Runspacing = diff(Endosome(i).Runloc)*0.1;
        
        % The diff function gives the diffference between indices of 2 succesive
        % run starts. But if we want the pause time then we have to
        % subtract the length of the run from the space between start of 2
        % runs. 
        
        for j = 1:length(Endosome(i).NonZeroRun)-1
            Endosome(i).Runspacing(j) = Endosome(i).Runspacing(j) - 0.1*abs(Endosome(i).NonZeroRun(j));
        end
            
        
    end 
    
    for i = 1:length(Endosome)
        for j = 1:length(Endosome(i).Runloc)-1
        Endosome(i).Runnext(j) = Endosome(i).Run(Endosome(i).Runloc(j+1));
        Endosome(i).Afterrun(j) = Endosome(i).Run(Endosome(i).Runloc(j)+1);
        end
    end   

       
    for i = 1:length(Endosome)
        for j = 1:length(Endosome(i).NonZeroRun)-1
        Endosome(i).SignRunnext(j) = sign(Endosome(i).NonZeroRun(j+1)) * sign(Endosome(i).NonZeroRun(j));
        Endosome(i).SignAdjRunnext(j) = Endosome(i).SignRunnext(j) * abs(Endosome(i).Runnext(j));
        end
    end
    
    
    % We can seperate the long runs. Nmin = Enter the minimum number of
    % displacements that can be counted as a run 
    
        for i = 1:length(Endosome)
            Endosome(i).LongRunLoc = find( abs(Endosome(i).Run) >= 4);
        end

        for i = 1:length(Endosome)
            Endosome(i).LongRuns      = Endosome(i).Run(Endosome(i).LongRunLoc);
            Endosome(i).LongRunLength = Endosome(i).Runlength (Endosome(i).LongRunLoc);
        end

        for i = 1:length(Endosome)
            Endosome(i).LongRunSpacing = diff(Endosome(i).LongRunLoc)*0.1;
            
            
        % The diff function gives the diffference between indices of 2 succesive
        % run starts. But if we want the pause time then we have to
        % subtract the length of the run from the space between start of 2
        % runs. 
        
        for j = 1:length(Endosome(i).LongRuns)-1
            Endosome(i).LongRunSpacing(j) = Endosome(i).LongRunSpacing(j) - 0.1*abs(Endosome(i).LongRuns(j));
        end
            
        end
              
        for i = 1:length(Endosome)
            for j = 1:length(Endosome(i).LongRuns)-1
                Endosome(i).SignLongRunNext(j) = sign(Endosome(i).LongRuns(j+1)) * sign(Endosome(i).LongRuns(j));
                Endosome(i).SignAdjLongRunNext(j) = Endosome(i).SignLongRunNext(j) * abs(Endosome(i).LongRuns(j+1));
            end
        end
        
        
        

end



            
            


