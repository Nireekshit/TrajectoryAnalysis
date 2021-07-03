% This function can calculate the properties of a particle moving inside a
% cell.

% Before running the function add the following to the MATLAB directory :
% 1. Folder contaning the m files required for tracking 
% 2. the path to the folder which stores individual images in which
% particles are tracked (named location).

% location    =  the path to the folder which stores individual images in which
%                particles are tracked. 
% fps         =  frames per second of data acquisition
% pixelfactor =  pixel size of the camera/objective magnification in um.

function [Dynein,Tracks,EnsembleMovement] = DyneinMovementTH(location,fps,pixelfactor);

% Extract the directory from the location. 

d  = dir(location);

% The directory contains some files named . , .. 
% These are created by MATLAB when accesing the
% folder. These have to be deleted from the directory. 

d(strncmp({d.name}, '.', 1)) = [];


counter = 0;

% Next store all the image names seperately. 

for i = 1:length(d)
    imagenames{i} = d(i).name;
end

% Images are stored in folder with the exact name as the image. 
% Since each of the folders contains multiple dynein tracks and also the
% location of a center, we save all the trcT files to a directory and save
% the center seperately as centrosome. 

for q = 1:length(imagenames)
    
files = struct();
files = dir(strcat(location,imagenames{q}));

files(strncmp({files.name}, '.', 1)) = [];

String = 'Center';

Name = {}; found = {};

           for i = 1:length(files)
                Name{i} = files(i).name;
                found{i} = strfind(Name{i},String);
           end


           Centrosome =  find(~cellfun(@isempty,found)); % Centrosome is that file in which found == 1, from the previous loop. 

           Centrosomedetails = load(files(Centrosome).name);

           Centrosomeposition = [Centrosomedetails(3),Centrosomedetails(4)];
           
           files(Centrosome) = [];
           
           
% We removed the centrosome from the list of tracks. 
           
           

           % The struct dynein info contains all the information of a
           % particular track
           

           Dyneininfo = struct();
           
           
           t = 0;

          % The image folder has files with different extensions. Select
          % the ones with .trcT extension. 
                   
          for i = 1:length(files) 
              
              
                  [path,filename,ext] = fileparts(files(i).name);
                  
                  switch ext
                  
                  case '.trcT'
                      
                      t = t + 1;
              
                  
                  Dyneininfo(t).name      = files(i).name;
                  Dyneininfo(t).data      = load(Dyneininfo(t).name);
                  [m,~]                   = size(Dyneininfo(t).data);
                  Dyneininfo(t).length    = m;
                  
                  end
                     
          end 
          
          Tracks(q) = length(Dyneininfo);
          
          disp(q)
     



for i = 1:length(Dyneininfo)   
    
        counter = counter + 1;
    
        Dynein(counter).name         = Dyneininfo(i).name;
        Dynein(counter).allintensity = Dyneininfo(i).data(:,11);
        Dynein(counter).intensity    = Dyneininfo(i).data(1,11);
        Dynein(counter).data         = Dyneininfo(i).data;
        Dynein(counter).length       = Dyneininfo(i).length;
        
        for j = 1:Dynein(counter).length % The data has the same length as the dynein track(l)
            Dynein(counter).centerdistance(j) = pixelfactor*(sqrt((Dynein(counter).data(j,3)-Centrosomeposition(1))^2 + (Dynein(counter).data(j,4)-Centrosomeposition(2))^2));
        end
        
        for j = 1:Dynein(counter).length-1 % The data has length l-1. The first data point corresponds to displacement 1. 
            Dynein(counter).Direction(j) = sign(Dynein(counter).centerdistance(j+1)-Dynein(counter).centerdistance(j));
        end
        
       
        for j = 1:Dynein(counter).length-1 % The data has length l-1. The first data point corresponds to displacement 1.   
            Dynein(counter).Movement(j) = Dynein(counter).Direction(j)* pixelfactor * (sqrt((Dynein(counter).data(j+1,3)-Dynein(counter).data(j,3))^2 + (Dynein(counter).data(j+1,4)-Dynein(counter).data(j,4))^2));
            if abs(Dynein(counter).Movement(j)) < 0.08
               Dynein(counter).Movement(j) = 0;
            end
        end
        
        for j = 1:Dynein(counter).length-1 % The data has length l-1. The first data point corresponds to frame 1.   
            Dynein(counter).Displacement(j) = pixelfactor * (sqrt((Dynein(counter).data(j+1,3)-Dynein(counter).data(1,3))^2 + (Dynein(counter).data(j+1,4)-Dynein(counter).data(1,4))^2));
        end
        
        Dynein(counter).NetVelocity = sign(Dynein(counter).centerdistance(end)-Dynein(counter).centerdistance(1))*Dynein(counter).Displacement(end)*fps/(Dynein(counter).length);
           
        Dynein(counter).Cumulative(1) = 0;
        
        for j = 1:length(Dynein(counter).Movement) % The data in field cumulative has length l. 
             Dynein(counter).Cumulative(j+1) = Dynein(counter).Cumulative(j) + Dynein(counter).Movement(j);
        end
        
        Dynein(counter).CumulativeVelocity = Dynein(counter).Cumulative(length(Dynein(counter).Cumulative))*fps/(Dynein(counter).length);
        
                    
           
end




end
        [EnsembleMovement] = Ensemble(Dynein);
               
        % Parsing cargo trajectories on the basis of net movement
        
        
        
        
     
end


      

