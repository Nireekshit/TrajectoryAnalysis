% Program to find the intensity of the tracks. 

function[intensity,meanintensity,Name] = Intensity(location) 
% In the location, enter the directory where tracking files are saved. 

files = dir(location);

files(strncmp({files.name}, '.', 1)) = [];

l = length(files);

meanintensity = zeros(1,l);
int = (cell(1,l));

for a = 1:l

% Extracts the name of the file whose serial number is a.    
Name{a} = files(a).name;  

% Loads the file 
A = load(Name{a});

int{a} = A(:,11);

intensity(a) = int{a}(1);

meanintensity(a) = mean(int{a});

end


end

