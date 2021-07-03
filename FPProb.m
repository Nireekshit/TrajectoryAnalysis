% First passage probability

% We intend to plot the first passage probabilities of a population of
% particles. 

% The run-lengths from our experience are exponentially distributed. We
% have to include pauses as well. So each iteration of the run loop, could
% be either stationary or moving. 

% If it is moving - then the runlength should be obtained from an
% exponential distribution. It should be stationary 80 % of the time.


% numpart is number of particle. 
% frames is number of movements observed
% framespacing is the timeinterval between frames in sec
% lambda is the mean run length/time from an exponential distribution 
% bias (0-1) is the fraction of time spent moving to minus end. 
% f (0-100) is percentage of time spent in stationary state 
% k is the vector with lengths at which we want to find FPP. 


function [Data] = FPProb(numpart,frames,lambda,bias,f,k)

for i = 1:numpart % For all the particles
    
        D(1,i) = 0;

        j =  1;
        
       
    while j < frames % This is a sequence of frames

        % Next we generate run times from an exponential distribution
        
        r = rand(1);
        
        if r > f/100 % This is to see that f % of the time is spent in stationary state 
            
                 c = rand(1);

                 p = (1/(lambda)) * log(1/(1-c)) ; % Generating runs from exponentail distribution
                                                
                         
                 vector = [];
                                 
                 p = round(p,1); % This wont work for all frame spacing. Only for 0.1
                 
                 if p >= 0.1
            
                 z = rand(1);
            
                       if z < bias/100 % Deciding the probability of minus end directed runs
                
                           p = -p;
                 
                       end
               
                 runlength = abs(p)/0.1;          % This is the number of 
                                                  % frames for which the
                                                  % particle runs
                                                  % continuously assuming 1
                                                  % um/s.
                  
                                                  
                                                  
                 if (j + runlength) <= frames                                  
            
                    vector = [sign(p)*0.1: sign(p)*0.1 : p];
                    
                 else
                     
                    remainingframes = (frames-j)*0.1;

                           
                    vector = [sign(p)*0.1: sign(p)*0.1 : sign(p)*remainingframes];
                    
                 end
                     
            
                     
                 % now we increment the displacement. Divide the total
                 % displacement obtained from p into 0.1 increments. 
                 % Assuking 1 um/s velocity. Else the code will change. 
            
                       for s = 1:length(vector)
                 
                           D(j+s,i) = D(j,i) + vector(s);
                
                       end
                 
                 
                 j = j + s;
                 
                 
                 else
                     
                     D(j+1,i) = D(j,i);
                     
                     j = j + 1;
                     
                 end
                 
             
        else 
            
                D(j+1,i) = D(j,i);
                
                j = j + 1;
            
        end
        
    end 
            
            
end

[m,n] = size(D); % Get the size of the displacement matrix. To make life easy. 

maximum = k(end); % To use later



for i = 1:n % The first loop runs over all the particles

for w = 1:length(k) % This one is for each particle trajectory. 

     
    q = find((D(:,i)) >= k(w)); % Find all the time points at which the particle's position is greater than k(w);
                                % Only in one direction though.  

                          
    t = isempty(q); % q could be empty if the particle never reached that l.
                    % Moreover, q could be multiple time points. 
                    % But we want to find the first time. 
    
    if t == 0

       FP(w,i) = 0.1*(q(1)); % converting frames to time. Use 0.1 for now. 

    else 

       FP(w,i) = (0); % if q is empty we assign it zero value. Can be removed later. 
       
    end

end 

end % end of loop to find the first time points at which positions defined in k are reached. 

minimum = 0.1; space = 0.5; BinEdges = [minimum:space:maximum]; % For the histogram counts

e = [(minimum + space/2) : space : (maximum-space/2)];

% The next loop calculates actual FPPs.

for i = 1:length(k) % The different length scales at which we want FPP

    Y  = FP(i,:);

    FPtimes{i} = Y(find(Y>0)); % Removing zero values. 

    [DistFPtimes{i},BinEdges] = histcounts(FPtimes{i},BinEdges); % counting the number of elements in each bin 

   FPP{i} = DistFPtimes{i}/length(FPtimes{i}); % Normalization 
     
    isNZ = ( FPP{i} > 0);  FPP{i} = FPP{i}(isNZ); E{i} = e(isNZ); % This is to remove bins with zero values. 
    
end

  
   Data.Numpart                  = numpart;
   Data.Frames                   = frames;
   Data.Lambda                   = lambda;
   Data.Bias                     = bias;
   Data.StatFrac                 = f;
   Data.ParticleDisplacement     = D;   
   Data.FirstPassProb            = FPP;
   Data.Positions                = E;


    
