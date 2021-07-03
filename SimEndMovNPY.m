% A function file to analyze the simulated endosome trajectories stored in
% files with extension 'npy'.

function [Trajectory,LowResTrajectory] = SimEndMovNPY(LatticeSeparation,FileName,Skip); % Write the FileName with extension. 
 
data = readNPY(FileName);

data = double(data);

[m,n] = size(data);

for i = 1:m    
    Trajectory(i).Location = data(i,:);
    LowResTrajectory(i).Location = Trajectory(i).Location(1:Skip:end);
end

for i = 1:length(Trajectory)   
    
        Trajectory(i).length         = length(Trajectory(i).Location);
        LowResTrajectory(i).length   = length(LowResTrajectory(i).Location);
        
        Trajectory(i).Movement       = LatticeSeparation*diff(Trajectory(i).Location);
        LowResTrajectory(i).Movement = LatticeSeparation*diff(LowResTrajectory(i).Location);
        
end
        
               
      
        
 for i = 1:length(Trajectory)
     
        for j = 1:length(Trajectory(i).Movement)
            Trajectory(i).Cumulative = cumsum(Trajectory(i).Movement);
        end
        
        for j = 1:length(LowResTrajectory(i).Movement)
            LowResTrajectory(i).Cumulative = cumsum(LowResTrajectory(i).Movement);
        end
        
        Trajectory(i).Direction = sign(Trajectory(i).Movement);
        LowResTrajectory(i).Direction = sign(LowResTrajectory(i).Movement);
        
 end
                
        

% Now call a function to calculate runlengths etc. 

[Trajectory] = TrajectoryRunAnalysis(Trajectory);
[LowResTrajectory] = TrajectoryRunAnalysis(LowResTrajectory);

% We also plot the cumulative movement in LowResTrajectories

for i = 1:length(LowResTrajectory)
Cum{i} = LowResTrajectory(i).Cumulative;
end

for i = 1:100
plot([0:1:length(Cum{i})-1]*0.1,Cum{i},'Color',[0.5,0.5,0.5])
hold on
end

set(gcf,'Position',[500 500 500 500])
axis square
set(gca,'FontSize',14)
ylim([-10 10])
xlabel('Time - s')
ylabel('Displacement - \mum')




