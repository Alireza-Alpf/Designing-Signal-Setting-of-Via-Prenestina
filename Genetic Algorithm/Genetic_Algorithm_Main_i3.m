%% Starting point, clear everything in matlab
tic; %Start the timer
clear all;
close all;
clc;

%% Problem Formulation

FitnessFunction=@(C,g,Q,S,T) TDi(C,g,Q,S,T); % FitnessFunction
T=0.25;                                      % [hour] %Time Interval of analysis
L=15.7;                                      % Lost Time
nphase=3;                                    % Number of phases
nIntersections=1;                            % Number of Intersections (static as 1 intersection)

VarSize=[1 nIntersections*nphase+1];         % Decision Chromosome genes based on number of Intersections


nGL=4;                                       % Number of group Lanes
Q=[915,311,1671,369];                        %[veh/h] %Flows for lane groups EB,WB,NB,SB
S=[3098,1417,3114,1807];                     %[veh/h] %Saturation flow corresponding lane groups
y1=[0.295];                                  %[adim] %Saturation degree of the critical movement for phase 1
y2=[0.242];                                  %[adim] %Saturation degree of the critical movement for phase 2
y4=[0.204];                                  %[adim] %Saturation degree of the critical movement for phase 3
Y=(y1+y2+y4);                                %[adim] %Total saturation degree
Cm=L/(1-Y);                                  %[s] %Minimum Cycle length
Cyclemin=Cm;                                 % Lower bound of CYCLE
Cyclemax=5*Cm;                               % Upper bound of CYCLE
%% Genetic Algorithm Parameters

Maxit=100;                                    % Max iteration
nPop=50;                                      % Population Size

pc=0.8;                                        % Crossover Percentage
nc=2*round(pc*nPop/2);                         % Number of Offsprings (parents)
 
pm=0.2;                                       % Mutation Percentage
nm=round(pm*nPop);                             % Number of Mutants
mu=0.1;                                        % Mutation Rate

pe=0.05;                                       % Elite percentage
ne=round(pe*nPop);                             % Number of Elites

beta=8;                                        % Selection Pressure
   
                                         % number of functions evaluated
%% Initialization

% Individual Structure
empty_individual.Green=[];
empty_individual.Cycle=[];
empty_individual.TotalDelay=[];

% Population Structure
pop=repmat(empty_individual,nPop,1);

% Initialize Population
i=1;
y4=0.212; %Calibrated saturation degree for phase 3

NFE=0;   % number of functions evaluated

while i<=nPop  
    
    % Initialize Individual
    pop(i).Cycle=Cyclemin+(Cyclemax-Cyclemin)*rand; % Generate random Cycle
    greenMin1=y1*pop(i).Cycle;                       % lower bound of GREEN LIGHT 1
    greenMin2=y2*pop(i).Cycle;                       % lower bound of GREEN LIGHT 2
    greenMin4=y4*pop(i).Cycle;                       % lower bound of GREEN LIGHT 3
    greenMax1=pop(i).Cycle-L-greenMin2-greenMin4;         % Upper bound of GREEN LIGHT 1
    greenMax2=pop(i).Cycle-L-greenMin1-greenMin4;         % Upper bound of GREEN LIGHT 2
    pop(i).Green(1)=greenMin1+(greenMax1-greenMin1)*rand;   %Generate radom Green light 1
    pop(i).Green(2)=greenMin2+(greenMax2-greenMin2)*rand;   %Generate radom Green light 2
    pop(i).Green(4)=pop(i).Cycle-L-pop(i).Green(1)-pop(i).Green(2);    %Generate radom Green light 3  
    
    if( pop(i).Green(4) < greenMin4 )
          continue;
    end
   
    pop(i).Green(3)= pop(i).Green(1)+ pop(i).Green(2);  % To compute delay we need g1 + g2 value because 
    %of one throw movement that exist during phase 1 and 2
    % Individual Evaluation from Fitness Function
    for j=1:nGL
        % Measure Delay for each traffic light with current congestion
        pop(i).TotalDelay(j)=FitnessFunction(pop(i).Cycle,pop(i).Green(j),Q(j),S(j),T);
        
        NFE=NFE+1;
    end   
    % Summation of Total Delays quotients
    pop(i).TotalDelay= sum(Q.*pop(i).TotalDelay)/sum(Q);
    i=i+1;
end

% Sort Population
TotalDelay=[pop.TotalDelay];
[TotalDelay, SortOrder]=sort(TotalDelay);
pop=pop(SortOrder);

% Store Best Solution
BestSol=pop(1);

% Store Best Fitness
BestDelay=pop(1).TotalDelay;

% Worst Fitness
WorstDelay=pop(end).TotalDelay;

disp(['FIRST Population..........Best TotalDelay = ' num2str(BestDelay)]);
    fprintf('\n')
    disp('Green Timings in seconds:');
    disp(['  Phase 1 Green Time = ' num2str(BestSol.Green(1))]);
    fprintf('\n')
    disp(['  Phase 2 Green Time = ' num2str(BestSol.Green(2))]);
    fprintf('\n')
    disp(['  Phase 1/2 Green Time = ' num2str(BestSol.Green(3))]);
    fprintf('\n')
    disp(['  Phase 3 Green Time = ' num2str(BestSol.Green(4))]);
    fprintf('\n')
    NFE

%% Loop For Number of Iterations
count=0;
for it=1:Maxit

    % Calculate Selection Probabilities
%     P=(TotalDelay/sum(TotalDelay));
    P=exp(-beta*TotalDelay/WorstDelay);
    P=P/sum(P);

    %% Crossover
    popc=repmat(empty_individual,nc/2,2);
    k=1;
    while k<=nc/2
        
        % Select Parents Indices from roulette wheel
        i1=RouletteWheelSelection(P);
        i2=RouletteWheelSelection(P);
              
        % Select Parents
        p1=pop(i1);
        p2=pop(i2);
 
        popc(k,1).Green=p1.Green;
        popc(k,2).Green=p2.Green;
        popc(k,1).Cycle=p1.Cycle;
        popc(k,2).Cycle=p2.Cycle;
        popc(k,1).TotalDelay=p1.TotalDelay;
        popc(k,2).TotalDelay=p2.TotalDelay;
        % Select random crossover point
        i=randi([1 4]);
        
        % crossover randomness
         if(i==1)
              
                popc1=popc(k,1).Green(2);
                popc(k,1).Green(2)= popc(k,2).Green(2);
                popc(k,2).Green(2)= popc1;

                popc1=popc(k,1).Green(1);
                popc(k,1).Green(1)= popc(k,2).Green(1);
                popc(k,2).Green(1)=popc1;
          
                popc(k,1).Green(3)= popc(k,1).Green(1)+popc(k,1).Green(2);
                popc(k,2).Green(3)= popc(k,2).Green(1)+popc(k,2).Green(2);  
                
                popc(k,1).Green(4)= popc(k,1).Cycle-L-popc(k,1).Green(1)-popc(k,1).Green(2);
                popc(k,2).Green(4)= popc(k,2).Cycle-L-popc(k,2).Green(1)-popc(k,2).Green(2);
                
                
        elseif(i==2)
    
                popc1=popc(k,1).Green(2);
                popc(k,1).Green(2)= popc(k,2).Green(2);
                popc(k,2).Green(2)= popc1;

                popc1=popc(k,1).Green(4);
                popc(k,1).Green(4)= popc(k,2).Green(4);
                popc(k,2).Green(4)=popc1;
                
                popc(k,1).Green(1)= popc(k,1).Cycle-L-popc(k,1).Green(4)-popc(k,1).Green(2);
                popc(k,2).Green(1)= popc(k,2).Cycle-L-popc(k,2).Green(4)-popc(k,2).Green(2);
        
                popc(k,1).Green(3)= popc(k,1).Green(1)+popc(k,1).Green(2);
                popc(k,2).Green(3)= popc(k,2).Green(1)+popc(k,2).Green(2);  
                
                
         elseif(i==3)
             
                popc1=popc(k,1).Green(1);
                popc(k,1).Green(1)= popc(k,2).Green(1);
                popc(k,2).Green(1)= popc1;

                popc1=popc(k,1).Green(4);
                popc(k,1).Green(4)= popc(k,2).Green(4);
                popc(k,2).Green(4)=popc1;
                
                popc(k,1).Green(2)= popc(k,1).Cycle-L-popc(k,1).Green(4)-popc(k,1).Green(1);
                popc(k,2).Green(2)= popc(k,2).Cycle-L-popc(k,2).Green(4)-popc(k,2).Green(1);
        
                popc(k,1).Green(3)= popc(k,1).Green(1)+popc(k,1).Green(2);
                popc(k,2).Green(3)= popc(k,2).Green(1)+popc(k,2).Green(2); 
                
               
         else
             
                popc1 = popc(k,1).Cycle;
                popc(k,1).Cycle = popc(k,2).Cycle;
                popc(k,2).Cycle = popc1;
                
                popc(k,1).Green(4)= popc(k,1).Cycle-L-popc(k,1).Green(1)-popc(k,1).Green(2);
                popc(k,2).Green(4)= popc(k,2).Cycle-L-popc(k,2).Green(1)-popc(k,2).Green(2);
                       
         end
          
        % check if new green times are out constraints. 
        
        if( popc(k,1).Green(1)<y1*popc(k,1).Cycle || popc(k,1).Green(2)<y2*popc(k,1).Cycle || popc(k,1).Green(4)<y4*popc(k,1).Cycle)
            continue;
        end
        if( popc(k,2).Green(1)<y1*popc(k,2).Cycle || popc(k,2).Green(2)<y2*popc(k,2).Cycle || popc(k,2).Green(4)<y4*popc(k,2).Cycle)
            continue;
        end
        
        % Evaluate Generated Offsprings for each traffic light according to
        % the corresponding traffic congestion
        for j=1:nGL
            popc(k,1).TotalDelay(j)=FitnessFunction(popc(k,1).Cycle,popc(k,1).Green(j),Q(j),S(j),T);
            popc(k,2).TotalDelay(j)=FitnessFunction(popc(k,2).Cycle,popc(k,2).Green(j),Q(j),S(j),T);
            NFE=NFE+2;
        end
        
        % TOTAL DELAY which correspongs to the 4 lane groups    
        popc(k,1).TotalDelay= sum(Q.*popc(k,1).TotalDelay)/sum(Q);
        popc(k,2).TotalDelay= sum(Q.*popc(k,2).TotalDelay)/sum(Q);
        
        k=k+1; %step
    end
    
    % Make 2 columns into 1
    popc=popc(:);
    % Sort popc matrix according to TotalDelay
    TotalDelay=[popc.TotalDelay];
    [TotalDelay, SortOrder]=sort(TotalDelay);
    popc=popc(SortOrder);
    
    %% Mutation
    % Create empty Matrix with length the number of mutants 
    popm=repmat(empty_individual,nm,1);
    k=1;
    while k<=nm

        % Select Parent population
        i=randi([1 nPop]); 
        p=pop(i);
        
        % Apply Mutation   
        nVar=2;
        nmu=ceil(mu*nVar);

        j=randi([1 nVar]);
        prosimo=randi([-1 1]);
        sigma=prosimo*pm*(pop(i).Cycle);  
        
        mutated=p.Green(j)+sigma;
        
        popm(k).Green = p.Green;
        popm(k).Green(j)=mutated;
        popm(k).Cycle = p.Cycle;
        
        if j==1 
            popm(k).Green(2)= p.Green(2);
            popm(k).Green(3)= popm(k).Green(1)+popm(k).Green(2);
            popm(k).Green(4)= popm(k).Cycle-L-popm(k).Green(1)-popm(k).Green(2);
            
        else 
            popm(k).Green(1)= p.Green(1);
            popm(k).Green(3)= popm(k).Green(1)+popm(k).Green(2);
            popm(k).Green(4)= popm(k).Cycle-L-popm(k).Green(1)-popm(k).Green(2);
        end
        
        if( popm(k).Green(1)<y1*popm(k).Cycle || popm(k).Green(2)<y2*popm(k).Cycle || popm(k).Green(4)<y4*popm(k).Cycle)
            continue;
        end
        
        for j=1:nGL
            % Evaluate Mutant 
            popm(k).TotalDelay(j)=FitnessFunction(popm(k).Cycle,popm(k).Green(j),Q(j),S(j),T);
            NFE=NFE+1;
        end
        % Summation of delay quotients
        popm(k).TotalDelay= sum(Q.*popm(k).TotalDelay)/sum(Q);
        k=k+1; %step
    end
    %% Merge Population

    pop=[pop
        popc
        popm]; %#ok
     
    % Sort New Population according to TotalDelay
    TotalDelay=[pop.TotalDelay];
    [TotalDelay, SortOrder]=sort(TotalDelay);
    pop=pop(SortOrder);
    
    % Update Worst Cost
    WorstDelay=max(WorstDelay,pop(end).TotalDelay);
    
    % Keep the Best Population from the given number
    pop=pop(1:nPop);
    TotalDelay=TotalDelay(1:nPop);
    
    % Store Best Solution Ever Found
    BestSol=pop(1);
    
    % Store Best Cost Ever Found
    BestDelay(it)=BestSol.TotalDelay;
    
        % Show Iteration Information
    disp(['                   Iteration ' num2str(it) ': Best TotalDelay = ' num2str(BestDelay(it))]);
    fprintf('\n')
    disp('Green Timings:');
    fprintf('\n')
    disp(['  Phase 1 Green Time = ' num2str(BestSol.Green(1))'' ' seconds']);
    fprintf('\n')
    disp(['  Phase 2 Green Time = ' num2str(BestSol.Green(2))'' ' seconds']);
    fprintf('\n')
    disp(['  Phase 1/2 Green Time = ' num2str(BestSol.Green(3))'' ' seconds']);
    fprintf('\n')
    disp(['  Phase 3 Green Time = ' num2str(BestSol.Green(4))'' ' seconds']);
    fprintf('\n')
    NFE
    %end of generation
    it=it+1;
    
end

%% Results / Plots
figure(1);
plot(BestDelay,'LineWidth',2);
% plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('AvgDelay');
grid on;
toc %stop the timer