
clear all
close all
clc
tic  % Start the timer

L=15.7; %[s] %Lost time [s]: Sum lost times for critical movements of different signal phases
T=0.25; %[hour] %Time Interval of analysis
Q=[915,311,1671,369]; %[veh/h] %Flows for lane groups
S=[3098,1417,3114,1807]; %[veh/h] %Saturation flow for lane groups
y1=[0.295]; %[adim] %Saturation degree of the critical movement for phase 1
y2=[0.242]; %[adim] %Saturation degree of the critical movement for phase 2
y3=[0.204]; %[adim] %Saturation degree of the critical movement for phase 3
Y=(y1+y2+y3); %[adim] %Total saturation degree
Cm=L/(1-Y); %[s] %Minimum Cycle length
Cyclemin=Cm;                                 % Lower bound of CYCLE
Cyclemax=5*Cm;                               % Upper bound of CYCLE

initial =[] ; % initial values of cycle length and green ratio
initial(1)= Cyclemin+(Cyclemax-Cyclemin)*rand; %Generate random number for Cycle
u1min=y1/(1-y1-L/initial(1)-y3); %[adim] %Determine the lowest value for the green ratio g1/g2 (notice that it depends on the cycle length value)
u1max=(1-y2-L/initial(1)-y3)/y2; %[adim] %Determine the highest value for the green ratio g1/g2 (notice that it depends on the cycle length value)
initial(2)=u1min+(u1max-u1min)*rand; %Generate random number for green ratio
y3=0.212; %y3 calibrated to get the best result 
g3=y3*initial(1);  %Calculating the green of phase 3
g1=initial(2)./(1+initial(2)).*(initial(1)-L); %Calculating the green of phase 1
g2=initial(1)-L-g1-g3; %Calculating the green of phase 2
G=[g1 g2 g1+g2 g3];  % Define a vector of green corresponding to group lanes
%Calculating the delay for initial values
for i=1:4
iDelay(i) = objectiveFunction(initial(1), G(i),Q(i),S(i),T);
end
initialDelay=sum(Q.*iDelay)/sum(Q); %initial Delay

stepSize = [0.1 , 0.01];  %Define the size of steps  

iteration = 1;
value.delay=[]; %Define a structure for storing delay of every iteration
value.greenratio=[]; %Define a structure for storing green ratio of every iteration
value.cycle=[]; %Define a structure for storing cycle of every iteration
bestcycle=initial(1);
bestgreenratio=initial(2);
minDelay1=initialDelay; %Define variable for while function when we are calculating neighbors of greenratio
minDelay2=initialDelay; %Define variable for while function when we are calculating neighbors of cycle
Neighbours=zeros(2,2); %Define a vector for neighbors
improvement1 = 0;
    while improvement1 == 0
       
        Neighbours=[bestcycle bestgreenratio+stepSize(2)    %Defining neighbours
                    bestcycle bestgreenratio-stepSize(2)
                    bestcycle+stepSize(1) bestgreenratio
                    bestcycle-stepSize(1) bestgreenratio
                    bestcycle+stepSize(1) bestgreenratio+stepSize(2)
                    bestcycle+stepSize(1) bestgreenratio-stepSize(2)
                    bestcycle-stepSize(1) bestgreenratio-stepSize(2)
                    bestcycle-stepSize(1) bestgreenratio+stepSize(2)];
        g3=y3*Neighbours(:,1);  %Calculating the green of phase 3
        g1=Neighbours(:,2)./(1+Neighbours(:,2)).*(Neighbours(:,1)-L);  %Calculating the green of phase 1
        g2=Neighbours(:,1)-L-g1-g3;  %Calculating the green of phase 2
        G=[g1(1) g2(1) g1(1)+g2(1) g3(1)
           g1(2) g2(2) g1(2)+g2(2) g3(2)
           g1(3) g2(3) g1(3)+g2(3) g3(3)
           g1(4) g2(4) g1(4)+g2(4) g3(4)
           g1(5) g2(5) g1(5)+g2(5) g3(5)
           g1(6) g2(6) g1(6)+g2(6) g3(6)
           g1(7) g2(7) g1(7)+g2(7) g3(7)
           g1(8) g2(8) g1(8)+g2(8) g3(8)];
       %Calculating the delay 
        for j=1:8
            for i=1:4
           iDelay(i) = objectiveFunction(Neighbours(j,1),G(j,i),Q(i),S(i),T); 
            end
           value(iteration).delay(j)= sum(Q.*iDelay)/sum(Q);
           value(iteration).greenratio(j)=Neighbours(j,2);
           value(iteration).cycle(j)=Neighbours(j,1);
        end
        for k=1:8
           if (value(iteration).delay(k)<minDelay1)
            minDelay1=value(iteration).delay(k)
            bestcycle=Neighbours(k,1);
            bestgreenratio=Neighbours(k,2);
           end
        end
            
        if (value(iteration).delay(1)>minDelay1 && value(iteration).delay(2)>minDelay1 && value(iteration).delay(3)>minDelay1 && value(iteration).delay(4)>minDelay1 && value(iteration).delay(5)>minDelay1 && value(iteration).delay(6)>minDelay1 && value(iteration).delay(7)>minDelay1 && value(iteration).delay(8)>minDelay1)
            improvement1=1;
        end
        bestDelay(iteration)=minDelay1;
        iteration=iteration+1;
    end
%% Results / Plots
figure(1);
plot(bestDelay,'LineWidth',2);
% plot(bestDelay,'LineWidth',2);
xlabel('Iteration');
ylabel('AvgDelay');
grid on;    
toc       % Stop the timer