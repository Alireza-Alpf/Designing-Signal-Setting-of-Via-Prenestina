%Enumerative method: M.File to analyze and minimize delays (HCM)
%The script applies a systematic search for optimal cycle length and green share by applying the two-term HCM delay formula.
tic;
%% Givens:
L=15.7; %[s] %Lost time [s]: Sum lost times for critical movements of different signal phases
T=0.25; %[hour] %Time Interval of analysis
temp=0; %This is a variable for storing the data
Q=[915,311,1671,369]; %[veh/h] %Flows for lane groups EB,WB,NB,SB
S=[3098,1417,3114,1807]; %[veh/h] %Saturation flow corresponding to the lane groups

%Preliminary computations
y1=0.295; %[adim] %Saturation degree of the critical movement for phase 1
y2=0.242; %[adim] %Saturation degree of the critical movement for phase 2
y3=0.204; %[adim] %Saturation degree of the critical movement for phase 3
Y=(y1+y2+y3); %[adim] %Total saturation degree
Cm=L/(1-Y); %[s] %Minimum Cycle length
%% Delay Computations
CC=Cm:1:5*Cm; %Definition of the Vector of cycle length (Here you can put the step size of cycle length for every iteration)
nc=numel(CC); %Number of elements in CC
for j=1:nc %Loop for the cycle length
    C=CC(j); %[s] %Select one value of the cycle length
    ming1=y1*C; ming2=y2*C; ming3=y3*C; % min green for every phase
    maxg1=C-L-ming2-ming3; %max green for phase 1 (Here you can put the step size of green for phase 1 for every iteration)
    maxg2=C-L-ming1-ming3; %max green for phase 2 (Here you can put the step size of green for phase 2 for every iteration)
    g1=ming1:1:maxg1; ng1=numel(g1); % definition of the vector of g1 and Number of elements in g1
    g2=ming2:1:maxg2; ng2=numel(g2); % definition of the vector of g2 and Number of elements in g2
    for i=1:ng1 %Loop for the green ratio g1/C
        for k=1:ng2 %Loop for the green ratio g2/C
            g3(k)=C-L-g1(i)-g2(k);  %calculating g3 based on g1 and g2  
           if( g3(k)<ming3 )   %If the value for green phase 3 become less than min the inner for loop should stop
               break
           end
           temp=temp+1;
        g(:,temp)=[g1(i);g2(k);g1(i)+g2(k);g3(k);C];%[s] %Collect the vector of greens corresponding to the lane groups
       
        for m=1:4   %Because we have 4 lane groups it should iterate 4 times to calculate the delays
        X=(Q(m)*C)/(g(m,temp)*S(m));  %Computing the V/C ratio
        c=(S(m)*g(m,temp))/C;   %Computing the lane group capacity
        d1(m,k)=(0.5*C*(1-(g(m,temp)/C))^2)/(1-(min(1,X)*(g(m,temp)/C))); %[s] %Compute the first term HCM Delay
        d2(m,k)=900*T*((X-1)+sqrt((X-1)^2+((8*0.5*X)/(c*T)))); %[s] %Compute the second term HCM Delay
        d(m,k)=d1(m,k)+d2(m,k); %Average delay for every lane group
        end
        avgd(:,temp)=sum(Q*d(:,k))/sum(Q); %here we store the average delay of every iteration for intersection
        end 
    end
end
[minDelay,idx]=min(avgd); %Display the Min delay
mD=['Min Delay = ',num2str(minDelay)];
disp(mD)
g1=['green for phase 1 = ',num2str(g(1,idx))];  %here you should put the row index corresponding to green phase 1 of green vector
disp(g1)
g2=['green for phase 2 = ',num2str(g(2,idx))];  %here you should put the row index corresponding to green phase 2 of green vector
disp(g2)
g3=['green for phase 3 = ',num2str(g(4,idx))];  %here you should put the row index corresponding to green phase 3 of green vector
disp(g3)
C=['Cycle length = ',num2str(g(5,idx))];  %here you should put the row index corresponding to cycle length of green vector
disp(C)
toc