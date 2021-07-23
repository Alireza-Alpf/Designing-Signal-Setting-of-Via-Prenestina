clear
clc
tic;
%OptimalCycleTwoTerms_v06.M: M.File to analyze and minimize delays.
%The script applies a systematic search for optimal cycle length and green share by applying the two-term HCM delay formula.
%% Givens:
L=15.7; %[s] %Lost time [s]: Sum lost times for critical movements of different signal phases
T=0.25; %[hour] %Time Interval of analysis
QQ=[915,311,1671,369]; %[veh/h] %Flows for lane groups
SS=[3098,1417,3114,1807]; %[veh/h] %Saturation flow for lane groups

%% Set Parameters for Signal Setting Optimization
nu=50; %Number of steps for the analysis of green ratio
nC=100; %Number of steps for the analysis of cycle length
color1=['c'; 'g'; 'b'; 'm'; 'r'; 'k']; %Predefinition of colors for the curves of delay for different values of the cycle length
color2=['c--'; 'g--'; 'b--'; 'm--'; 'r--'; 'k--']; %Predefinition of colors for the curves of delay for different values of the cycle length
figuredelaygreenratio=1; %Set: figuredelaygreenratio=1, to plot the function of delay w.r.t. to the green ratio
figuredelaycycle=1; %Set: figuredelaygcycle=1, to plot the function of delay w.r.t. to the cycle length
figure3D=1; %Set: figure3D=1, to draw three-dimension contour and surface figures
nz=50; %Number of contour curves in the contour figure

%Preliminary computations
y1=[0.295]; %[adim] %Saturation degree of the critical movement for phase 1
y2=[0.242]; %[adim] %Saturation degree of the critical movement for phase 2
y3=[0.204]; %[adim] %Saturation degree of the critical movement for phase 3
Y=(y1+y2+y3); %[adim] %Total saturation degree
Cm=L/(1-Y); %[s] %Minimum Cycle length

%% Delay Computations
CC=linspace(1.01*Cm,5*Cm,nC); %Definition of the Vector of cycle length
    figure(1) %Open a new figure to plot the graph
    clf %Clear the graph
for j=1:nC %Loop for the cycle length
    C=CC(j); %[s] %Select one value of the cycle length
    g3=y3*C;
    u1min=y1/(1-y1-L/C-y3); %[adim] %Determine the lowest value for the green ratio g1/g2 (notice that it depends on the cycle length value)
    u1max=(1-y2-L/C-y3)/y2; %[adim] %Determine the highest value for the green ratio g1/g2 (notice that it depends on the cycle length value)
    u=linspace(1.02*u1min,0.98*u1max,nu); %[adim] %Definition of the Vector of green ratio g1/g2
    for i=1:nu %Loop for the green ratio g1/g2
        g1=u(i)./(1+u(i)).*(C-L-g3); %[s] %Determine the value of the green for phase 1. The following condition is applied: g1=u*g2, that is: g1=u*(C-L-g1), and: (1+u)*g1=u*(C-L).
        g2=C-L-g1-g3; %[s] %Determine the value of the green for phase 1. The green g2 is obtained as the complementary of g1 to C-L
        G=[g1 g2 g1+g2 g3]; %Reorganize the green times for each lane group
        d1=delay1(QQ,SS,G,C); %[s] %Compute the first term (deterministic) of delay per approach
        %d2=delay2(QQ,SS,G,C); %[s] %Compute the first term (probabilistic) of delay per approach
        d2=delay2(QQ,SS,G,C,T); %[s] %Compute the first term (probabilistic) of delay per approach
        d=d1+d2; %[s] %Determine the unitary delay (determinsitic plus probabilistic) per approach
        avgd=nansum(QQ.*d)/sum(QQ); %[s] %Average unitary junction delay. nansum(X) returns the sum of X, treating NaNs as missing values.
        dd(i,j)=avgd; %[s] %Fill the matrix of delay in function of C and green ratio (indexes i and j)
        [dminu(j),umin(j)]=nanmin(dd(:,j)); %Minimum delay dminu(j) w.r.t green ratio u for every value j-th of the cycle length; umin(j) is the position of the best green ratio u in the vector u(j)
    end %End of Green ratio variation
    
    %% Computations for the Graph of the the delay-green ratio function
    if figuredelaygreenratio==1 %Plot the curve with the lowest cycle in the delay-green ratio plane
        if j==2 %Check if the cycle is close to Cmin to draw a first curve
            plot(u,dd(:,2),color1(1,:)) %Plot the delay-green ratio function for a cycle very close to the minimum cycle
            plot(u(umin(2)),dminu(2),'*m'); %Plot for every curve the optimal point as a magenta asterisk ('*'))
            xlabel('Green Ratio'); ylabel('delay [s]') %Add x and y labels
            grid on %Plot a grid in the graph
        elseif rem(j/20,1)==0 %Condition to select the cycle length values (one every 20) whose corresponding unitary delay are to be plotted
            plot(u,dd(:,j),color1(j/20+1,:)) %Plot one delay-green ratio curve for every selected value of the cycle length
            plot(u(umin(j)),dminu(j),'*m'); %Plot for every curve the optimal point as a magenta asterisk ('*'))
        end %End of the condition to select the cycle length values whose corresponding unitary delay are to be plotted
        hold on %Keep the graphs previously computed within the loop for saturation degree variation
    end %End of the condition for a cycle close to Cmin
end %End of Cycle length variation

%% Find the Optimal Signal Settings
dmin=min(min(dd));  %[s] %Find the minimum average unitary junction delay
[uo,Co]=find(dd==dmin); %Find the positions of the optimal signal settings that minimize the average unitary junction delay
uopt=u(uo); %[adim] %Find the the optimal green ratio g1/g2
Copt=CC(Co); %[s] %Find the optimal cycle length
%% Display and plot the optimal signal settings and minimum delay-green ratio plane
disp('Optimal cycle=') %Command that displays a text in the Command window
Copt %[s] %Optimal Cycle length
disp('Optimal green ratio g1/g2=') %Command that displays a text in the Command window
uopt %[adim] %Optimal green ratio g1/g2
disp('Minimum unitary delay=') %Command that displays a text in the Command window
dmin %[s] %Mimimum average unitary delay of the junction
%% Draw 3D Graphs of the delay(green ratio, cycle length)
if figure3D==1 %If the option "figure3D" is equal to 1, one contour figure and a surface are plotted
    [GGG,CCC]=meshgrid(u,CC); %Fill a (nu,nC) grid of green with the values of ratio and cycle length (points of the independent variables in the 3D figure)
    figure(2) %Open a new figure to plot a contour figure
    clf %Clear the graph
    contour(u,CC,dd',nz,'k:') %Plot a contour figure with nz contour curves with a black color ('k')
    hold on %Keep the graphs previously computed within the loop for saturation degree variation
    xlabel('Green ratio'); ylabel('Cycle [s]')
    plot(uopt,Copt,'+m'); plot(uopt,Copt,'om') %Plot the optimal point in the green ratio-cycle plane %Plot the optimal point as a target ('+''o') in magenta color ('r'))
    figure(3) %Open a new figure to draw a 3D surface figure
    clf %Clear the graph
    mesh(u,CC,dd') %Draw a 3D surface figure in the space of green ratio, cycle length, and unitary junction delay
    xlabel('Green ratio'); ylabel('Cycle [s]'); zlabel ('delay [s]')
    grid on %Plot a grid in the graph
end %End of the condition for 3D surface figure

%Plot of delay as a function of cycle length with the green split as a
%parameter (code lines are very similar to the plot delay-green ratio
if figuredelaycycle==1 %If the option "figuredelaycycle" is equal to 1, a cycle-delay figure for differnt green ratio is plotted
    figure(4) %Open a new figure
    clf %Clear the graph
    p=1; %Index to count the curves with different values of green ratio that are plotted in the figure
    for k=1:nu %Loop for the green ratio changes
        plot(CC,dd(2,:),color1(1,:)) %Plot one delay-cycle curve for a small value of the green ratio for phase one
        if rem(k/10,1)==0 %Condition to select the green ratio values (one every 10) whose corresponding unitary delay are to be plotted
            plot(CC,dd(k,:),color1(k/10+1,:))  %Plot one delay-cycle curve for every selected value of the green ratio
            GreenRatio(p)=u(k); %Different values of green ratio that are plotted in the figure
            p=p+1; %Update the index p of the vector GreenRatio
            hold on %Holds the current graph to add new curves to the same figure
            grid on %Plot a grid in the graph
        end %End of the condition to select the curves to plot
    end %End of the loop for the green ratio changes
    %Plot the optimal signal setting in the delay-cycle plane
    plot(Copt,dmin,'*m'); %Plot the optimal point as a magenta cross
    xlabel('Cycle [s]'); ylabel('delay [s]') %Add x and y labels to the graph
end %End of the condition to plot the cycle-delay figure
toc