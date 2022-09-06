%% Modeling the Spread of an Infection through a Small Community
N = 5000; %population
Days = 50; %number of simulation days
mI = 4; %number of days between infection and being infectious
mS = 4; %number of days between infection and showing symptoms (sick)
mR = 14; %number of days to recovery/death (immune)
pd = 0.03; %probability of death
pIO = 0.02; %fraction of N infected at t=0
pI = 0.085; %probability of infection transmission

homemean=3;
homesdv=1;
stayshome = 100*rand(N,1)<=80; %Logical vector, true if person stays home when 

schoolmean=5;
schoolsdv=3;
goestoschool = 100*rand(N,1)<=17; %Logical vector true if person attends school

workmean=5;
worksdv=3;
goestowork = 100*rand(N,1)<=45; %Logical vector true if person works

smallrecmean=4;
smallrecsdv=2;
insmallgroup = 100*rand(N,1)<=80; %Logical vector true if person is part of small group

commutemean=10;
commutesdv=5;
commutes=100*rand(N,1)<=30; %Logical vector true if person commutes

largerecmean=10;
largerecsdv=10;
inlargegroup=100*rand(N,1)<=30; %Logical vector true if attends large group gatherings

errandmean=5;
errandsdv=4;
runserrands=100*rand(N,1)<=14.3; %Logical vector true if person does errands on a given day

%Data Initialization
infectionday = randi(50,Days,1);%integer vector specifying day each person infected
infected = rand(N,1)<=pIO;%Logical vector identifying infected people
sick = false(N,1);%Logical vector identifying people with symptoms
immune = false(N,1); %Logical vector identifying recovered/immune people
dead = false(N,1); %Logical vector identifying people who have died

%Behavior Initialization
stayshome=initialize_trait(stayshome,N); %Members who stay home when sick (are good citizens)
commutes=initialize_trait(commutes,N); %Members who commute
inlargegroup=initialize_trait(inlargegroup,N); %Members who attend large recreational groups

homegroups=generate_groups(homemean,homesdv,stayshome,N);%integer vector specifying home # occupied by each person

goestoschool=initialize_trait(goestoschool,N);
schoolgroups=generate_groups(schoolmean,schoolsdv,goestoschool,N);%integer vector specifying a school group for each person

goestowork=initialize_trait(goestowork,N);
workgroups=generate_groups(workmean,worksdv,goestowork,N); %integer vector specifying a work group for each person

insmallgroup=initialize_trait(insmallgroup,N);
smallrecgroups=generate_groups(smallrecmean,smallrecsdv,insmallgroup,N); %integer vector specifying small social groups

infectedatschool = zeros(Days,1);
infectedatwork = zeros(Days,1);
infectedathome = zeros(Days,1);
infectedamiderrands = zeros(Days,1);
infectedamidsmallrecgroup = zeros(Days,1);
infectedamidlargerecgroup = zeros(Days,1);
infectedamidcommute = zeros(Days,1);

n_infected = zeros(Days,1);
n_dead = zeros(Days,1);
n_sick = zeros(Days,1);
n_immune = zeros(Days,1);




for day=1:Days
    day;
    if (commutes>0)
        commutergroups=generate_groups(commutemean,commutesdv,commutes,N);%commutegroups: integer vector - groups in contact during commute
      
        ni = sum(infected);
        [infected,infectionday] = group_progression(commutergroups,infected,sick,immune,infectionday,day,stayshome,pI,mI);
        infectedamidcommute(day) = sum(infected)-ni;
    end 
    
    if (goestoschool>0)
        ni = sum(infected);
        [infected,infectionday] = group_progression(schoolgroups,infected,sick,immune,infectionday,day,stayshome,pI,mI);
        infectedatschool(day) = sum(infected)-ni;
    end
    
    if (goestowork>0)
        ni = sum(infected);
        [infected,infectionday] = group_progression(workgroups,infected,sick,immune,infectionday,day,stayshome,pI,mI);
        infectedatwork(day) = sum(infected)-ni;
    end 

    
    runserrands = initialize_trait(round(N/7),N); %Select at random 1/7 of N that run errands on any given day
    errandgroups = generate_groups(errandmean,errandsdv,runserrands,N); %integer vector - groups in contact during errands
    
    ni = sum(infected);
    [infected,infectionday] = group_progression(errandgroups,infected,sick,immune,infectionday,day,stayshome,pI,mI);
    infectedamiderrands(day) = sum(infected)-ni;

    
    if (insmallgroup>0)
       ni = sum(infected);
       [infected,infectionday] = group_progression(smallrecgroups,infected,sick,immune,infectionday,day,stayshome,pI,mI);
       infectedamidsmallrecgroup(day) = sum(infected)-ni;
    end 


    if (inlargegroup>0)
        largerecgroups=generate_groups(largerecmean,largerecsdv,inlargegroup,N); %integer vector specifying large social group attendance
        
        ni = sum(infected);
        [infected,infectionday] = group_progression(largerecgroups,infected,sick,immune,infectionday,day,stayshome,pI,mI);
        infectedamidlargerecgroup(day) = sum(infected)-ni;
    end 

    ni = sum(infected);
    [infected,infectionday] = group_progression(homegroups,infected,sick,immune,infectionday,day,stayshome,pI,mI);
    infectedathome(day) = sum(infected)-ni;
    
    [infected,sick,immune,dead] = disease_progression(infected,sick,immune,dead,infectionday,day,mS,mR,pd,N);
    mask = (day-infectionday == mS);
    sick_today=day-infectionday >= mS;
    sick(mask)=true;
    n_sick(day) = sum(sick);
    mask = (day-infectionday == mR);
    sick_or_dead=day-infectionday >= mR;
    immune(mask) = true;
    n_immune(day) = sum(immune);
    dead(mask)=true;
    n_dead(day)=sum(dead);
 
end

ni(day)=sum(infected);

days = 1:day;
plot(days,n_immune)

hold on;
plot(days,n_sick)
hold on
plot(days,n_dead)
legend('People Immune per day','People Sick per day',"People dead per day")

hold off
%average_sick=mean(n_sick)

%n_immune-n_dead) / N

%plot(n_dead/N)
%hold off


   
function [infected_new,infectionday_new] = group_progression(groupnumbers,infected,sick,immune,infectionday,day,stayshome,pI,mI) %Function 1
n_groups = max(groupnumbers);
infected_new = infected;
infectionday_new = infectionday;
rollthedice=randi(100,1); %Decides whether a person coming into contact with an infected person gets infected
canbeinfected = rollthedice<pI; %Selects people who get infected
infectious = (day - mI >= 0) & infected & ~(sick & stayshome);%Logical vector that identifies whether someone is infectious
for i = 1:n_groups
    mask = groupnumbers==i;
    if (sum(mask & infectious)>0)
        mask2 = mask & ~infected & canbeinfected & ~immune;
        infected_new(mask2) = true;
        infectionday_new(mask2) = day;
    end
end 
infected_new(infected)=true;
%infectionday_new(infected)=infectionday(infected);
end 

function [infected_new,sick_new,immune_new,dead_new]= disease_progression(infected,sick,immune,dead,infectionday,day,mS,mR,pd,N) %Function 2
infected_new = infected;
sick_new = sick;
immune_new = immune;
dead_new = dead;
mask = (day-infectionday == mS);
sick_new(mask)=true;
mask = (day-infectionday == mR);
immune_new(mask)=true;
sick_new(mask)=false;
infected_new(mask)=false;
n = sum(mask==true); %Counting number of members
 if (n>0)
    rollthedice = 100*ones(N,1);
    rollthedice(mask)=randi(100,n,1);
    mask2 = rollthedice<pd;
    dead_new(mask2) = true;
    infected_new(mask2) = false; 
 end 
end


function flagtraits = initialize_trait(percentwithattribute,N)
choosemembers = randi(100,N,1);
flagtraits = choosemembers<=percentwithattribute;
end 

function groupnumbers = generate_groups(groupmeansize,groupsdv,participation,population)
% Function to sort population into distinct groups
% The size of the groups is approximately Gaussian, with prescribed mean and standard deviation
%       groupmeansize = average size of the groups
%       groupsdv = standard deviation of group size
%       participation(i=1:population) Logical flag indicating whether ith member of the population should be placed in a group
%       Population = total # people in the population
%
%       groupnumbers(i)   The group occupied by the ithe member of the population (zero if member is not in a group)
%
nmembers = sum(participation); % Number of people to be placed in groups
ngroups = round(nmembers/groupmeansize); % Number of groups
groupsize = round(normrnd(groupmeansize,groupsdv,ngroups,1)); % group sizes (approx Gaussian distribution)
mask = groupsize>0;
groupsize = groupsize(mask); % Remove groups with zero or fewer members (distribution no longer perfectly Gaussian)
groupnumbers = zeros(population,1); % Initialize group number for whole population
nassigned = 0;
assignmentorder = 3*population*ones(population,1); % Initialize vector used to place participants in queue
assignmentorder(participation) = randi(2*population,nmembers,1); % Each member of the population is put in the queue in random order
[~,index] = sort(assignmentorder);   % Generate index table pointing to group members in order of priority
% Fill the groups
for i = max(groupsize):-1:1
   availablegrouplist = find(groupsize>=i);  % Add a member to all groups with size >= i
   navailablegroups = length(availablegrouplist); % No. groups to receive a new member
   if (navailablegroups<1) break; end % Break out of loop if no groups left
   groupnumbers(index(nassigned+1:min(nassigned+navailablegroups,nmembers))) = availablegrouplist(1:min(navailablegroups,nmembers-nassigned)); % Assign people to groups
   nassigned = nassigned + min(navailablegroups,nmembers-nassigned);
   if (nassigned >= nmembers) break; end % Break out of loop if everyone has been assigned to a group
end


end

function [adj] = update_adjacency(groupnumbers,adj)

n_groups = max(groupnumbers);
for i=1:n_groups
    indices = find(groupnumbers==i);
    adj(indices,indices) = true;
end
% Remove self interactions
[n,~] = size(adj);
adj(1:n+1:end) = false;

end