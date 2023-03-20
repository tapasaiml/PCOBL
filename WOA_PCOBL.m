%___________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) with Partial Centroid Opposition Based
%  Learning (PCOBL)                                                         %
%  Source codes                                                                             %
%  Developed in MATLAB R2018b                                               %
%                                                                           %
%  Author and programmer: Dr. Tapas Si                                      %
%                                                                           %
%         e-Mail: tapassi.aiml@gmail.com                                    %
%                                                                           %
%                                                                           %
%   %% Some code parts were used from orginal WOA code: 
%  Source link: https://seyedalimirjalili.com/woa                                                                        %
%                                                                           %
%   WOA paper: S. Mirjalili, A. Lewis                                      %
%               The Whale Optimization Algorithm,                           %
%               Advances in Engineering Software , in press,                %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008     %
%                                                                           %
%___________________________________________________________________________%
clear all;
clc;

fhd=str2func('cec13_func');

load('extremum2013');


N=30;
dim=input('Dimension:');


target_err=1e-8;%when target error is greater thn 10^-8 then prog will terminate.
TEST_RUN=51;
fes=zeros(TEST_RUN,1);%initialize fes by zeros 
best_it=zeros(TEST_RUN,1);%best iteration is a matrix where we store the best result

FES=10000*dim;

Max_iter=round(FES/N);

lb=-100;
ub=100;

Start_fun=input('Starting function#');
End_fun=input('Ending function#');
%Initialize the positions of moths
%Moth_pos=initialization(N,dim,ub,lb);

for func_num=Start_fun:1:End_fun

tic;
for test=1:1:TEST_RUN

rand('twister', sum(100*clock));
filename=sprintf('population\\%dD\\POPfun%d_run%d_D%d',dim,func_num,test,dim);    

load(filename);

SearchAgents_no=N;
% The Whale Optimization Algorithm
%function [Leader_score,Leader_pos,Convergence_curve]=WOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions=X;
fitcount=0;

for i=1:1:SearchAgents_no
   fitness(i)=feval( fhd, Positions(i,:)',func_num);
   fitcount=fitcount+1; 
    
end
dynXmin=min(Positions);
dynXmax=max(Positions);

centroid=mean(Positions);


for i=1:1:N
    
    Opp_X(i,:)=2.*centroid-Positions(i,:);
    
    for j=1:1:dim
        
        if Opp_X(i,j) > dynXmax(j) 
             
          d=centroid(j);
          c=dynXmax(j);
          Opp_X(i,j)=d+(c-d).*rand;
        end
        
        if Opp_X(i,j) < dynXmin(j) 
             
          d=dynXmin(j);
          c=centroid(j);
          Opp_X(i,j)=d+(c-d).*rand;
        end
        
     end
        
    Opp_Fitness(1,i)=feval( fhd,Opp_X(i,:)',func_num);
    fitcount=fitcount+1;
end

   
   AllFitness=[fitness Opp_Fitness];
   All_X=[Positions;  Opp_X];
   
   [sorted_fitness index]=sort(AllFitness);
   Positions=All_X(index(1:N),:);
   fitness=sorted_fitness(1:N);
   Leader_pos=Positions(1,:);
   Leader_score=fitness(1,1);
    
Convergence_curve=zeros(1,Max_iter);

t=0;% Loop counter

% Main loop
%while t<Max_iter

err=1;

num_pobl=3;

Pgj=0.3;

err=abs(extremum(func_num) - Leader_score);

while fitcount < FES && err > target_err 

  if rand < Pgj
        
        dynXmin=min(Positions);
        dynXmax=max(Positions);

        centroid=mean(Positions);

        index=0;  
for i=1:1:N
    
    Opp_X(i,:)=2.*centroid-Positions(i,:);
    
    for j=1:1:dim
        
        if Opp_X(i,j) > dynXmax(j) 
             
          d=centroid(j);
          c=dynXmax(j);
          Opp_X(i,j)=d+(c-d).*rand;
        end
        
        if Opp_X(i,j) < dynXmin(j) 
             
          d=dynXmin(j);
          c=centroid(j);
          Opp_X(i,j)=d+(c-d).*rand;
        end
        
     end
        
     Opp_Fitness(1,i)=feval( fhd,Opp_X(i,:)',func_num);
     fitcount=fitcount+1;
    
    
    
            for j=1:1:num_pobl
             %---degree of pobl
                degree=ceil(rand.*(dim-1));
                index=index+1;
                rand_idx=ceil(rand(1,degree).*dim);
                pobl_Positions(index,:)=Opp_X(i,:);
                pobl_Positions(index,rand_idx)=Positions(i,rand_idx);
                pobl_Fitness(1,index)=feval( fhd,pobl_Positions(index,:)',func_num);
                fitcount=fitcount+1;
              end
end

   
     AllFitness=[fitness Opp_Fitness pobl_Fitness];
     All_X=[Positions;  Opp_X; pobl_Positions];
    
   
    [sorted_fitness index]=sort(AllFitness);
    Positions=All_X(index(1:N),:);
    fitness=sorted_fitness(1:N);
    Leader_pos=Positions(1,:);
    Leader_score=fitness(1,1);
        
        
    else      
    
    
    
    
    a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                
            end
            
        end
    end
    
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        %fitness=fobj(Positions(i,:));
        fitness(i)=feval( fhd, Positions(i,:)',func_num);
        fitcount=fitcount+1;
        % Update the leader
        if fitness(i)<Leader_score % Change this to > for maximization problem
            Leader_score=fitness(i); % Update alpha
            Leader_pos=Positions(i,:);
        end
        
    end
    
    
    
  end
    t=t+1;
    Convergence_curve(t)=Leader_score;
    %[t Leader_score]
    
    
   err=abs(extremum(func_num) - Leader_score);  
    
end

fprintf('global best=%f error=%f\n',Leader_score,err);
 
 best_it(test,1)=err;
 
 
 fprintf('RUN %d is completed=================\n',test);
 
 fes(test,1)=fitcount;
end

cpu_time=toc;
T2=cpu_time/TEST_RUN

best_run=min( best_it)
worst_run=max( best_it )
mean_run=mean( best_it)
median_run=median(best_it);
std_run=std( best_it)

mean_fes=mean(fes)
std_fes=std(fes)

Succes_Ratio=(sum(best_it < 10e-8)/TEST_RUN)*100;
filename='WOA_PCOBL.xlsx';
sheet=sprintf('%dD',dim);
if func_num>=26
    column=sprintf('A%c','A'+func_num-26);
else
     column=sprintf('%c','A'+func_num);
end

range=sprintf('%s2:%s51',column,column);

xlswrite(filename,  best_it, sheet, range);
%xlswrite(filename,  best_it, sheet, range);

range=sprintf('%s53',column);
xlswrite(filename, best_run, sheet, range);

range=sprintf('%s54',column);
xlswrite(filename, worst_run, sheet, range);

range=sprintf('%s55',column);
xlswrite(filename, median_run, sheet, range);

range=sprintf('%s56',column);
xlswrite(filename, mean_run, sheet, range);

range=sprintf('%s57',column);
xlswrite(filename, std_run, sheet, range);

range=sprintf('%s58',column);
xlswrite(filename, mean_fes, sheet, range);

range=sprintf('%s59',column);
xlswrite(filename, std_fes, sheet, range);

range=sprintf('%s61',column);
xlswrite(filename, Succes_Ratio, sheet, range);


end%func

