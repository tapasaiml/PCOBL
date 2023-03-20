%  Sine Cosine Algorithm (SCA)  with Partial Centroid Opposition Based
%  Learning (PCOBL)
%
%  Source codes                                                                     
%                                                                                                     
%  Developed in MATLAB R2018b                                                                
%                                                                                                     
%  Author and programmer: Dr. Tapas Si                                                          
%                                                                                                     
%         e-Mail: tapassi.aiml@gmail.com                                                             
%                                                                
%                                                                                                     
%                                                              
%  %% Some code parts were used from orginal SCA code: 
%  Source link: https://seyedalimirjalili.com/sca
%  SCA paper:                                                                                        
%  S. Mirjalili, SCA: A Sine Cosine Algorithm for solving optimization problems
%  Knowledge-Based Systems, DOI: http://dx.doi.org/10.1016/j.knosys.2015.12.022
%______________________________________________________________________________________________
clear all;
clc;

fhd=str2func('cec13_func');

load('extremum2013');


N=30;
dim=input('Dimension:');


target_err=1e-8;%when target error is greater than 10^-8 then prog will terminate.
TEST_RUN=51;
fes=zeros(TEST_RUN,1);%initialize fes by zeros 
best_it=zeros(TEST_RUN,1);%best iteration is a matrix where we store the best result

FES=10000*dim;
%FES=200000;
Max_iteration=round(FES/N);

lb=-100;
ub=100;

Start_fun=input('Starting function#');
End_fun=input('Ending function#');

%disp('SCA is now tackling your problem')

for func_num=Start_fun:1:End_fun

tic;
for test=1:1:TEST_RUN

rand('twister', sum(100*clock));
filename=sprintf('population\\%dD\\POPfun%d_run%d_D%d',dim,func_num,test,dim);    

load(filename);



Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));

% Calculate the fitness of the first set and find the best one
fitcount=0;
for i=1:size(X,1)
    Objective_values(1,i)=feval(fhd, X(i,:)',func_num);
    fitcount=fitcount+1;
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
     
end
  
dynXmin=min(X);
dynXmax=max(X);

centroid=mean(X);


for i=1:1:N
    
    Opp_X(i,:)=2.*centroid-X(i,:);
    
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

   
   AllFitness=[Objective_values Opp_Fitness];
   All_X=[X;  Opp_X];
   
   [sorted_fitness index]=sort(AllFitness);
   X=All_X(index(1:N),:);
   AllFitness=sorted_fitness(1:N);
    Destination_position=X(1,:);
    Destination_fitness=AllFitness(1,1);
    
    %All_objective_values(1,i)=Objective_values(1,i);


%Main loop
t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitnes
err=abs(extremum(func_num) - Destination_fitness);
Pj=0.2;
num_pobl=3;
while fitcount < FES && err > target_err 
     if rand < Pj
        dynXmin=min(X);
        dynXmax=max(X);

        centroid=mean(X);

        index=0;  

for i=1:1:N
    
    Opp_X(i,:)=2.*centroid-X(i,:);
    
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
                pobl_Positions(index,rand_idx)=X(i,rand_idx);
                pobl_Fitness(1,index)=feval( fhd,pobl_Positions(index,:)',func_num);
                fitcount=fitcount+1;
            end
    
    
    
    
end


   
   AllFitness=[Objective_values Opp_Fitness pobl_Fitness];
   All_X=[X;  Opp_X ; pobl_Positions];
   
   [sorted_fitness index]=sort(AllFitness);
   X=All_X(index(1:N),:);
   AllFitness=sorted_fitness(1:N);
   Destination_position=X(1,:);
   Destination_fitness=AllFitness(1,1);   
     
     else
    % Eq. (3.4)
    a = 2;
%     Max_iteration = Max_iteration;
   % r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
    r1=a-fitcount*((a)/FES);
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            
            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
            
        end
    end   
    
    for i=1:size(X,1)
         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %AllFitness(i) = feval( fhd,X(i,:)',func_num);
        
        % Calculate the objective values
        Objective_values(1,i)=feval( fhd,X(i,:)',func_num);
        fitcount=fitcount+1;
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    
    end
        
   end
    
    
    
    
    % Display the iteration and best optimum obtained so far
%     if mod(t,50)==0
%         display(['At iteration ', num2str(t), ' the optimum is ', num2str(Destination_fitness)]);
%     end
%     
    Convergence_curve(t)=Destination_fitness;
    % Increase the iteration counter
    t=t+1;
    
     err=abs(extremum(func_num) - Destination_fitness);
end
 
   
    
   fprintf('global best=%f error=%f\n',Destination_fitness,err);
 
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

Succes_Ratio=(sum(best_it < 10e-8)/50)*100;
filename='SCAPCOBL.xlsx';
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

