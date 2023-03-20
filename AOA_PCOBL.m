%  Arithmetic Optimization Algorithm (AOA)  with Partial Centroid Opposition-Based
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
%  %% Some code parts were used from orginal AOA code: 
%  Source link: https://www.mathworks.com/matlabcentral/fileexchange/84742-the-arithmetic-optimization-algorithm-aoa
%  AOA paper:                                                                                        
%  Laith Abualigah, Ali Diabat, Seyedali Mirjalili, Mohamed Abd Elaziz, Amir H. Gandomi,
%The Arithmetic Optimization Algorithm,Computer Methods in Applied Mechanics and Engineering,
%Volume 376, 2021, 113609, ISSN 0045-7825, https://doi.org/10.1016/j.cma.2020.113609.
%(https://www.sciencedirect.com/science/article/pii/S0045782520307945)
%______________________________________________________________________________________________

close all; clear all; clc


LB=-100;
UB=100;

N=30;

dim=input('Dimension:');

fhd=str2func('cec13_func');

load('extremum2013');


target_err=1e-8;%when target error is greater thn 10^-8 then prog will terminate.
TEST_RUN=51;
fes=zeros(TEST_RUN,1);%initialize fes by zeros 
best_it=zeros(TEST_RUN,1);%best iteration is a matrix where we store the best result

FES=10000*dim;
%FES=200000;
M_Iter=round(FES/N);

Start_fun=input('Starting function#');
End_fun=input('Ending function#');
for func_num=Start_fun:1:End_fun

tic;
for test=1:1:TEST_RUN

rand('twister', sum(100*clock));
filename=sprintf('population\\%dD\\POPfun%d_run%d_D%d',dim,func_num,test,dim);    

load(filename);

%function [Best_FF,Best_P,Conv_curve]=AOA(N,M_Iter,LB,UB,Dim,F_obj)
%display('AOA Working');
%Two variables to keep the positions and the fitness value of the best-obtained solution

Best_P=zeros(1,dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);

%Initialize the positions of solution
%X=initialization(N,Dim,UB,LB);
Xnew=X;
Ffun=zeros(1,size(X,1));% (fitness values)
Ffun_new=zeros(1,size(Xnew,1));% (fitness values)

MOP_Max=1;
MOP_Min=0.2;
C_Iter=1;
Alpha=5;
Mu=0.499;


for i=1:size(X,1)
    Ffun(1,i)=feval(fhd, X(i,:)',func_num);    %Calculate the fitness values of solutions
    if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
    end
end
 
fitcount=N;
 
dynXmin=min(X);
dynXmax=max(X);

centroid=mean(X);
for i=1:1:N
    %----Calculate opposite 

    
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
       
    Opp_Objective_values(1,i) = feval( fhd,Opp_X(i,:)',func_num);
    fitcount=fitcount+1; 
 end
   AllFitness=[Ffun Opp_Objective_values];
   All_X=[X; Opp_X];
   
   [sorted_fitness index]=sort(AllFitness);
   X=All_X(index(1:N),:);
   Ffun=sorted_fitness(1:N);
   Best_P=X(1,:);
   Best_FF=Ffun(1,1);
C_Iter=1;
%--Diversity------------
    mean_X=median(X);
  
    for i=1:1:N
       
       sum1(i)=mean(abs(X(i,:)-mean_X));
       
   end
  
   diversity(C_Iter)=sum(sum1)/N;

 
Pgj1=0.3;
Pgj2=Pgj1;
eta=1.0;
mu=1;
sigma=0.25;


num_pobl=3;
err=abs(extremum(func_num) - Best_FF);
%while C_Iter<M_Iter+1  %Main loop
 while fitcount < FES && err > target_err 
 

     
  if rand < Pgj2
   
  A=min(X);
  B=max(X);    
      
     
   index=0;  
    for i=1:1:N
     if fitcount < FES
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
       
    Opp_Fitness(1,i) = feval( fhd,Opp_X(i,:)',func_num);
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
        
    end

   
   AllFitness=[Ffun Opp_Fitness pobl_Fitness];
   All_X=[X;  Opp_X ; pobl_Positions];

   
   
   
   [sorted_fitness index]=sort(AllFitness);
   X=All_X(index(1:N),:);
   Ffun=sorted_fitness(1:N);
   Best_P=X(1,:);
   Best_FF=Ffun(1,1);
    
   Xnew=X;   
   
      
  else    
     
     
     
    MOP=1-((C_Iter)^(1/Alpha)/(M_Iter)^(1/Alpha));   % Probability Ratio 
    MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter); %Accelerated function
   
    %Update the Position of solutions
    for i=1:size(X,1)   % if each of the UB and LB has a just value 
        for j=1:size(X,2)
           r1=rand();
            if (size(LB,2)==1)
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=Best_P(1,j)/(MOP+eps)*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,j)=Best_P(1,j)*MOP*((UB-LB)*Mu+LB);
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=Best_P(1,j)-MOP*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,j)=Best_P(1,j)+MOP*((UB-LB)*Mu+LB);
                    end
                end               
            end
            
           
            if (size(LB,2)~=1)   % if each of the UB and LB has more than one value 
                r1=rand();
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=Best_P(1,j)/(MOP+eps)*((UB(j)-LB(j))*Mu+LB(j));
                    else
                        Xnew(i,j)=Best_P(1,j)*MOP*((UB(j)-LB(j))*Mu+LB(j));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=Best_P(1,j)-MOP*((UB(j)-LB(j))*Mu+LB(j));
                    else
                        Xnew(i,j)=Best_P(1,j)+MOP*((UB(j)-LB(j))*Mu+LB(j));
                    end
                end               
            end
            
        end
        
        Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
 
        Ffun_new(1,i)=feval(fhd, Xnew(i,:)',func_num); % calculate Fitness function 
        fitcount=fitcount+1;
        if Ffun_new(1,i)<Ffun(1,i)
            X(i,:)=Xnew(i,:);
            Ffun(1,i)=Ffun_new(1,i);
        end
        if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
        end
       
    end
  end

  
  
    %Update the convergence curve
    Conv_curve(C_Iter)=Best_FF;
    C_Iter=C_Iter+1;  % incremental iteration
    %--Diversity------------
    mean_X=median(X);
  
    for i=1:1:N
       
       sum1(i)=mean(abs(X(i,:)-mean_X));
       
   end
  
   diversity(C_Iter)=sum(sum1)/N;
    %Print the best solution details after every 50 iterations
%     if mod(C_Iter,50)==0
%         display(['At iteration ', num2str(C_Iter), ' the best solution fitness is ', num2str(Best_FF)]);
%     end
%      
    
   
   err=abs(extremum(func_num) - Best_FF);
    
end


fprintf('global best=%f error=%f\n',Best_FF,err);
 
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
filename='AOA_PCOBL.xlsx';
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

% divmax=max(diversity);
% AOA_PCOBL_XPL=(diversity./divmax).*100;
% AOA_PCOBL_XPT=(abs(diversity-divmax)./divmax).*100;
% trade_off=sqrt(AOA_PCOBL_XPL.*AOA_PCOBL_XPT);
% avg_AOA_PCOBL_XPL(func_num,1)=mean(AOA_PCOBL_XPL);
% avg_AOA_PCOBL_XPT(func_num,1)=mean(AOA_PCOBL_XPT);
% avg_trade_off(func_num,1)=mean(trade_off);
% AOA_pcobl=diversity;
%  fname=sprintf('XPL_XPT\\AOA_PCOBL_DIV_GRAPH_Func%d_%dD', func_num,dim);
%  save(fname,'AOA_pcobl','AOA_PCOBL_XPL','AOA_PCOBL_XPT');


end

