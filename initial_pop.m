%---Generate initial Population
function initial_pop()
clear all
clc;
D_POP=100;% Dimension of each individual
N_POP=30;

Xmin=-100;
Xmax=100; 


for func_No=1:1:28
for run=1:1:51
rand('twister', sum(100*clock));


X=zeros(N_POP,D_POP);
X=Xmin+(Xmax-Xmin).*rand(N_POP,D_POP);
filename=sprintf('population\\%dD\\POPfun%d_run%d_D%d',D_POP,func_No,run,D_POP);
save(filename, 'X');

end

end