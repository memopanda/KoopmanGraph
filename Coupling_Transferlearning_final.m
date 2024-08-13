clear all
close all
% for examples from section 6
% This code is adapted from codes used in the paper M. Korda, M. Putinar, I. Mezic (2020)
addpath('Resources\')
addpath('Datensatz\');
addpath('..\..\Toolbox_share\');
% rng shuffle
rng(2141444)
%% preprocessing(ensure X0 cent_f cent cent_c Nrbf are same)

example= '3'; % 1: learning between identical subsystems 2: learning between partically identical subsystems 3: learning within identical & partially subsystems

input_ind = 1;  % local input index of 1-th agent for 2-th agent
input_ind_2 = 2;  % local input index of 3-th agent for 4-th agent
n = 6; %num of full states of orginal system
N=3;   % num of agents

switch example
    case '1'
        N_n = N+1;
        input_ind_transfer =  [1];  % global input index  for 3-th agent  
    case '2'
        N_n = N;
        input_ind_transfer =  [1 4];  % global input index  for 3-th agent  
        input_ind_new_trans = 2;  % local input index of 4-th agent for 3-th agent 
    case '3'
        N_n = N+1;
        input_ind_transfer =  [1 8];  % global input index  for 3-th agent  
        input_ind_new_trans = 2;   %local input index of 4-th agent for 3-th agent, which is computed from global input index .
end
 
printing=1;
folder='Figs/';
%% *************************** Dynamics ***********************************

mu1 = 0.2;
lamda1 = 0; 
lamda1_2 = 0.1;
mu2 =  0.06;  
lamda2 = 0; 
lamda2_2 = 0.08;
beta2 =  0.05;  
mu3 = 0.004; %
lamda3 = 0; 
lamda3_2 = 0.03; 
beta3 = 0.001; 
ga3 =0.08;    

q=3;


f_0 =  @(t,x)([ mu1*x(2,:); lamda1*x(1,:)-lamda1_2*x(1,:).^3; ...
        mu2*x(4,:); lamda2*x(3,:)-lamda2_2*x(3,:).^q+ beta2*(x(4,:).^0).*x(input_ind,:); ...
        mu3*x(6,:); lamda3*x(5,:)-lamda3_2*x(5,:).^q + beta3*(x(6,:).^0).*x(input_ind_transfer(1),:)] ); 

if length(input_ind_transfer)>1
id_input_3 = ceil(input_ind_transfer(2)*N/n);  % transfer global index to index of subsystem
% system with N=4 agents (for partialy identical + (identical) transfer
% learning). Note: for partially identical case, the subsystem 4 is
% connected but does not influence the other 3 subsystems due to topology
f_0_new =  @(t,x)([mu1*x(2,:); lamda1*x(1,:)-lamda1_2*x(1,:).^3; ...
        mu2*x(4,:); lamda2*x(3,:)-lamda2_2*x(3,:).^q+ beta2*(x(4,:).^0).*x(input_ind,:); ...
         mu3*x(6,:); lamda3*x(5,:)-lamda3_2*x(5,:).^q + beta3*(x(6,:).^0).*x(input_ind_transfer(1),:) + ga3*(x(6,:).^1).*x(input_ind_transfer(2),:);...
        mu2*x(8,:); lamda2*x(7,:)-lamda2_2*x(7,:).^q+ beta2*(x(8,:).^0).*x(4+input_ind_2,:)] ); 
N_o=4;
else 
    id_input_3 = [];
    % system with N=4 agents (for identical transfer learning)
    f_0_new =  @(t,x)([mu1*x(2,:); lamda1*x(1,:)-lamda1_2*x(1,:).^3; ...
        mu2*x(4,:); lamda2*x(3,:)-lamda2_2*x(3,:).^q+ beta2*(x(4,:).^0).*x(input_ind,:); ...
         mu3*x(6,:); lamda3*x(5,:)-lamda3_2*x(5,:).^q + beta3*(x(6,:).^0).*x(input_ind_transfer(1),:);...
        mu2*x(8,:); lamda2*x(7,:)-lamda2_2*x(7,:).^q+ beta2*(x(8,:).^0).*x(4+input_ind_2,:)] ); 
    N_o=4;
end



%% ************************** Discretization ******************************

deltaT = 0.01;
deltaTg = 0.005;
%Runge-Kutta 4

f_0d = discrete_runge_kutta4(f_0, deltaT);
f_0dg = discrete_runge_kutta4(f_0, deltaTg);


f_0_newd = discrete_runge_kutta4(f_0_new, deltaT);
f_0_newdg = discrete_runge_kutta4(f_0_new, deltaTg);

%% ************************** Basis functions *****************************


% load('Design_data_450_20data.mat','X0','cent_f','cent','cent_c','Nrbf','Xcurrent'); %paper
load('Design_data_450_50data.mat','X0', 'cent_f','cent','cent_c','Nrbf','Xcurrent');   % paper
% load('Design_data_450_2000data.mat','X0', 'cent_f','cent','cent_c','Nrbf','Xcurrent');

upb=3;


Nlift_i = 10;

for i=1:N
    if Nlift_i==10

     liftFun{i} = @(xx)( [xx;xx(:,:).^2;xx(:,:).^3;  xx(1,:).*xx(2,:); (xx(1,:).^2).*xx(2,:);  xx(1,:).*(xx(2,:).^2);   ones(1, size(xx,2))] );
     dliftFun{i} = @(xx,yy)([yy; 2*xx(:,:).*yy(:,:); 3*(xx(:,:).^2).*yy(:,:); xx(1,:).*yy(2,:)+xx(2,:).*yy(1,:); xx(1,:).*xx(1,:).*yy(2,:)+2*xx(1,:).*xx(2,:).*yy(1,:);...
                             2*xx(1,:).*xx(2,:).*yy(2,:)+xx(2,:).*xx(2,:).*yy(1,:);   zeros(1,size(yy,2)) ]);

    else
     liftFun{i} = @(xx)( [xx;xx(1,:).^q] );
     dliftFun{i} = @(xx,yy)([yy; q*(xx(1,:).^(q-1)).*yy(1,:)]);
    end
end


%% ************************** Collect data ********************************
tic
disp('Starting data collection')
Nsim = 1;  
% Random initial conditions
if exist("Xcurrent") && N<4
    Ntraj=size(Xcurrent,2);
else
Ntraj = size(Xcurrent,2); % num of x(0)
Xcurrent= rand(n,Ntraj)*upb - upb/2;
end

X = []; X_f=[]; Y = [];  Y_f=[];  
X_f_g = []; Y_f_g=[];


Xcurrent_int=Xcurrent; % for transfer learning
%% one-step (Nsim=1) data generation
for i = 1:Nsim
% X_f= [x_1(0)..x_Ntraj(0), x_1(1)...x_Ntraj(1),...,x_1(Nsim-1)...x_Ntraj(Nsim-1)]

 
       X_current_temp = Xcurrent;
       for ii=1:deltaT/deltaTg
           X0next_temp=f_0dg(0, X_current_temp);
           X_current_temp = X0next_temp;
       end
    X_f = [X_f Xcurrent];
    Y_f= [Y_f X_current_temp];

     
     X0next= f_0(0,Xcurrent);   % for generator
     X_f_g = [X_f_g Xcurrent];
    Y_f_g = [Y_f_g X0next]; 


end

%% local-then-global
for i=1:N
X{i}=X_f(2*i-1:2*i,:);
Y{i}=Y_f(2*i-1:2*i,:);

X_g{i}= X_f_g(2*i-1:2*i,:);
Y_g{i}=Y_f_g(2*i-1:2*i,:);

tt=1:n;
tt2=tt(tt~=2*i-1);
tt3=tt2(tt2~=2*i);
ind_topo{i}=tt3;
U{i}=X_f(ind_topo{i},:);  
end


U_t{1}=U{1};  % Not Used
U_t{2}=X_f(input_ind,:);
U_t{3}=X_f(input_ind_transfer(1),:); 


U_t_g{1}=X_f_g(ind_topo{1},:); 
U_t_g{2}=X_f_g(input_ind,:);  
U_t_g{3}=X_f_g(input_ind_transfer(1),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Data collection DONE, time = %1.2f s \n', toc);

%% ******************************* Lift ***********************************

disp('Starting LIFTING')
tic

for i=1:N
Xlift{i} = liftFun{i}(X{i});
Ylift{i} = liftFun{i}(Y{i}); 


Xlift_g{i} = liftFun{i}(X_g{i});
Ylift_g{i} = dliftFun{i}(X_g{i},Y_g{i}); 


end


fprintf('Lifting DONE, time = %1.2f s \n', toc);

%% ********************** Build predictor *********************************

disp('Starting REGRESSION')
tic

%% m(gEDMD)

%%% (mEDMD)
for i=1:N
Wn = [Ylift{i}];
Vn_loop=Xlift{i};
if i>1
for j=1:size(U_t{i},1)  % nj
Vn_loop = [Vn_loop ; Xlift{i}*diag(U_t{i}(j,:))];
end
end
VVtn = Vn_loop*Vn_loop';
WVtn = Wn*Vn_loop';
Mn_t{i} = WVtn * pinv(VVtn);  % [A,B1,...Bn]
Alift_t{i}=Mn_t{i}(:, 1:Nlift_i);
    if i>1
     Blift_t{i}=Mn_t{i}(:, Nlift_i+1:end);
    else
     Blift_t{i}=[];   
     end
end

%%%   (mgEDMD)
for i=1:N
Wn = [Ylift_g{i}];
Vn_loop=Xlift_g{i};
if i>1
for j=1:size(U_t_g{i},1)   % nj
Vn_loop = [Vn_loop ; Xlift_g{i}*diag(U_t_g{i}(j,:))];
end
end
VVtn = Vn_loop*Vn_loop';
WVtn = Wn*Vn_loop';
Mn_t_g{i} = WVtn * pinv(VVtn);  % [A,B1,...Bn]
Alift_t_g{i}=Mn_t_g{i}(:, 1:Nlift_i);
    if i>1
     Blift_t_g{i}=Mn_t_g{i}(:, Nlift_i+1:end);
    else
     Blift_t_g{i}=[];   
     end
end

for ii=1:N
    if ii==1
       f_g_md{ii} = @(t,x)(Alift_t_g{ii}*x);
       
           f_g_md{ii}=discrete_runge_kutta4(f_g_md{ii}, deltaTg);
    
    else
       f_g_md{ii} = @(t,x,u)(Alift_t_g{ii}*x+Blift_t_g{ii}*kron(eye(size(U_t_g{ii},1)), x)*u);
       
           f_g_md{ii}=discrete_runge_kutta4_u(f_g_md{ii}, deltaTg);
       
    end
end



%% lEDMD

for i=1:N
Wn_lin = [Ylift{i}];
Vn_lin=Xlift{i};
if i>1
Vn_lin = [Xlift{i} ; Xlift{1}];
end
VVtn_lin = Vn_lin*Vn_lin';
WVtn_lin = Wn_lin*Vn_lin';
Mn_lin_t{i} = WVtn_lin * pinv(VVtn_lin);  % [A,B1,...Bn]
Alift_lin_t{i}=Mn_lin_t{i}(:, 1:Nlift_i);
    if i>1
     Blift_lin_t{i}=Mn_lin_t{i}(:, Nlift_i+1:end);
    else
     Blift_lin_t{i}=[];   
     end
end



fprintf('Regression done, time = %1.2f s \n', toc);


%% **********************transfer learning *****************************
if 1
close all

disp('Starting transfer learning')
tic



Tmax = 0.5;
Nsim = Tmax/deltaT;
Nsim_g = Tmax/deltaTg;

Nrun = 1; % Nsim;

Nsample=500;

X0=rand(n/N+n, Nsample)*1-0.5; % add #4

EMAX_f=[];
EMAX =[];
EMAX_t=[];
EMAX_t_g = [];
EMAX_c=[];
EMAX_lin=[];
EMAX_lin_t=[];


% Lifted initial condition
% copy liftfct of #2 for #4
liftFun{4} = liftFun{2};
% copy koopman model of #2 for #4
Alift_t{4}=Alift_t{2};
Blift_t{4}=Blift_t{2};
Alift_t_g{4}=Alift_t_g{2};
Blift_t_g{4}=Blift_t_g{2};


Alift_lin_t{4}=Alift_lin_t{2};
Blift_lin_t{4}=Blift_lin_t{2};

if ~isempty(id_input_3) 

Xcurrent_4=rand(n/N,Ntraj)*3-3/2;
Xcurrent=[Xcurrent_int;Xcurrent_4];
X_f=zeros(n+n/N, Ntraj*Nrun);
Y_f=zeros(n+n/N, Ntraj*Nrun);
X_f_g=zeros(n+n/N, Ntraj);
Y_f_g=zeros(n+n/N, Ntraj);
for i = 1:Nrun
% X_f= [x_1(0)..x_Ntraj(0), x_1(1)...x_Ntraj(1),...,x_1(Nsim-1)...x_Ntraj(Nsim-1)]
    X0next= f_0_newd(0,Xcurrent); 
    X_f(:, Ntraj*(i-1)+1:i*Ntraj) = Xcurrent;
    Y_f(:, Ntraj*(i-1)+1:i*Ntraj) = X0next;     
    if Nrun>1
    Xcurrent = X0next;
    end  
 
     if  Nrun==1
     X0next= f_0_new(0,Xcurrent);   % for generator
    X_f_g(:, Ntraj*(i-1)+1:i*Ntraj) =  Xcurrent;
    Y_f_g(:, Ntraj*(i-1)+1:i*Ntraj) =  X0next; 
     end

end
for i=1:N+1
X{i}=X_f(2*i-1:2*i,:);
Y{i}=Y_f(2*i-1:2*i,:);

X_g{i}=X_f_g(2*i-1:2*i,:);
Y_g{i}=Y_f_g(2*i-1:2*i,:);
end
U_t{3}=X_f(input_ind_transfer,:); 
U_t_g{3}=X_f_g(input_ind_transfer,:);
Xlift_3=liftFun{3}(X{3});
Ylift_3=liftFun{3}(Y{3});
Xlift_3_g = liftFun{3}(X_g{3});
Ylift_3_g = dliftFun{3}(X_g{3}, Y_g{3});
Xlift_1=liftFun{1}(X{1});
Xlift_4=liftFun{4}(X{4});
Xlift_1_g=liftFun{1}(X_g{1});
Xlift_4_g=liftFun{4}(X_g{4});

% mEDMD

Wn=Ylift_3-Alift_t{3}*Xlift_3-Blift_t{3}(:,1:Nlift_i)*Xlift_3.*U_t{3}(1,:);
Vn=[Xlift_3.*U_t{3}(2,:)];
VVtn=Vn*Vn';
WVtn=Wn*Vn';
Blift_3_p=WVtn*pinv(VVtn);
Blift_t{3}=[Blift_t{3} Blift_3_p];

%mgEDMD
Wn_g=Ylift_3_g-Alift_t_g{3}*Xlift_3_g-Blift_t_g{3}(:,1:Nlift_i)*Xlift_3_g.*U_t_g{3}(1,:);
Vn_g=[Xlift_3_g.*U_t_g{3}(2,:)];
VVtn_g=Vn_g*Vn_g';
WVtn_g=Wn_g*Vn_g';
Blift_3_p_g=WVtn_g*pinv(VVtn_g);
Blift_t_g{3}=[Blift_t_g{3} Blift_3_p_g];


f_g_md{3} = @(t,x,u)(Alift_t_g{3}*x+Blift_t_g{3}*kron(eye(size(U_t_g{3},1)), x)*u);
f_g_md{3}=discrete_runge_kutta4_u(f_g_md{3}, deltaTg);

 

% lEDMD
Wn=Ylift_3-Alift_lin_t{3}*Xlift_3-Blift_lin_t{3}*Xlift_1;
Vn=liftFun{id_input_3}(X{id_input_3});
VVtn=Vn*Vn';
WVtn=Wn*Vn';
Blift_lin_3_p=WVtn*pinv(VVtn);
Blift_lin_t{3}=[Blift_lin_t{3} Blift_lin_3_p];


end

f_g_md{4}=f_g_md{2};




fprintf('learning finished, time = %1.2f s \n', toc);



%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%

disp('----------------------------------------------------');
disp('----------------------------------------------------');
disp('Starting simulation');


% % % evaluate assumption in Prop.4.7 with 2000 Samples 
%  i=2,j=3 ;  i=3,j=2
if N_n >N && ~isempty(id_input_3) && Ntraj>1500
  boundx1= 0.5*2;
  boundx2 = 0.5*2;
  nz= size(Alift_t_g{3},1);
  z0_max_i= [max(sum(abs(liftFun{3}(X0(2*3-1:2*3,:)))));  max(sum(abs(liftFun{4}(X0(2*4-1:2*4,:)))))   ]; 
  x_max_j = [(boundx1/2+boundx2/2) ; (boundx1/2+boundx2/2)  ];
  Lv0_i = [norm(Alift_t_g{3},1) ; norm(Alift_t_g{4},1) ];
  sum_Lij = [ norm(Blift_t_g{3}(:,nz+1:end),1); norm(Blift_t_g{4}(:,1:nz),1) ];
  sum_Lr_ij = [ norm(Blift_t_g{3}(:,nz+1:end)+Alift_t_g{3},1); norm(Blift_t_g{4}(:,1:nz)+Alift_t_g{4},1) ];
  vi = [2*Lv0_i(1)+2*x_max_j(1)*(Lv0_i(1)+sum_Lr_ij(1));  2*Lv0_i(2)+2*x_max_j(2)*(Lv0_i(2)+sum_Lr_ij(2))  ];
 E34 = z0_max_i(1)*Tmax*exp(vi(1)*Tmax)*sum_Lij(1);
 E43 = z0_max_i(2)*Tmax*exp(vi(2)*Tmax)*sum_Lij(2);
 disp('-----------------------------------------------------');
 disp(['E34=', num2str(E34)]);
 disp(['E43=', num2str(E43)]);
 disp(['E34*E43=', num2str(E34*E43)]);
disp('------------------------------------------------------');
end

for Ni=1:Nsample
x0 = X0(:,Ni);
x_true = x0(1:N_o*2,:);
x_true_g = x0(1:N_o*2,:);


for i=1:N_o
xlift{i}= liftFun{i}(x0(2*i-1:2*i));  % 
xlift_t{i}=xlift{i};  % mEDMD
xlift_t_g{i} = xlift{i};  %mgEDMD
xlift_lin_t{i}=xlift{i};  % lEDMD
end


%% sim for lEDMD and mEDMD
for i = 0:Nsim-1



x_n_old_t{2}=xlift_t{1}(input_ind,end);  
if  ~isempty(id_input_3)
x_n_old_t{3}=[xlift_t{1}(input_ind_transfer(1), end); xlift_t{id_input_3}(input_ind_new_trans,end)];
else
 x_n_old_t{3}=[x_n_old_t{2}];   
end
x_n_old_t{4}=xlift_t{3}(input_ind_2,end);

x_lift_old_t{2}=[xlift_lin_t{1}(:,end)]; 
if ~isempty(id_input_3)
x_lift_old_t{3}=[x_lift_old_t{2}; xlift_lin_t{id_input_3}(:,end)];
else
 x_lift_old_t{3}=[x_lift_old_t{2}]; 
end
x_lift_old_t{4}=xlift_lin_t{3}(:,end);


    for ii=1:N_o  % add agent #4

    if ii==1
    xlift_t{ii} = [xlift_t{ii}, Alift_t{ii}*xlift_t{ii}(:,end)];  % mEDMD
   
    xlift_lin_t{ii}=[xlift_lin_t{ii}, Alift_lin_t{ii}*xlift_lin_t{ii}(:,end)];  % lEDMD
    else
     xlift_t{ii} = [xlift_t{ii}, Alift_t{ii}*xlift_t{ii}(:,end) + Blift_t{ii}*kron(eye(length(x_n_old_t{ii})), xlift_t{ii}(:,end))*x_n_old_t{ii}];  % mEDMD

    xlift_lin_t{ii}=[xlift_lin_t{ii}, Alift_lin_t{ii}*xlift_lin_t{ii}(:,end)+Blift_lin_t{ii}*x_lift_old_t{ii}];   % lEDMD
    end
    end


    % True dynamics
   
     x_true = [x_true, f_0_newd(0,x_true(:,end)) ];   

    
end




%% sim for mgEDMD 
for i=1:Nsim_g

x_n_old_t_g{2}=xlift_t_g{1}(input_ind,end);  
if  ~isempty(id_input_3)
x_n_old_t_g{3}=[xlift_t_g{1}(input_ind_transfer(1),end); xlift_t_g{id_input_3}(input_ind_new_trans,end)];
else
x_n_old_t_g{3}=[xlift_t_g{1}(input_ind_transfer(1),end)]; 
end
x_n_old_t_g{4}=xlift_t_g{3}(input_ind_2,end);

for ii=1:N_o
if ii==1
 %%% mgEDMD
     xlift_t_g{ii} = [xlift_t_g{ii}, f_g_md{ii}(0,xlift_t_g{ii}(:,end))];
else

     %%% mgEDMD
     xlift_t_g{ii} = [xlift_t_g{ii}, f_g_md{ii}(0,xlift_t_g{ii}(:,end), x_n_old_t_g{ii})];
end
end
    % True dynamics
  
       x_true_g = [x_true_g, f_0_newdg(0,x_true_g(:,end)) ];
   
end


%% use x_true_g as true traj.

x_true = x_true_g (:, 1:deltaT/deltaTg:Nsim_g+1);



for ii=1:N_o

x_koopn_t{ii}=xlift_t{ii}(1:n/N,:);   % mEDMD
x_koopn_t_g{ii}=xlift_t_g{ii}(1:n/N,:);  %mgEDMD
x_koopn_lin_t{ii}=xlift_lin_t{ii}(1:n/N,:);  %lEDMD

end

emax_t = max(abs([x_koopn_t{1};x_koopn_t{2};x_koopn_t{3};x_koopn_t{4}]-x_true),[],2);
emax_t_g = max(abs([x_koopn_t_g{1};x_koopn_t_g{2};x_koopn_t_g{3};x_koopn_t_g{4}]-x_true_g),[],2);
emax_lin_t = max(abs([x_koopn_lin_t{1};x_koopn_lin_t{2};x_koopn_lin_t{3};x_koopn_lin_t{4}]-x_true),[],2);

EMAX_t = [EMAX_t emax_t];
EMAX_t_g = [EMAX_t_g emax_t_g];
EMAX_lin_t = [EMAX_lin_t emax_lin_t];

end

     file_name= ['Datensatz\EMAX_all_', num2str(Ntraj),'_inv_transfer.mat'];
    save(file_name, 'EMAX_t', 'EMAX_lin_t', 'EMAX_t_g');
   

%% overall plot 


psize = [14.65, 9.6];
ppos  = [0.3, 0.2, psize(1)-0.3 psize(2)-0.2];

gcaPos = [0.1 0.1 0.85 0.75];

FTsize_axisnumbers = 14;
Legend_size = 18;
Title_size = 18;
Temp_comb=[ EMAX_lin_t;EMAX_t; EMAX_t_g];
Nsample=size(EMAX_t,2);
if ~isempty(id_input_3) && N_n==N
     N=3; % plot only the first 3 agents
    col=repmat({'g','c','w'}, 1,3);
    Temp_comb=[ EMAX_lin_t(1:2*N,:);EMAX_t(1:2*N,:); EMAX_t_g(1:2*N,:)];
else
    N=4; % plot all 4 agents
    col=repmat({'y','g','c','w'}, 1,3);
end

EMAX_all=zeros(N*3, size(Temp_comb,2));
Tmax=0.5;
for i=1:N*3
%     EMAX_all(i,:)=log(sqrt(Temp_comb(2*i-1,:).^2+Temp_comb(2*i,:).^2));
     EMAX_all(i,:)=log( abs(Temp_comb(2*i-1,:))+abs(Temp_comb(2*i,:)) );
end
g = kron([1:N*3]', ones(Nsample,1));
if N==4
pos=[2 2.4 2.8 3.2  3.8+0.2 4.2+0.2 4.6+0.2 5+0.2  3.8+2.2 4.2+2.2 4.6+2.2 5+2.2     ];
posmean=[ mean(pos(1:4)) mean(pos(5:8))  mean(pos(9:end)) ];
pos2 = mean([pos(4:5);pos(8:9)],2);
else
    pos=[2 2.6 3.2   4  4.6 5.2  4+2  4.6+2  5.2+2     ];
    posmean=[ mean(pos(1:3)) mean(pos(4:6))  mean(pos(7:end)) ];
    pos2 = mean([pos(3:4);pos(6:7)],2);
end

% col=repmat({cm(6,:),'g','c','w'}, 1,2);

row1 = { 'lEDMD', 'mEDMD', 'mgEDMD' };
tickLabels = row1;
figure(11)
boxplot([reshape(EMAX_all', [numel(EMAX_all),1])],g, 'positions',pos );
% set(gca,'xtick',posmean,'YLim', [-1/21 1]*0.4);
set(gca,'xtick',posmean);  %, ,'YLim', [-12 22]
set(gca,'xtick',posmean,'YLim', [-18 -3]); 
set(gca,'xticklabel',tickLabels, 'TickLabelInterpreter', 'latex');  
set(gca,'YScale','linear');
grid on
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),col{j},'FaceAlpha',.6);
end
c = get(gca, 'Children');
xline(pos2,':');
set(gcf,'PaperUnits','centimeters','PaperSize',psize,'PaperPosition',ppos)
set(gca,'Position',gcaPos,'fontsize',FTsize_axisnumbers);
 title(['$\ln(\|\Delta x_i\|_{\infty})$, ', num2str(Ntraj), ' snapshots'],'interpreter','latex', 'FontSize', Title_size); 
 if N==4
hleg1 = legend(c(1:4), ' $ x_1 $', '$x_2$', '$x^{\sharp}_3$', '$x_4$', 'Interpreter', 'latex','FontSize', Legend_size, 'Orientation','horizontal'); %
 else
   hleg1 = legend(c(1:3), ' $ x_1 $', '$x_2$', '$x^{\sharp}_3$', 'Interpreter', 'latex','FontSize', Legend_size,'Orientation','horizontal');
 end
set(gca, 'XGrid','off', 'YGrid','on');

if printing
    if ~isempty(id_input_3)
        if N_n>3
    name = sprintf(['Overall_Duffing_', num2str(Ntraj),'_transfer_new_part_ident']);
        else
          name = sprintf(['Overall_Duffing_', num2str(Ntraj),'_transfer_new_part']);
        end
    else
     name = sprintf(['Overall_Duffing_', num2str(Ntraj),'_transfer_new_ident']);  
    end
    print([folder name],'-depsc','-painters')
    print([folder name],'-dpdf','-painters')
    print([folder name],'-dbmp')
    savefig([folder name]) 

end


end



