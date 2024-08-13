clear all
% close all
% for example of fig. 4 from section 5.2
% This code is adapted from  codes used in the paper M. Korda, M. Putinar, I. Mezic (2020)
addpath('./Resources');
addpath('Datensatz\');
addpath('..\..\Toolbox_share\');
rng(2141444)

%% preprocessing(ensure X0 cent_f cent Nrbf are same)

debug_x0 = 1; % 1: use default x0, otherwise generate new x0
input_ind =2;  % index of agent 2 (3) for agent 3 (2)
printing =1;
folder = 'Figs/';

%% *************************** Dynamics ***********************************
n = 6; %num of full states\dot
N=3;   % num of agents
the30=deg2rad(32);
the3h=deg2rad(-29);
the9h=deg2rad(50);
the12o=deg2rad(-13);
the12h=deg2rad(10);
c123=0.001;
c312=0.1;
c39=0.001;
c129=c39;
ep3=0.01;
ep9=0.1;
ep12=0.01;
w=1;
p3=(4*c312*(the3h-the12h)+4*the3h*ep3+the9h^2*c39)/(the3h^3*ep3);
Om3=sqrt(4*w);
p9=4/(the9h^2);
Om9=w;
p12=(4*c123*(the12h-the3h)+4*the12h*ep12+the9h^2*c129)/(the12h^3*ep12);
Om12=Om3;
%  (#1, #3)->#2, (#1, #2)->#3
        f_0 =  @(t,x)([ x(2,:) ; 0.1*(1-p9*(x(1,:)).^2.*x(2,:))-Om9^2*x(1,:); ...
        x(4,:); 0.01*(1-p3*(x(3,:)).^2.*x(4,:))-Om3^2*(x(3,:))+0.001*x(1,:).*x(2,:)+0.1*(x(4,:)-x(6,:));...
        x(6,:); 0.01*(1-p12*x(5,:).^2.*x(6,:))-Om12^2*x(5,:)+0.001*x(1,:).*x(2,:)+0.1*(x(6,:)-x(4,:))]);  %%(Max S. Dutra, Modeling of a bipedal.... (2022))




%% ************************** Discretization ******************************

deltaT = 0.01;
deltaTg = 0.005;
%Runge-Kutta 4
f_0d = discrete_runge_kutta4(f_0, deltaT);
f_0dg = discrete_runge_kutta4(f_0, deltaTg);



%% ************************** Basis functions *****************************

basisFunction = 'rbf';  
 % RBF centers
rbf_type = 'thinplate';


% load('Design_data_450_2500data_Pol_new.mat'); % paper
load('Design_data_450_8000data_Pol_new.mat'); %paper


boundx1=pi;
boundx2=2;

if exist("Nrbf")
else Nrbf = 450;
end
Nrbf_i=Nrbf/N;

if exist("cent")
else
    for i=1:N
     cent{i} = diag([boundx1, boundx2])*rand(n/N,Nrbf_i) -[boundx1/2; boundx2/2];   

    end
    cent_f=kron(eye(N), diag([boundx1, boundx2]))*rand(n,Nrbf) - kron(ones(N,1), [boundx1/2; boundx2/2]) ;
end

for i=1:N

liftFun{i} = @(xx)( [xx;rbf(xx,cent{i},rbf_type)] );  
dliftFun{i} = @(xx,yy)([yy; drbf(xx,cent{i},yy,rbf_type)]);
liftFun_lin{i} = liftFun{i}; 

if i==1  
cent_1=cent{i};
cent_1 = cent_1(:,1:end-1);
liftFun{i} = @(xx)( [xx; xx(1,:).*xx(2,:); rbf(xx,cent_1,rbf_type)] ); 
end

end


% Lifting mapping - RBFs + the state itself
liftFun_f = @(xx)( [xx;rbf(xx,cent_f,rbf_type)] );  
Nlift_f = Nrbf + n;
Nlift_i = Nlift_f/N;  




%% ************************** Collect data ********************************
tic
disp('Starting data collection')
Nsim = 1;  %200  length of traj 

if exist("Xcurrent_int")
    Ntraj=size(Xcurrent_int,2);
else
Ntraj = 2500; % num of x(0)
end

% Random initial conditions
if exist('Xcurrent_int')
else
Xcurrent_int = kron(eye(N), diag([boundx1, boundx2]))*rand(n,Ntraj) - kron(ones(N,1), [boundx1/2; boundx2/2]) ;
end

X = []; X_f=[]; X_f_e=[]; Y = [];  Y_f=[]; Y_f_e=[];  
X_f_g = []; Y_f_g =[];



%% EDMD
Xcurrent=Xcurrent_int;
for i = 1:Nsim
% X_f= [x_1(0)..x_Ntraj(0), x_1(1)...x_Ntraj(1),...,x_1(Nsim-1)...x_Ntraj(Nsim-1)]

   
        X_current_temp = Xcurrent;
        for ii=1:deltaT/deltaTg
           X0next_temp=f_0dg(0, X_current_temp);
           X_current_temp = X0next_temp; 
        end
    Xcurrent_1=Xcurrent(1,:).*Xcurrent(2,:);
    X_f = [X_f Xcurrent];
    X_f_e= [X_f_e Xcurrent_1 ];
    Y_f= [Y_f X0next_temp]; 

     X0next= f_0dg(0,Xcurrent);   % for generator approx
     X_f_g = [X_f_g Xcurrent];
   
      Y_f_g = [Y_f_g X0next]; 
    
end

%% m(g)EDMD lEDMD
for i=1:N
 
X{i}=X_f(2*i-1:2*i,:);
Y{i}=Y_f(2*i-1:2*i,:);

X_g{i}=X_f_g(2*i-1:2*i,:);
Y_g{i}=Y_f_g(2*i-1:2*i,:);
    
end


U_t{1}=X_f_g(1:2,:);  % Not Used

U_t{2}=[X_f_e; X{3}(input_ind,:)]; 
U_t{3}=[X_f_e; X{2}(input_ind,:)];




U_t_g{1}=X_f_g(1:2,:);   %not used

U_t_g{2}=[X_f_e; X_g{3}(input_ind,:)];
U_t_g{3}=[X_f_e; X_g{2}(input_ind,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ******************************* Lift ***********************************

disp('Starting LIFTING')
tic
Xlift_f = liftFun_f(X_f);
Ylift_f = liftFun_f(Y_f);

for i=1:N
Xlift{i} = liftFun{i}(X{i});
Ylift{i} = liftFun{i}(Y{i});  
Xlift_lin{i} = liftFun_lin{i}(X{i});  
Ylift_lin{i} = liftFun_lin{i}(Y{i});  

Xlift_g{i} = liftFun{i}(X_g{i});
Ylift_g{i} = liftFun{i}(Y_g{i}); 
end
fprintf('Lifting DONE, time = %1.2f s \n', toc);

%% ********************** Build predictor *********************************

disp('Starting REGRESSION')
tic
%% EDMD


W = [Ylift_f];
V = [Xlift_f];
VVt = V*V';
WVt = W*V';
Alift_f = WVt * pinv(VVt); 
 


%% mEDMD

for i=1:N
Wn = [Ylift{i}];
Vn_loop=Xlift{i};
if i>1
    if ~isempty(U_t{i})
for j=1:size(U_t{i},1)  % ni
Vn_loop = [Vn_loop ; Xlift{i}*diag(U_t{i}(j,:))];
end
    end
end
VVtn = Vn_loop*Vn_loop';
WVtn = Wn*Vn_loop';
Mn_t{i} = WVtn * pinv(VVtn);  % [A,B1,...Bn]
Alift_t{i}=Mn_t{i}(:, 1:Nlift_i);
    if i>1 && ~isempty(U_t{i})
     Blift_t{i}=Mn_t{i}(:, Nlift_i+1:end);
    else
     Blift_t{i}=[];   
    end
end

%% mgEDMD
for i=1:N
Wn = [Ylift_g{i}];
Vn_loop=Xlift_g{i};
if i>1
     if ~isempty(U_t_g{i})
for j=1:size(U_t_g{i},1)  % nj
Vn_loop = [Vn_loop ; Xlift_g{i}*diag(U_t_g{i}(j,:))];
end
     end
end
VVtn = Vn_loop*Vn_loop';
WVtn = Wn*Vn_loop';
Mn_t_g{i} = WVtn * pinv(VVtn);  % [A,B1,...Bn]
Alift_t_g{i}=Mn_t_g{i}(:, 1:Nlift_i);
    if i>1 && ~isempty(U_t_g{i})
     Blift_t_g{i}=Mn_t_g{i}(:, Nlift_i+1:end);
    else
     Blift_t_g{i}=[];   
     end
end

for ii=1:N
    if ii==1|| isempty(U_t_g{ii})
       f_g_md{ii} = @(t,x)(Alift_t_g{ii}*x);
    else
       f_g_md{ii} = @(t,x,u)(Alift_t_g{ii}*x+Blift_t_g{ii}*kron(eye(size(U_t_g{ii},1)), x)*u);
    end
end



%% lEDMD

for i=1:N
Wn_lin = [Ylift_lin{i}];
Vn_lin=Xlift_lin{i};
if i>1
ind=1:N;
ind(i)=[];

Vn_lin = [Xlift_lin{i} ; Xlift_lin{ind(1)}; Xlift_lin{ind(2)}];

end
VVtn_lin = Vn_lin*Vn_lin';
WVtn_lin = Wn_lin*Vn_lin';
Mn_lin_t{i} = WVtn_lin * pinv(VVtn_lin);  % [A,B1,...Bn]
Alift_lin_t{i}=Mn_lin_t{i}(:, 1:Nlift_i);
    if i>1 && ~isempty(U_t_g{i})
     Blift_lin_t{i}=Mn_lin_t{i}(:, Nlift_i+1:end);
    else
     Blift_lin_t{i}=[];   
     end
end


fprintf('Regression done, time = %1.2f s \n', toc);


%% ********************** error comparision (boxplot)
% EDMD vs lEDMD vs mEDMD vs mgEDMD
if 1
close all
clf
Tmax = 0.5;
Nsim = Tmax/deltaT;
Nsim_g = Tmax/deltaTg;
if exist("X0") && debug_x0
    Nsample=size(X0,2);
else
Nsample=500;
scale=0.3;
X0 = kron(eye(N), diag([scale*boundx1, scale*boundx2]))*rand(n,Nsample) - kron(ones(N,1), scale*[boundx1/2; boundx2/2]) ;
end

EMAX_f=[];
EMAX =[];
EMAX_t=[];
EMAX_t_g=[];
EMAX_c=[];
EMAX_lin=[];
EMAX_lin_t=[];
for Ni=1:Nsample
x0 = 0.2*X0(:,Ni);
x_true = x0;
x_true_g = x0;
% Lifted initial condition

for i=1:N
xlift{i}= liftFun{i}(x0(2*i-1:2*i)); 
xlift_t{i}=xlift{i};  % mEDMD
xlift_t_g{i}=xlift{i};  % mgEDMD

xlift_lin_t{i}=liftFun_lin{i}(x0(2*i-1:2*i));  % lEDMD
end
xlift_f =  liftFun_f(x0);  % EDMD



%%
% Simulate
x_n_old=[];
x_n_old_t=[];
x_lift_old = [];
x_lift_old_t=[];
x_n_old_t_g = [];
for i = 0:Nsim-1
  

%% EDMD
    xlift_f=[xlift_f, Alift_f*xlift_f(:,end)];
%% mEDMD ,  lEDMD


      x_n_old_t{2}=[xlift_t{1}(1,end).*xlift_t{1}(2,end); xlift_t{3}(input_ind,end)];  % mEDMD 
       x_n_old_t{3}=[xlift_t{1}(1,end).*xlift_t{1}(2,end); xlift_t{2}(input_ind,end)]; 
      x_lift_old_t{2}=[xlift_lin_t{1}(:,end); xlift_lin_t{3}(:,end)];       % lEDMD  
      x_lift_old_t{3}=[xlift_lin_t{1}(:,end); xlift_lin_t{2}(:,end)];

   for ii=1:N

    if ii< 2 
    xlift_t{ii} = [xlift_t{ii}, Alift_t{ii}*xlift_t{ii}(:,end)];  % mEDMD
   
    xlift_lin_t{ii}=[xlift_lin_t{ii}, Alift_lin_t{ii}*xlift_lin_t{ii}(:,end)];  % lEDMD
    else   
  
     xlift_t{ii} = [xlift_t{ii}, Alift_t{ii}*xlift_t{ii}(:,end) + Blift_t{ii}*kron(eye(size(x_n_old_t{ii},1)), xlift_t{ii}(:,end))*x_n_old_t{ii}];  % mEDMD
  
    xlift_lin_t{ii}=[xlift_lin_t{ii}, Alift_lin_t{ii}*xlift_lin_t{ii}(:,end)+Blift_lin_t{ii}*x_lift_old_t{ii}];   % lEDMD 
    end

    end

    % True dynamics
    x_true = [x_true, f_0d(0,x_true(:,end)) ];

    
end

%% mgEDMD
for i=1:Nsim_g
  
   x_n_old_t_g{2}=[xlift_t_g{1}(1,end).*xlift_t_g{1}(2,end); xlift_t_g{3}(input_ind,end)];  
   x_n_old_t_g{3}=[xlift_t_g{1}(1,end).*xlift_t_g{1}(2,end); xlift_t_g{2}(input_ind,end)];


for ii=1:N
if ii<2  

     xlift_t_g{ii} = [xlift_t_g{ii}, f_g_md{ii}(0,xlift_t_g{ii}(:,end))];
else

     xlift_t_g{ii} = [xlift_t_g{ii}, f_g_md{ii}(0,xlift_t_g{ii}(:,end), x_n_old_t_g{ii})];
end
end
    % True dynamics
    x_true_g = [x_true_g, f_0dg(0,x_true_g(:,end)) ];

end




%%


x_true = x_true_g (:, 1:deltaT/deltaTg:Nsim_g+1);


vec=1:n;
x_koop_f = xlift_f(vec,:); %  EDMD


for ii=1:N
x_koopn_t{ii}=xlift_t{ii}(1:n/N,:);  % mEDMD
x_koopn_lin_t{ii}=xlift_lin_t{ii}(1:n/N,:);  % lEDMD
x_koopn_t_g{ii}=xlift_t_g{ii}(1:n/N,:);  % mgEDMD
end

emax_f = max(abs(x_koop_f-x_true),[],2);
emax_t = max(abs([x_koopn_t{1};x_koopn_t{2};x_koopn_t{3}]-x_true),[],2);
emax_t_g = max(abs([x_koopn_t_g{1};x_koopn_t_g{2};x_koopn_t_g{3}]-x_true_g),[],2);
emax_lin_t = max(abs([x_koopn_lin_t{1};x_koopn_lin_t{2};x_koopn_lin_t{3}]-x_true),[],2);


EMAX_f = [EMAX_f  emax_f];
EMAX_t = [EMAX_t emax_t];
EMAX_t_g = [EMAX_t_g emax_t_g];
EMAX_lin_t = [EMAX_lin_t emax_lin_t];

end



    file_name= ['Datensatz\EMAX_all_', num2str(Ntraj),'_', basisFunction , '_Pol.mat'];

    save(file_name,'EMAX_t', 'EMAX_f', 'EMAX_lin_t', 'EMAX_t_g');
   

end
%%  comparison plots 
% close all
psize = [14.65, 9.6];
ppos  = [0.3, 0.2, psize(1)-0.3 psize(2)-0.2];

gcaPos = [0.1 0.1 0.85 0.75];

FTsize_axisnumbers = 14;
Legend_size = 18;
Title_size = 18;
load(file_name);

Temp_comb=[EMAX_f;   EMAX_lin_t; EMAX_t; EMAX_t_g];
Nsample=size(EMAX_f,2);
ni=2;
N=3;
row1 = {'EDMD',  'lEDMD','mEDMD', 'mgEDMD'};
nr= length(row1);
EMAX_all=zeros(N*nr, size(Temp_comb,2));
Tmax=0.5;
for i=1:N*nr
%     EMAX_all(i,:)=log(sqrt(Temp_comb(2*i-1,:).^2+Temp_comb(2*i,:).^2));
     EMAX_all(i,:)=log( abs(Temp_comb(2*i-1,:))+abs(Temp_comb(2*i,:)) );
end
g = kron([1:N*nr]', ones(Nsample,1));
pos=[1 1.4 1.8  2.4 2.8 3.2  3.8 4.2 4.6  5.2 5.6 6   ];
pos=pos(1:3*nr);
pos2 = mean([pos(3:4);pos(6:7); pos(9:10)],2);
col=repmat({'g','c','w'}, 1,nr);
posmean=pos([2:3:3*nr-1 ]);
tickLabels = row1;
figure(11)
boxplot([reshape(EMAX_all', [numel(EMAX_all),1])],g, 'positions',pos );
set(gca,'xtick',posmean,'YLim', [-12, 0]);  %, 'YLim', [-12, 8]
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
hleg1 = legend(c(1:3), ' $ x_1 $', '$x_2$', '$x_3$', 'Interpreter', 'latex','FontSize', Legend_size,'Orientation','horizontal');

set(gca, 'XGrid','off', 'YGrid','on');

if printing
    name = sprintf(['Overall_VanderPol_', num2str(Ntraj),'_new_rbf']);
    print([folder name],'-depsc','-painters')
    print([folder name],'-dpdf','-painters')
    print([folder name],'-dbmp')
    savefig([folder name]) 

end


