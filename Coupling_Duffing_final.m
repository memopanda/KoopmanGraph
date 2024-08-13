clear all
close all
% for example of fig. 3 from section 5.2
% This code is adapted from codes used in the paper:M. Korda, M. Putinar, I. Mezic (2020)
addpath('Resources\')
addpath('Datensatz\');
addpath('..\..\Toolbox_share\');
% rng shuffle
rng(2141444)
%% preprocessing(ensure X0 cent_f cent cent_c Nrbf are same)

input_ind = 1; % index of agent 1 for input of agent 2 and 3
printing=1;
folder='Figs/';
%% *************************** Dynamics ***********************************
n = 6; %num of full states\dot
N=3;   % num of agents
nj=2; % number of inputs for each subsysetm

theta=0.5;
beta1=0.5;
beta2=1;
beta3=-0.5;
ga2=0.5;
ga3=1;
ga4=0;
f_0 =  @(t,x)([ 0.5*x(2,:); -theta*x(2,:)-2*(beta1)*x(1,:).^3;...
    0.5*x(4,:); -theta*x(4,:)-2*(beta3)*x(3,:).^3+0.5*ga2*x(1,:); ...
    0.5*x(6,:); -theta*x(6,:)-2*(beta3)*x(5,:).^3+0.5*ga3*x(1,:)] );  %%(Corbinain Schlosser, Sparsity Structures.... (2022))



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

load('Design_data_450_1500data.mat','X0', 'cent_f','cent','cent_c','Nrbf','Xcurrent'); %
% load('Design_data_450_5000data.mat','X0', 'cent_f','cent','cent_c','Nrbf','Xcurrent');


upb=3;
switch basisFunction
    case 'rbf'
        if exist("Nrbf")
        else Nrbf = 450;
        end
        Nrbf_i=Nrbf/N;

        if exist("cent")
        else
            for i=1:N
                cent{i} = rand(n/N,Nrbf_i)*ubp - upb/2;
            end
            cent_f=rand(n,Nrbf)*upb - upb/2;
        end

        for i=1:N
            liftFun{i} = @(xx)( [xx;rbf(xx,cent{i},rbf_type)] );

            dliftFun{i} = @(xx,yy)([yy; drbf(xx,cent{i},yy,rbf_type)]);
        end


        % For C. Schlosser
        if exist("cent_c")
        else
            cent_c{1}=cent{1};
            cent_c{2}=rand(4,Nrbf_i)*upb - upb/2;
            cent_c{3}=rand(4,Nrbf_i)*upb - upb/2;
        end
        %

        liftFun_c{1} = liftFun{1};
        liftFun_c{2} = @(xx)( [[zeros(2) eye(2)]*xx;rbf(xx,cent_c{2},rbf_type)] );
        liftFun_c{3} = @(xx)( [[zeros(2) eye(2)]*xx;rbf(xx,cent_c{3},rbf_type)] );


        % Lifting mapping - RBFs + the state itself
        liftFun_f = @(xx)( [xx;rbf(xx(1:2,:),cent{1},rbf_type);rbf(xx(3:4,:),cent{2},rbf_type);rbf(xx(5:6,:),cent{3},rbf_type)] );
        Nlift_f = Nrbf + n;
        Nlift_i = Nrbf_i + n/N;

        %%%%% monomial observables
    case 'mono'
        for i=1:N
            liftFun{i} = @(xx)( [xx;mono_fct(xx,2,3)] );
            dliftFun{i} = @(xx,yy)([yy; mono_fct(xx,2,3,yy) ]);



        end
        liftFun_c{1} = liftFun{1};
        liftFun_c{2} = @(xx)( [[zeros(2) eye(2)]*xx;mono_fct(xx(1:2,:),2,3); mono_fct(xx(3:4,:),2,3)] );
        liftFun_c{3} = @(xx)( [[zeros(2) eye(2)]*xx;mono_fct(xx(1:2,:),2,3); mono_fct(xx(3:4,:),2,3)] );

        liftFun_f = @(xx)( [xx;mono_fct(xx(1:2,:),2,3); mono_fct(xx(3:4,:),2,3); mono_fct(xx(5:6,:),2,3)] );
        Nrbf_i = size(mono_order_gen(2,3,0),1);
        Nlift_f = Nrbf_i*N;
        Nlift_i = Nrbf_i;

end


%% ************************** Collect data ********************************
tic
disp('Starting data collection')
Nsim = 1;  %200  length of traj
% Random initial conditions
if exist("Xcurrent")
    Ntraj=size(Xcurrent,2);
else
    Ntraj = 1500; % num of x(0)
    Xcurrent= rand(n,Ntraj)*upb - upb/2;
end

X = []; X_f=[]; Y = [];  Y_f=[];
X_f_g = []; Y_f_g=[];


%% one-step (Nsim=1) data generation
for i = 1:Nsim

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

%% mEDMD
for i=1:N
    X{i}=X_f(2*i-1:2*i,:);
    Y{i}=Y_f(2*i-1:2*i,:);

    X_g{i}= X_f_g(2*i-1:2*i,:);
    Y_g{i}=Y_f_g(2*i-1:2*i,:);

    % Topology is unkonwn
    tt=1:n;
    tt2=tt(tt~=2*i-1);
    tt3=tt2(tt2~=2*i);
    ind_topo{i}=tt3;
    U{i}=X_f(ind_topo{i},:);
end


U_t{1}=U{1};  % Not Used
U_t{2}=X_f(input_ind,:);
U_t{3}=X_f(input_ind,:);



U_t_g{1}=X_f_g(ind_topo{1},:);
U_t_g{2}=X_f_g(input_ind,:);
U_t_g{3}=X_f_g(input_ind,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% sEDMD (Toplogy is known)  (1, 1--2, 1--3)

X_c{1}=X_f(1:2,:);
Y_c{1}=Y_f(1:2,:);

X_c{2}=X_f(1:4,:);
Y_c{2}=Y_f(1:4,:);

X_c{3}=X_f([1 2 5 6],:);
Y_c{3}=Y_f([1 2 5 6],:);

fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ******************************* Lift ***********************************

disp('Starting LIFTING')
tic
Xlift_f = liftFun_f(X_f);
Ylift_f = liftFun_f(Y_f);

for i=1:N
    Xlift{i} = liftFun{i}(X{i});
    Ylift{i} = liftFun{i}(Y{i});


    Xlift_g{i} = liftFun{i}(X_g{i});
    Ylift_g{i} = dliftFun{i}(X_g{i},Y_g{i});

end

%%% sEDMD
Xlift_c{1}=liftFun_c{1}(X_c{1});
Xlift_c{2}=liftFun_c{2}(X_c{2});
Xlift_c{3}=liftFun_c{3}(X_c{3});
Ylift_c{1}=liftFun_c{1}(Y_c{1});
Ylift_c{2}=liftFun_c{2}(Y_c{2});
Ylift_c{3}=liftFun_c{3}(Y_c{3});


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

%%% topology is known
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

%% mgEDMD
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
%%% topology is known
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


%% sEDMD
for i=1:N
    Wn_c = Ylift_c{i};
    Vn_c = Xlift_c{i};
    VVtn_c = Vn_c*Vn_c';
    WVtn_c = Wn_c*Vn_c';
    Mn_c{i} = WVtn_c * pinv(VVtn_c);  % [A,B1,...Bn]
    Alift_c{i}=Mn_c{i};

end


fprintf('Regression done, time = %1.2f s \n', toc);




%% ********************** error comparision
% EDMD vs sEDMD vs lEDMD vs mEDMD vs mgEDMD
if 1
    close all
    clf
    Tmax = 0.5;
    Nsim = Tmax/deltaT;
    Nsim_g = Tmax/deltaTg;

    if exist("X0")
        Nsample=size(X0,2);
    else
        Nsample=500;
        X0=rand(n, Nsample)*0.6-0.6/2;
    end

    EMAX_f=[];
    EMAX_t=[];
    EMAX_t_g=[];
    EMAX_c=[];
    EMAX_lin_t=[];
    for Ni=1:Nsample
        x0 = X0(:,Ni);
        x_true = x0;
        x_true_g = x0;
        % Lifted initial condition

        for i=1:N
            xlift{i}= liftFun{i}(x0(2*i-1:2*i));  
            xlift_t{i}=xlift{i};  % mEDMD
            xlift_t_g{i}=xlift{i};  % mgEDMD
            xlift_lin_t{i}=xlift{i};  % lEDMD
        end
        xlift_f =  liftFun_f(x0);  % EDMD
        % sEDMD
        xlift_c{1}= liftFun_c{1}(x0(1:2));
        xlift_c{2}= liftFun_c{2}(x0(1:4));
        xlift_c{3}= liftFun_c{3}(x0([1,2,5,6]));
   

        % Simulate
        x_n_old=[];
        x_n_old_t=[];
        x_n_old_t_g = [];
        x_lift_old = [];
        x_lift_old_t=[];
        for i = 1:Nsim
            % Koopman predictor

            %% EDMD
            xlift_f=[xlift_f, Alift_f*xlift_f(:,end)];
            %% mEDMD and lEDMD


            x_n_old_t=xlift_t{1}(input_ind,end);  % mEDMD
            x_lift_old_t=xlift_lin_t{1}(:,end);  % lEDMD


            for ii=1:N

 
                if ii==1
                    xlift_t{ii} = [xlift_t{ii}, Alift_t{ii}*xlift_t{ii}(:,end)]; %mEDMD
                    xlift_lin_t{ii}=[xlift_lin_t{ii}, Alift_lin_t{ii}*xlift_lin_t{ii}(:,end)];  % lEDMD
                else
                    xlift_t{ii} = [xlift_t{ii}, Alift_t{ii}*xlift_t{ii}(:,end) + Blift_t{ii}*kron(eye(length(input_ind)), xlift_t{ii}(:,end))*x_n_old_t];  % mEDMD
                    xlift_lin_t{ii}=[xlift_lin_t{ii}, Alift_lin_t{ii}*xlift_lin_t{ii}(:,end)+Blift_lin_t{ii}*x_lift_old_t];  % lEDMD
                end
           
                xlift_c{ii} =[xlift_c{ii}, Alift_c{ii}*xlift_c{ii}(:,end)]; % sEDMD
            end
            % True dynamics
            x_true = [x_true, f_0d(0,x_true(:,end)) ];


        end
        %% sim for mgEDMD
        for i=1:Nsim_g
            x_n_old_t_g = xlift_t_g{1}(input_ind,end);
            for ii=1:N
                if ii==1
                    xlift_t_g{ii} = [xlift_t_g{ii}, f_g_md{ii}(0,xlift_t_g{ii}(:,end))];
                else
                    xlift_t_g{ii} = [xlift_t_g{ii}, f_g_md{ii}(0,xlift_t_g{ii}(:,end), x_n_old_t_g)];
                end
            end
            % True dynamics
            x_true_g = [x_true_g, f_0dg(0,x_true_g(:,end)) ];

        end
        %% use x_true_g as true traj.

        x_true = x_true_g (:, 1:deltaT/deltaTg:Nsim_g+1);

        vec=1:n;
        x_koop_f = xlift_f(vec,:); % Koopman predictions  


        for ii=1:N
    
            x_koopn_t{ii}=xlift_t{ii}(1:n/N,:);   % mEDMD
            x_koopn_lin_t{ii}=xlift_lin_t{ii}(1:n/N,:);  % lEDMD
            x_koopn_t_g{ii}=xlift_t_g{ii}(1:n/N,:);  %mgEDMD
            x_koopn_c{ii}=xlift_c{ii}(1:n/N,:);  % sEDMD
        end



        emax_f = max(abs(x_koop_f-x_true),[],2);
        emax_t = max(abs([x_koopn_t{1};x_koopn_t{2};x_koopn_t{3}]-x_true),[],2);
        emax_t_g = max(abs([x_koopn_t_g{1};x_koopn_t_g{2};x_koopn_t_g{3}]-x_true_g),[],2);
        emax_c = max(abs([x_koopn_c{1};x_koopn_c{2};x_koopn_c{3}]-x_true),[],2);
        emax_lin_t = max(abs([x_koopn_lin_t{1};x_koopn_lin_t{2};x_koopn_lin_t{3}]-x_true),[],2);


        EMAX_f = [EMAX_f  emax_f];
        EMAX_t = [EMAX_t emax_t];
        EMAX_t_g = [EMAX_t_g emax_t_g];
        EMAX_c = [EMAX_c emax_c];
        EMAX_lin_t = [EMAX_lin_t emax_lin_t];

    end


    switch basisFunction
        case 'rbf'
            file_name= ['Datensatz\EMAX_all_', num2str(Ntraj),'_rbf','.mat'];
        case 'mono'
            file_name= ['Datensatz\EMAX_all_', num2str(Ntraj),'_mono','.mat'];
    end
    save(file_name,'EMAX_c','EMAX_t', 'EMAX_f', 'EMAX_lin_t', 'EMAX_t_g');


    %% % overall plot

    psize = [14.65, 9.6];
    ppos  = [0.3, 0.2, psize(1)-0.3 psize(2)-0.2];
    gcaPos = [0.1 0.1 0.85 0.75];

    FTsize_axisnumbers = 14;
    Legend_size = 18;
    Title_size = 18;

    load(file_name);
    Temp_comb=[EMAX_f; EMAX_c;  EMAX_lin_t; EMAX_t; EMAX_t_g]; % without EMAX(for unknown topology)
    Nsample=size(EMAX_t,2);
    ni=2;
    N=3;
    EMAX_all=zeros(N*5, size(Temp_comb,2));
    Tmax=0.5;
    for i=1:N*5
        %     EMAX_all(i,:)=log(sqrt(Temp_comb(2*i-1,:).^2+Temp_comb(2*i,:).^2));
        EMAX_all(i,:)=log ( abs(Temp_comb(2*i-1,:))+abs(Temp_comb(2*i,:)) );
    end
    g = kron([1:N*5]', ones(Nsample,1));
    pos=[1 1.4 1.8  2.4 2.8 3.2  3.8 4.2 4.6  5.2 5.6 6  6.6 7 7.4];
    pos2 = mean([pos(3:4);pos(6:7); pos(9:10); pos(12:13)],2);
    col=repmat({'g','c','w'}, 1,5);
    posmean=pos([2:3:3*5-1 ]);
    row1 = {'EDMD', 'sEDMD',  'lEDMD','mEDMD', 'mgEDMD'};
    tickLabels = row1;
    figure(22)
    boxplot([reshape(EMAX_all', [numel(EMAX_all),1])],g, 'positions',pos );
    set(gca,'xtick',posmean,'Ylim', [-11, 0]);  %, , 'Ylim', [-11, 1]
    set(gca,'xticklabel',tickLabels, 'TickLabelInterpreter', 'latex');
    set(gca,'YScale','linear');
    grid on
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),col{j},'FaceAlpha',.6);
    end
    c = get(gca, 'Children');
    hold on
    xline(pos2,':');
    set(gcf,'PaperUnits','centimeters','PaperSize',psize,'PaperPosition',ppos)
    set(gca,'Position',gcaPos,'fontsize',FTsize_axisnumbers);

    title(['$\ln(\|\Delta x_i\|_{\infty})$, ', num2str(Ntraj), ' snapshots'],'interpreter','latex', 'FontSize', Title_size);
    hleg1 = legend(c(1:3), ' $ x_1 $', '$x_2$', '$x_3$', 'Interpreter', 'latex','FontSize', Legend_size, 'Orientation','horizontal');

    set(gca, 'XGrid','off', 'YGrid','on');


    if printing
        name = sprintf(['Overall_Duffing_', num2str(Ntraj), '_new_wg_rbf']);
        print([folder name],'-depsc','-painters')
        print([folder name],'-dpdf','-painters')
        print([folder name],'-dbmp')
        savefig([folder name])

    end




end

