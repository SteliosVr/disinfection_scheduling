fclose all;clc;clear;clear class;
try
d.unload
catch exception
end
clear exception
addpath(genpath(pwd)); 

%% EPANET Input File
% inpname='networks/Net1_Rossman2000_S_CMH.inp';
% inpname='networks/Net1_CMH.inp';
inpname = enterNetwork([]);
d=epanet(inpname);

%% Change simulation time and simulation time steps
% sim_time= 5*24*60*60; % in seconds
% th = 3*60; % hydraulic step in seconds
% tq = 3*60; % quality step in seconds
% d.setTimeSimulationDuration(sim_time)
% d.setTimeHydraulicStep(th); %Hydraulic time step
% d.setTimeQualityStep(tq); %Quality Step

sim_time=double(d.getTimeSimulationDuration);
th=double(d.getTimeHydraulicStep); %Hydraulic time step
tq=double(d.getTimeQualityStep); %Quality Step

%% Set uncertainty and input/output nodes
alfa = 0.02;  % water flow uncertainty
bita = 0.02;  % decay rate uncertainty
inputN1 = 1;  % define input node (only one)
inputN2 = 1;  % define input node (only one)
Na = 9;      % define output node (only one)

% d.plot('nodes','yes','fontsize',12)

%% Change quantities related to Hydraulic simulation
% d.setTimePatternStep(1*60*60); %Demand pattern change time step.
% d.setOptionsPatternDemandMultiplier(2*0.2*rand(1)+0.8)       

%% Solve Hydraulics
d.solveCompleteHydraulics
%d.saveHydraulicFile([pwd,'\hydraulics.hyd'])
%d.useHydraulicFile([pwd,'\hydraulics.hyd'])

%% Initialize Quality simulation - Define chlorine input node
d.setQualityType('chem','mg/L')
d.openQualityAnalysis
d.initializeQualityAnalysis
tleft=1; T=[]; Q=[]; yq=[]; H=[];
for i=1:d.getNodeCount %all nodes
d.setNodeSourceType(i,'CONCEN')
%if you don't set all nodes' type then it gives NAN
end
d.setNodeSourceType(inputN1,'SETPOINT') %input node type
d.setNodeSourceType(inputN2,'SETPOINT') %input node type
% value=d.getNodeInitialQuality;
% value=value+0.5;
% d.setNodeInitialQuality(value);

%% Set decay rate coefficients
%Bulk reaction coefficients
% bita=0.5; % decay coefficient uncertainty
value=d.getLinkBulkReactionCoeff;
value= value + (2*bita*rand(size(value))-bita).*value;
d.setLinkBulkReactionCoeff(value);
%Wall reaction coefficients
value=d.getLinkWallReactionCoeff;
value=zeros(size(value));
d.setLinkWallReactionCoeff(value);
%Tank reaction coefficient
KTank=d.getNodeTankBulkReactionCoeff;
% d.setNodeTankBulkReactionCoeff(KTank);

%% Solve Quality Step-by-step
i=1;
while (tleft>0)
    %Add code which changes something related to quality
    t=d.runQualityAnalysis;
    Q=[Q; d.getLinkFlows]; %water flows in pipes
    yq=[yq; d.getNodeActualQuality]; %water quality at nodes
    H=[H; d.getNodeHydaulicHead];
    values = d.getNodeSourceQuality;
% change input value here
%     values(inputN1)=values(inputN1)+0.1*sin(0.1*i) +0.1*cos(0.2*i); 
    values(inputN1)=1;
    values(inputN2)=1;
    d.setNodeSourceQuality(values);
    T=[T; t]; 
    tleft = d.stepQualityAnalysisTimeLeft;
    i=i+1;
end
d.closeQualityAnalysis;

%% Convert simulation time steps from seconds to hours 
th = double(d.getTimeHydraulicStep)/3600; %Hydraulic time step in hours
tq = double(d.getTimeQualityStep) /3600; %Quality Step in hours

%% Pipe Area Ar, Pipe length L
diam=d.getLinkDiameter/1000; %in meters
Ar = pi*(diam./2).^2; %pipe area= pi*r^2 in (meters^2)
L = d.getLinkLength; %pipe lenght in meters

%% Tank and reservoir parameters
TankID=d.getNodeTankIndex;
ReservoirID=d.getNodeReservoirIndex;
TankV=(H(:,TankID)- ones(size(H(:,1)))*d.NodeElevations(TankID)).*(pi*ones(size(H(:,1)))*(d.NodeTankDiameter(TankID)/2).^2); %

%% Create incidence matrix (topology graph)
%row index are nodes, column index are pipes
%a value of -1 means that the pipe of that column brings water into the
%node of that row
%Note that flows with negative signs means that the flow is in the opposite
%direction to the direction in which the pipe was drawn initially.
node_node_link_index=[d.NodesConnectingLinksIndex d.LinkIndex'];
node_node_link_name=[d.NodesConnectingLinksID d.LinkNameID'];
%initial flow direction convention: 
%first column = 1 (water leaves the node),
%second column = -1 (water enters the node).
Ain=zeros(d.NodeCount,d.LinkCount); %size nodes x pipes
for i=1:length(d.NodesConnectingLinksIndex)
    Ain(node_node_link_index(i,1),node_node_link_index(i,3))=1;
    Ain(node_node_link_index(i,2),node_node_link_index(i,3))=-1;
end

%% Water Flow Bounds Qu, Ql
% Generate upper Qu and lower Ql bound of flows.
% d.LinkFlowUnits
% d.LinkVelocityUnits

% d.getLinkVelocity
% alfa = 0.1; %flow uncertainty
Ql = Q - alfa.*abs(Q); 
Qu = Q + alfa.*abs(Q);
% randomness=rand(size(Q));
% Ql = Q - alfa*randomness.*abs(Q);
% Qu = Q + alfa*(1-randomness).*abs(Q);
% Ql = Q - alfa.*abs(Q);
% Qu = Q + alfa.*abs(Q);
% Ql = Q - alfa* rand(size(Q,1),1)*(mean(Q)/2);
% Qu = Q + alfa* rand(size(Q,1),1)*(mean(Q)/2) ;
% load Qbounds; Ql=Qmin; Qu=Qmax;
% Generate flow estimate by adding random error on the real flows keeping 
% it between the bounds.
% Qe = Q + (2*alfa*variability-alfa).*Q;


%% Kb, Kw, K, Ke, Ku, Kl - Reaction rates
%%%% Kb - Bulk reaction coefficient in hours
Kb = d.getLinkBulkReactionCoeff./24; % Bulk reaction coefficient (1/hour) 
d.QualityReactionCoeffBulkUnits; % Bulk reaction coefficient units

%%%% Kw - Wall reaction coefficient in hours
Kw = d.getLinkWallReactionCoeff./24; %Wall reaction coefficient (m/hour)
d.QualityReactionCoeffWallUnits; % Wall reaction coefficient units

%%%% Ktank - Tank reaction coefficient in hours
KTank = KTank./24;
KlTank=KTank - bita*KTank;
KuTank=KTank + bita*KTank;

%%%% Kf - Calculation of mass tranfer coefficient
% NOTE: mass transfer coefficient depends on flow!
visc=8.9*10^(-3)*3600;%kinematic viscosity of water (Pascal/hour at 25degress Celsious)
D=1.21*10^(-9)*3600;%molecular diffusivity of chlorine in water (m^2/hour)
Sc=visc/D;
Kf=zeros((length(L)-1),1);
for i=1:(length(L)-1)
R=sqrt(norm(Q(:,i)))/Ar(i)*diam(i)/visc;
if R>=2300
    Sh=0.0149*R^0.88*Sc^0.333;
else
    Sh=3.65+(0.0668*(diam(i)./L(i)*(R*Sc))/(1+0.04*(diam(i)/L(i)*R*Sc)^(0.667)));
end
Kf(i)=Sh*D/diam(i); %Mass transfer coefficient [Rossman 1994]
end
Kf=[Kf' 0];

%%%% K - Total chlorine decay rate
K = Kb ;%- (Kf.*Kw)./(Ar.*(Kf+Kw)); %All pipe reaction rate

%%%% Kl, Ku, Ke - Bounds on total chlorine decay rate 
Kl = K + bita.*K; %negative plus negative = most negative
Ku = K - bita.*K; %negative minus negative = less negative
Ke = K + (2*bita*rand(size(K))-bita).*K;

%% Set active node Na, simulation time
start_step=10;
end_step=(sim_time/3600/tq); %Calculate last time step

% Na=9; %Set output node

% fault=real_conc.*[zeros(length(real_conc)-600,1); -0.2*ones(600,1)];
% noise=-0.05+0.1*rand(length(real_conc),1);
% sens_conc=real_conc +noise +fault;

%% Initialize calculated bounds and set as global variable
global c_calc
c_calc=NaN(d.NodeCount,2,end_step);
c_calc(:,:,1)=[yq(1,:)' yq(1,:)']; %initial conditions
% c_calc(:,:,1)=zeros; %initial conditions
c_calc(inputN1,1,1:end_step)=yq(1:end_step,inputN1); %input node minimum
c_calc(inputN1,2,1:end_step)=yq(1:end_step,inputN1); %input node maximum
c_calc(inputN2,1,1:end_step)=yq(1:end_step,inputN2); %input node minimum
c_calc(inputN2,2,1:end_step)=yq(1:end_step,inputN2); %input node maximum
% c_calc(11,1,1:end_step)=yq(1:end_step,11); %tank
% c_calc(11,2,1:end_step)=yq(1:end_step,11); %tank

%% Find output node min-max concentration
% for Na=1:9
real_conc=yq(1:end_step,Na);
i=1;
starttime=tic;
Max_Loop_Time=0;
cnmin=zeros(1,end_step);
cnmax=zeros(1,end_step);

%%%%%%% Constant Parameter and Calculated hydraulic states struct
times = struct('QualityStep',tq, 'HydraulicStep',th, 'SimulationTime',sim_time);
links = struct('FlowUpper',Qu, 'FlowLower',Ql, 'Length',L, 'Area',Ar, 'DecayRateUpper',Ku, 'DecayRateLower',Kl);
tanks = struct('ID',TankID, 'Volume',TankV, 'DecayRateUpper',KuTank, 'DecayRateLower',KlTank);
nodes = struct('TankInfo',tanks, 'ReservoirID',ReservoirID, 'ActualClConcen',yq);
sim_param = struct('time',times, 'links',links, 'tanks',tanks, 'nodes',nodes, 'IncidenceMat',Ain);

for kc=start_step:end_step
    loopstarttime=tic;
    [cnmin(kc), cnmax(kc)] = Algorithm_5(Na,kc,sim_param);
    clc
    %%%%%%%%%%%%%%%%%% Print times  %%%%%%%%%%%%%%
    Elapsed_Time=toc(starttime);
    loop_Elapsed_Time(i)=toc(loopstarttime);
    if i>10
    Max_Loop_Time=max(Max_Loop_Time,loop_Elapsed_Time(i));
    end
    fprintf('The program elapsed time is:  %i:%.2i  minutes:seconds \n\n', (floor(Elapsed_Time/60)),floor(mod(Elapsed_Time,60)));
    fprintf('The current  loop   time is:  %.2f  seconds\n\n', loop_Elapsed_Time(i));
    Average_Loop_Time=mean(loop_Elapsed_Time);
    fprintf('The average  loop   time is:  %.2f  seconds\n\n', Average_Loop_Time);
    fprintf('The maximum  loop   time is:  %.2f  seconds\n\n', Max_Loop_Time);
    fprintf('Time step is: --> %i <-- of   %i   steps\n\n', kc,end_step);
    est_rem_time=Average_Loop_Time*(end_step-kc);
    est_rem_time_min=floor(est_rem_time/60);
    est_rem_time_sec=floor(mod(est_rem_time,60));
    fprintf('The estimated remaining time is:  %i:%.2i  minutes:seconds \n\n', est_rem_time_min,est_rem_time_sec);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
i=i+1;
end
fprintf('\nSimulation Finished\n\n');

%% plot results
minutes=1:(tq*60):(sim_time/60);
hours=1-59/60:(tq):(sim_time/3600);
figure
plot(hours,cnmax,'r','LineWidth',1.5)
hold all
plot(hours,real_conc,'Color',[0.9 0.9 0],'LineWidth',1.5)
plot(hours,cnmin,'b','LineWidth',1.5)
legend('max Cl concentration','EPANET Cl concentration','min  Cl concentration','Location','northeast')
xlabel('Time (hours)')
ylabel('Cl (mg/L)')
axis([-inf inf -inf 1])
grid on
% end

%% save data
% savefile = sprintf('Net_%.0f-Na_%d-alpha_%.2f-bita_%.2f-step_%.0f-SimTime_%.1f.mat',net,Na,alfa,bita,tq*60,sim_time/24/3600 )
% save(savefile);
% save(['saved_simulations/' savefile '.mat']);
d.unload