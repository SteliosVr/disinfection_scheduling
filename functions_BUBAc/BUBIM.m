function [Hb] = BUBIM(d,inputNodesCLid,sensorNodesCLid,stateEstim,Kdunc)

%% Network quality setting check:
% NetworkQualitySettingsCheck(inpname);

%%
bita=Kdunc; 

Q = stateEstim.Flow;
Ql = stateEstim.FlowLB;
Qu = stateEstim.FlowUB;
H = stateEstim.Head;
CL = stateEstim.NodeQuality;

% d=epanet(inpname);
sim_time=double(d.getTimeSimulationDuration);
sensorNodesCL = double(d.getNodeIndex(sensorNodesCLid)); %sensor nodes

%% Pipe Area Ar, Pipe length L
diam=d.getLinkDiameter/1000; %in meters
radius = diam/2;
Ar = pi*(diam./2).^2; %pipe area= pi*r^2 in (meters^2)
L = d.getLinkLength; %pipe lenght in meters

%% Tank and reservoir parameters
TankID=d.getNodeTankIndex;
ReservoirID=d.getNodeReservoirIndex;
TankV=(H(:,TankID)- ones(size(H(:,1)))*d.NodeElevations(TankID)).*(pi*ones(size(H(:,1)))*(d.NodeTankDiameter(TankID)/2).^2); %

%% Create incidence matrix (topology graph)
node_node_link_index=[d.NodesConnectingLinksIndex d.LinkIndex'];
node_node_link_name=[d.NodesConnectingLinksID d.LinkNameID'];
Ain=zeros(d.NodeCount,d.LinkCount); %size nodes x pipes
for i=1:length(d.NodesConnectingLinksIndex)
    Ain(node_node_link_index(i,1),node_node_link_index(i,3))=1;
    Ain(node_node_link_index(i,2),node_node_link_index(i,3))=-1;
end

%% Kb, Kw, K, Ke, Ku, Kl - Reaction rates
%%%% Kb - Bulk reaction coefficient in hours
Kb = d.getLinkBulkReactionCoeff./24; % Bulk reaction coefficient (1/hour) 
d.QualityReactionCoeffBulkUnits; % Bulk reaction coefficient units

%%%% Kw - Wall reaction coefficient in hours
Kw = d.getLinkWallReactionCoeff./24; %Wall reaction coefficient (m/hour)
d.QualityReactionCoeffWallUnits; % Wall reaction coefficient units

%%%% Ktank - Tank reaction coefficient in hours
KTank=d.getNodeTankBulkReactionCoeff;
KTank = KTank./24;
KlTank=KTank - bita*KTank;
KuTank=KTank + bita*KTank;

%%%% Kf - Calculation of mass tranfer coefficient
Kfw = (4./diam).*Kw; % page 44 epanet manual

%%% K - Total chlorine decay rate:
K = (Kb + Kfw);

%%%% Kl, Ku - Bounds on total chlorine decay rate 
Kl = K + bita.*K; %negative plus negative = most negative
Ku = K - bita.*K; %negative minus negative = less negative

%% Set simulation steps and Convert simulation time steps from seconds to hours 
th = double(d.getTimeHydraulicStep)/3600; %Hydraulic time step in hours
tq = double(d.getTimeQualityStep) /3600; %Quality Step in hours
start_step=1;
end_step=(sim_time/3600/tq)+1; %Calculate last time step

%% Replicate missing hydraulic steps
% if th~=tq
if size(CL,1)~=size(Ql,1)
n=th/tq;
TankV=[repelem(TankV(1:end-1,:)',1,n) TankV(end,:)']';
n=(end_step-1)/(size(Ql,1)-1);
Ql=[repelem(Ql(1:end-1,:)',1,n) Ql(end,:)']';
Qu=[repelem(Qu(1:end-1,:)',1,n) Qu(end,:)']';
% Ql=[Ql(1,:)' repelem(Ql(2:end,:)',1,n) ]';
% Qu=[Qu(1,:)' repelem(Qu(2:end,:)',1,n) ]';
end

%% Find input nodes:
inputN=d.getNodeIndex(inputNodesCLid);
% type=d.getNodeSourceType;
% for i=1:length(type)
%     if strcmp(type{i},'SETPOINT')
%         inputN = [inputN i];
%     end
% end

%% Constant Parameter and Calculated hydraulic states struct:
times = struct('QualityStep',tq, 'HydraulicStep',th, 'SimulationTime',sim_time);
links = struct('FlowUpper',Qu, 'FlowLower',Ql, 'Length',L, 'Area',Ar, 'DecayRateUpper',Ku, 'DecayRateLower',Kl);
tanks = struct('ID',TankID, 'Volume',TankV, 'DecayRateUpper',KuTank, 'DecayRateLower',KlTank);
nodes = struct('TankInfo',tanks, 'ReservoirID',ReservoirID, 'ActualClConcen',CL,'inputN',inputN);
sim_param = struct('time',times, 'links',links, 'tanks',tanks, 'nodes',nodes, 'IncidenceMat',Ain);

%% Find output node min-max concentration
i=1;
starttime=tic;
Max_Loop_Time=0;
nn=d.getNodeCount;
m=length(Q);
Hb = zeros(m,m);
for Na=sensorNodesCL
for kc=start_step:end_step
%     kc
%     loopstarttime=tic;
    [Nin kin impin ratioin pipein] = Algorithm_5c(Na,kc,sim_param); 
    Hb(kc,kin)=ratioin.*exp(impin);% impin contains: K*delay*tq
%     Hb(unique([kc kin'], 'rows'))=1;
%     clc
%     %%%%%%%%%%%%%%%%%% Print times  %%%%%%%%%%%%%%
%     Elapsed_Time=toc(starttime);
%     loop_Elapsed_Time(i)=toc(loopstarttime);
%     if i>10
%     Max_Loop_Time=max(Max_Loop_Time,loop_Elapsed_Time(i));
%     end
%     fprintf('The program elapsed time is:  %i:%.2i  minutes:seconds \n\n', (floor(Elapsed_Time/60)),floor(mod(Elapsed_Time,60)));
%     fprintf('The current  loop   time is:  %.2f  seconds\n\n', loop_Elapsed_Time(i));
%     Average_Loop_Time=mean(loop_Elapsed_Time);
%     fprintf('The average  loop   time is:  %.2f  seconds\n\n', Average_Loop_Time);
%     fprintf('The maximum  loop   time is:  %.2f  seconds\n\n', Max_Loop_Time);
%     fprintf('Node number : --> %i <-- of   %i   nodes\n\n', Na,nn);
%     fprintf('Time step is: --> %i <-- of   %i   steps\n\n', kc,end_step);
%     est_rem_time=Average_Loop_Time*(end_step-kc)*(nn-Na);
%     est_rem_time_min=floor(est_rem_time/60);
%     est_rem_time_sec=floor(mod(est_rem_time,60));
%     fprintf('The estimated remaining time is:  %i:%.2i  minutes:seconds \n\n', est_rem_time_min,est_rem_time_sec);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
i=i+1;
end
end
% fprintf('\nBUBA Simulation Finished\n\n');

%%


%% reshape results:
% cnmax=reshape(c_calc(:,2,:),size(c_calc,1),size(c_calc,3));
% cnmin=reshape(c_calc(:,1,:),size(c_calc,1),size(c_calc,3));

%% get only sensor nodes:
% cnmax = cnmax(sensorNodesCL,:);
% cnmin = cnmin(sensorNodesCL,:);

%%
% d.unload
end

