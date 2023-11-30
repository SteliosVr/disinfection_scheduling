function [ Nin, kin, impin, ratioin, pipein] = Algorithm_5c(Na, kc, sim_param)
%ALGORITHM_5 

%% if time step exceeds initial time step:
if kc<=1
    Nin=[];
    kin=[];
    impin=[];
    ratioin=[];
    pipein=[];
    return
end

%% terminal condition reached - input node
if ismember(Na,sim_param.nodes.inputN)
    Nin=Na;
    kin=kc+1;
    impin=0;
    ratioin=1;
    pipein=[];
    return
end

%% if node is a reservoir:
% ResID = sim_param.nodes.ReservoirID;
% %concentration is the same as the initial concentration
% if ismember(Na, ResID(:))
%     cnmin=c_calc(Na,1,1);
%     cnmax=c_calc(Na,2,1);
%     return
% end

%% if node is a tank:
% TankID = sim_param.tanks.ID;
% if ismember(Na, TankID(:))
%     [cnmin, cnmax]= Algorithm_4( Na, kc, sim_param );
%     return
% end

%% Find all pipes that bring water into node Na:
Ain = sim_param.IncidenceMat;
Qu = sim_param.links.FlowUpper;
Ql = sim_param.links.FlowLower;
pin=[];
pin1=find(Ain(Na,:)==-1); %according to convention
for l=pin1
    if (Qu(kc,l)>0 || Ql(kc,l)>0)
        pin=[pin; l];
    end
end
pin2=find(Ain(Na,:)==1); %opposite to convention
for l=pin2
    if (Qu(kc,l)<0 || Ql(kc,l)<0)
        pin=[pin; -l];
    end
end


%% Find possible upstream node Naus or Na and water detention times
Nin=[]; kin=[]; impin=[]; ratioin=[]; pipein=[];
tq = sim_param.time.QualityStep;
Kl = sim_param.links.DecayRateLower;
Ku = sim_param.links.DecayRateUpper;
Q  = (Ql+Qu)/2;
for l = pin'
    if l<0;keyboard; end
    [ Naus, delay] = Algorithm_1c( Na, abs(l), kc, sim_param ); 
    temp = unique([Naus delay], 'rows');
    Naus=temp(:,1); delay=temp(:,2); 
    impactL = Kl(l)*tq*(delay);
    impactU = Ku(l)*tq*(delay);
    impact=(impactL+impactU)/2;
    ratio=(Q(kc,l)/sum(Q(kc,abs(pin))));
    pipe = l;
    
    for i= 1:length(Naus)
        Ntemp=Naus(i);d=delay(i);
        [ Nintemp, kintemp, impactemp, ratiotemp, pipetemp] = Algorithm_5c(Ntemp, kc-d, sim_param);
        Nin = [Nin Nintemp];
        kin = [kin kintemp];
        impin = [impin impactemp+impact(i)];
        ratioin = [ratioin ratiotemp*ratio];
        pipein = [pipein pipetemp repmat(pipe,size(Nintemp))];
    end
    
end





end

