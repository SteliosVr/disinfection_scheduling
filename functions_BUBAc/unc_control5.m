function [Ugurobi] = unc_control5(Hn,H1,H2,m,r)
H1diag=zeros(m*m,m);
H2diag=zeros(m*m,m);
T1=zeros(m*m,m);
T2=zeros(m*m,m);
row=1;
for i=1:m
    H1diag(row:i*m,1:m)=diag(H1(i,:));
    T1(row:i*m,i)=1;
    H2diag(row:i*m,1:m)=diag(H2(i,:));
    T2(row:i*m,i)=1;
    row=row+m;
end

H1zind=find(all(H1diag'==0));
H1diag(H1zind,:)=[];
T1(H1zind,:)=[];
H2zind=find(all(H2diag'==0));
H2diag(H2zind,:)=[];
T2(H2zind,:)=[];

%% Formulate Gurobi
% reference and output constraints
R=r;
ymin = r-0.0001;

% smooth input
Itemp = [zeros(m,1) eye(m)]; 
Itemp(:,end)=[];
Itilda=eye(m)-Itemp; % with \Delta u penalization

% objective
QHz  = [Hn zeros(m,2*m)];
DelU = 0.5*[Itilda zeros(m,2*m)];
Qobj = (QHz'*QHz) + (DelU'*DelU);
bobj = 2*QHz'*(-R);
% bobj = 2*(-R)'*QHz;

% inequality constraints
m1=size(H1diag,1);
m2=size(H2diag,1);

Aineq = [-H1diag    T1            zeros(m1,m);
         -H2diag    zeros(m2,m)   T2        ;   
         zeros(m)   -eye(m)       -eye(m)];
    
bineq =[zeros(m1,1);
        zeros(m2,1);
        -ymin];

    
% lower and upper bounds
xl = [zeros(m,1); zeros(size(Aineq,2)-m,1)];
xu = [ones(m,1);  max(H1')'; max(H2')'];

% gurobi model
% names = {'u', 't11', 't12', 't21', 't22'};
% model.varnames = names;
model.Q = sparse(Qobj);
model.obj = bobj;
model.A = sparse([Aineq]);
model.rhs = [bineq];
model.sense = '<';
model.vtype = 'C';
model.lb = xl;
model.ub = xu;
% params.LogFile = 'gurobiLogFile';
% params.method = 1;
% params.OutputFlag = 1;
% gurobi_write(model, 'qp.lp'); % Writes the model to a file.
result = gurobi(model);

Ugurobi=result.x(1:m);

t1=(result.x(m+1:2*m));
t2=(result.x(2*m+1:3*m));
figure
plot(Hn*Ugurobi,'Linewidth',1.2)
hold all
plot(r)
plot(Ugurobi)
plot(t1);plot(t2)
plot(t1+t2)
legend('y','r','u','t1','t2','t1+t2')


end