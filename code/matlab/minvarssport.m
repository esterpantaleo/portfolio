%minimum variance short selling portfolio
function [PortRisk,PortWts,mineigenvalue]=minvarssport(ExpCovariance,optimization_algorithm)
azioni=length(ExpCovariance);

switch optimization_algorithm
    case 'quadprog'
        w0=ones(azioni, 1)/azioni;
        localmaxiter = size([],1) + size([],1);      % <- change
        localmaxiter = max(localmaxiter,azioni);
        localmaxiter = 10 * max(localmaxiter,20);
        Aineq=[];
        Bineq=[];
        LB=[];
        UB=[];
        options=optimset('largescale', 'off', 'MaxIter', localmaxiter);

        v0=zeros(1,azioni)';
        v1=ones(1,azioni);
        PortWts=quadprog(ExpCovariance,v0,Aineq,Bineq,v1,1,LB,UB,w0,options);
        PortWts=PortWts';
        clear w0 localmaxiter Aineq Bineq LB UB options v0 v1
        
    case 'tola'
        covarianzaorlata=[(2*ExpCovariance) ones(azioni,1);ones(1,azioni) 0];
        b=zeros(azioni+1,1);
        CC=pinv(covarianzaorlata);
        b(end,1)=1;
        pes=CC*b;
        PortWts(1,:)=pes(1:azioni,1);
%         [V,D]=eig(ExpCovariance);
%         d=diag(D);
%         d(find(d>-0.0000001 & d<0))=0;
%         D=diag(d);
%         ExpCovariance=V*D*V';
        clear covarianzaorlata D V b pes CC
    case 'markowitz'
        A=pinv(ExpCovariance)*ones(azioni,1);%pinv
        L=pinv(ExpCovariance);
        save PInvExpCovariance2007novemberS L -ascii
        PortWts=A/(ones(1,azioni)*A);
        PortWts=PortWts';
        clear A
end

mineigenvalue=min(min(eig(ExpCovariance)));
%save ExpCovariance2007novemberS ExpCovariance -ascii
    
% matrix=V'*PortWts';
% PortRisk=sqrt(matrix'*D*matrix);
PortRisk=sqrt(PortWts(1,:)*ExpCovariance*PortWts(1,:)');