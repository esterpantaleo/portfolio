function [rischiopredetto,rischiorealizzato,alpha,rischiopercentuale,pesi,minimoautovalore,rendimentorealizzato] = efficient_frontier(algoritmo,rend,shortselling,giorno_iniziale,giorniA,giorniB)

%determina il portafoglio di minima varianza e i rischi con un portafoglio selezionato secondo
%l'algoritmo:
%'Markowitz','PseudoinverseTola','PseudoinverseMarkowitz','RMTStanley','RMTBouchaud','RMTCovarianceStanley','RMTCovarianceBouchaud','SI',
%'SIStanleyFirstEigenvalue','MSTSingle','CLCA_c','CLCA_s','MSTAverage','MSTWeighted','MSTHausdorff_s','MSTHausdorff_h','MSTAverage_correlation',
%'MSTWeighted_correlation','ShrinkageSI','ShrinkagePerfectPositiveCorrelation','ConstantCorrelation','ShrinkageConstantCorrelation','DiagonalCommonVariance',
%'ShrinkageDiagonalCommonVariance','ShrinkageCommonCoVariance','ShrinkageDiagonalUnequalVariance','ShrinkageSIS_Pw','ShrinkageSIS_Pu','SIS_Pw','SIS_Pu'
%'ShrinkageKConstantCorrelation','DiagonalUnequalVariance','ShrinkageKDiagonalUnequalVariance'
%I M P O R T O   I   D A T I
%Il set A (passato) ha un numero di giorni pari a giorni A,partendo dal
%giorno iniziale giorno_iniziale, il set B ha un numero di giorni
%pari a giorniB,partendo dall'ultimo giorno del set A.






alpha=NaN;
[noth,azioni]=size(rend); 
rendimentiA=rend(giorno_iniziale:giorno_iniziale+giorniA-1,:);               
rendimentiB=rend((giorno_iniziale+giorniA):(giorno_iniziale+giorniA+giorniB-1),:);       

%C A L C O L O   I   P O R T A F O G L I   O T T I M A L I (DI MINIMA VARIANZA)
correlazioneA=corrcoef(rendimentiA);
rendimentiattesiB=mean(rendimentiB);
varB=(std(rendimentiB)).^2;

%FILTRAGGIO MATRICE DI CORRELAZIONE
switch algoritmo
	case {'Markowitz','Naive','PseudoinverseTola','PseudoinverseMarkowitz'}
        	correlazionefiltrataA=correlazioneA;
    case {'RMTStanley','RMTBouchaud','RMTCovarianceStanley','RMTCovarianceBouchaud'}
        	[V,D]=eig(correlazioneA);
        	d=diag(D);
        	sigmaquadro=1-max(real(d))/azioni;
        	Q=giorniA/azioni;
        	lmax=sigmaquadro*(1+1/Q+2*sqrt(1/Q));
        	numautograndi=length(d(find(real(d)>lmax)));
        	switch algoritmo
            		case 'RMTStanley'
                		D(find(real(D)<=lmax))=0;
                		correlazionefiltrataA=V*D*V';
                		correlazionefiltrataA=correlazionefiltrataA.*not(eye(azioni))+eye(azioni);
            		case 'RMTBouchaud' 
                		d(find(real(d)<=lmax))=mean(d(find(real(d)<=lmax)));
                		D=diag(d);
                		correlazionefiltrataA=V*D*V';
                		correlazionefiltrataA=correlazionefiltrataA./sqrt(diag(correlazionefiltrataA)*diag(correlazionefiltrataA)');
        	end
    case 'SI'
        	rendimentiattesiA=mean(rendimentiA);
        	%CARICO L'INDICE
        	load 'index';
        	rendimentiS_PA=logreturn(giorno_iniziale:(giorno_iniziale+giorniA-1),1);
        	rendimentiS_PB=logreturn((giorno_iniziale+giorniA):(giorno_iniziale+giorniA+giorniB-1),1);
        	clear logreturn
        	rendimentiattesiS_PA=mean(rendimentiS_PA);
        	varS_PA=(std(rendimentiS_PA)).^2;
        	xx=rendimentiA-repmat(rendimentiattesiA,giorniA,1);
        	xi=rendimentiS_PA-repmat(rendimentiattesiS_PA,giorniA,1);
        	covarianzaAS_PA=xx'*xi/(giorniA-1);
        	betaA=covarianzaAS_PA/varS_PA;
        	varA=(std(rendimentiA)).^2;
        	correlazionefiltrataA=varS_PA*(betaA*betaA')./(sqrt(varA'*varA));
        	correlazionefiltrataA=not(eye(azioni)).*correlazionefiltrataA+eye(azioni);
    case 'SIStanleyFirstEigenvalue'
        	[V,D]=eig(correlazioneA);
        	d=diag(D);
        	theindex=find(real(d)==max(real(d)));
        	betaA=sqrt(d(theindex))*V(:,theindex);
        	correlazionefiltrataA=betaA*betaA';
        	correlazionefiltrataA=not(eye(azioni)).*correlazionefiltrataA+eye(azioni);
    case {'MSTSingle','CLCA_c','CLCA_s','MSTAverage','MSTWeighted','MSTHausdorff_s','MSTHausdorff_h','MSTAverage_correlation','MSTWeighted_correlation'}
      		switch algoritmo
            		case 'MSTSingle'
                		method='single';
            		case {'CLCA_s','CLCA_c'}
            			method='complete';
            		case {'MSTAverage','MSTAverage_correlation'}
                		method='average';
            		case {'MSTWeighted','MSTWeighted_correlation'}
                		method='weighted';
            		case {'MSTHausdorff_s','MSTHausdorff_h'}
                		method='hausdorff';
            
            %correlazioneA(find(correlazioneA<0))=0;
            end

            switch algoritmo
                    case {'MSTAverage_correlation','MSTWeighted_correlation'}
                        tempcorrelazioneA=correlazioneA.*not(eye(azioni));
                        correlazionevec=squareform(tempcorrelazioneA);
                        Z=linkage_corr(correlazionevec,method);
                otherwise
                        distanze=sqrt(2*(1-correlazioneA));
                        distanzevec=squareform(distanze);
                        Z=linkage_P(distanzevec,method);
            end
        Zmod=Z; 
        v=[];
        for i=1:azioni-1
            	if Z(i,1)<=azioni
                	v(i,1:azioni)=[Z(i,1),zeros(1,azioni-1)];
            	else
                	v(i,1:azioni)=v(Z(i,1)-azioni,:);
            	end
            	temp1=v(i,find(v(i,:)));
            	if Z(i,2)<=azioni 
                	temp2=Z(i,2);
            	else
                	temp2=v(Z(i,2)-azioni,find(v(Z(i,2)-azioni,:)));
            	end
            	v(i,1:azioni)=[temp1,temp2,zeros(1,azioni-length(temp1)-length(temp2))];
            
            	switch algoritmo
                	case {'MSTSingle','MSTAverage','MSTWeighted','MSTHausdorff_s','CLCA_s'}
                            [ll,mm]=find(distanze==min(min(distanze(temp1,temp2))));
                	case'CLCA_c'
                    		[ll,mm]=find(distanze==max(max(distanze(temp1,temp2))));
                	case {'MSTAverage_correlation','MSTWeighted_correlation'}
                    		[ll,mm]=find(correlazioneA==max(max(correlazioneA(temp1,temp2))));
                	case 'MSTHausdorff_h'
                            [ll,mm]=find(distanze==Z(i,3));    
            	end
          
            	ll=intersect(ll,temp1);
            	mm=intersect(mm,temp2);
            	if (length(ll)==1 && length(mm)==1)
            	else 
            		warning('degenerate tree')
            	end
            	Zmod(i,1)=ll(1);
        	    Zmod(i,2)=mm(1);
        end
        G=graph;
		add(G,Zmod(:,1:2));
		H=zeros(azioni);
		for (i=1:azioni-1)
			H(Zmod(i,1),Zmod(i,2))=Zmod(i,3);
		end
		H=H+H';
																
		correlazionefiltrataA=eye(azioni);
        for i=1:azioni-1
		for j=i+1:azioni
                	P=find_path(G,i,j);
					edges=zeros(length(P)-1,1);											
					for k=1:length(P)-1											
						edges(k)=H(P(k),P(k+1));	
					end
																
                	switch algoritmo
                    		case {'MSTAverage_correlation','MSTWeighted_correlation'}
                        		correlazionefiltrataA(i,j)=min(edges); 
                    		otherwise
                        		correlazionefiltrataA(i,j)=1-((max(edges))^2)/2;
                	end
                	correlazionefiltrataA(j,i)=correlazionefiltrataA(i,j);
            	end
        end
        clear P G Z Zmod edges ll mm v temp1 temp2
    case {'ShrinkageSI','ShrinkagePerfectPositiveCorrelation','ConstantCorrelation','ShrinkageConstantCorrelation','DiagonalCommonVariance','ShrinkageDiagonalCommonVariance','ShrinkageCommonCoVariance','ShrinkageDiagonalUnequalVariance','ShrinkageSIS_Pw','ShrinkageSIS_Pu','SIS_Pw','SIS_Pu'}
    %do nothing
    case {'ShrinkageKConstantCorrelation','DiagonalUnequalVariance','ShrinkageKDiagonalUnequalVariance'}
        %matrice target
        switch algoritmo
            case 'ShrinkageKConstantCorrelation'
                r=sum(sum(triu(correlazioneA,1)))/(azioni*(azioni-1)/2);
            case {'ShrinkageKDiagonalUnequalVariance','DiagonalUnequalVariance'}
                r=0;
        end
        T=not(eye(azioni))*r+eye(azioni);
        %intensita' dello Shrinkage
        switch algoritmo
            case {'ShrinkageKConstantCorrelation','ShrinkageKDiagonalUnequalVariance'}
                alpha=(azioni+1)/(azioni+1+giorniA+((giorniA-azioni)*(giorniA-azioni-3)*trace(T*pinv(correlazioneA)*T*pinv(correlazioneA))-(giorniA-azioni-1)*trace(pinv(correlazioneA)*T)*trace(T*pinv(correlazioneA))-(2*(giorniA-azioni-1))/giorniA*trace(T*pinv(correlazioneA)))/(giorniA*azioni));
            case 'DiagonalUnequalVariance'
                alpha=1;
        end
        correlazionefiltrataA=alpha*T+(1-alpha)*correlazioneA;
    otherwise
        error('wrong argument: %s',algoritmo)
end



%FILTRAGGIO DELLA MATRICE DI COVARIANZA
switch algoritmo	
    case {'Markowitz','Naive','PseudoinverseTola','PseudoinverseMarkowitz','RMTStanley','RMTBouchaud','SI','SIStanleyFirstEigenvalue','MSTSingle','CLCA_s','CLCA_c','MSTAverage','MSTWeighted','MSTHausdorff_s','MSTHausdorff_h','MSTAverage_correlation','MSTWeighted_correlation','ShrinkageKConstantCorrelation','DiagonalUnequalVariance','ShrinkageKDiagonalUnequalVariance'}
        varA=(std(rendimentiA)).^2;
        covarianzafiltrataA=correlazionefiltrataA.*sqrt(varA'*varA);
    case {'RMTCovarianceStanley','RMTCovarianceBouchaud'}
        covarianzaA=cov(rendimentiA);
        [V,C]=eig(covarianzaA);
        dd=diag(C);
        [ddd,sortindex]=sort(real(dd));
        cc=dd(sortindex);
        ccc=cc(azioni-numautograndi);    
        if length(find(C==ccc))>1
            error('the filtered covariance matrix is not uniquely defined in this algoritm')
        else
            switch algoritmo
                case 'RMTCovarianceStanley'
                    dd(find(real(dd)<=ccc))=0;
                case 'RMTCovarianceBouchaud'
                    dd(find(real(dd)<=ccc))=mean(real(dd(find(real(dd)<=ccc))));
            end
            covarianzafiltrataA=V*diag(dd)*V';
            covarianzafiltrataA=covarianzafiltrataA.*not(eye(azioni))+diag(diag(covarianzaA));           
        end
    case {'ShrinkagePerfectPositiveCorrelation','ConstantCorrelation','ShrinkageConstantCorrelation','DiagonalCommonVariance','ShrinkageDiagonalCommonVariance','ShrinkageCommonCoVariance','ShrinkageDiagonalUnequalVariance','ShrinkageSIS_Pw','ShrinkageSIS_Pu','SIS_Pw','SIS_Pu','SIS_Pw','SIS_Pu'}
        covarianzaA=cov(rendimentiA);
        bb=1;
        switch algoritmo
            case {'ConstantCorrelation','DiagonalCommonVariance'}
                %do nothing
            otherwise
                rendimentiattesiA=mean(rendimentiA); 
                xx=rendimentiA-repmat(rendimentiattesiA,giorniA,1);
        end
               
        %matrice Target
        switch algoritmo
            case {'DiagonalCommonVariance','ShrinkageDiagonalCommonVariance'}
                r=0;%questa costante viene utilizzata successivamente
                ni=mean(diag(covarianzaA));
                Target=eye(azioni)*ni;
            case 'ShrinkageCommonCoVariance'
                r=0;%questa costante viene utilizzata successivamente
                ni=mean(diag(covarianzaA));
                csi=sum(sum(triu(covarianzaA,1))/(azioni*(azioni-1)/2));
                Target=(not(eye(azioni))*csi+eye(azioni)*ni);
            case 'ShrinkageDiagonalUnequalVariance'
                Target=diag(diag(covarianzaA));
            case 'ShrinkagePerfectPositiveCorrelation'
                r=1;%questa costante viene utilizzata successivamente
                Target=sqrt(diag(covarianzaA)*diag(covarianzaA)');
            case {'ConstantCorrelation','ShrinkageConstantCorrelation'}
                r=sum(sum(triu(correlazioneA,1)))/(azioni*(azioni-1)/2);
                Target=sqrt(diag(covarianzaA)*diag(covarianzaA)');
                Target=Target.*(r*not(eye(azioni))+eye(azioni));
            case {'ShrinkageSIS_Pw','ShrinkageSIS_Pu','SIS_Pw','SIS_Pu'}
                bb=(giorniA-1)/giorniA;
                switch algoritmo 
                    case {'ShrinkageSIS_Pw','SIS_Pw'}
                        load 'index';
                    case {'ShrinkageSIS_Pu','SIS_Pu'}
                        logreturn=mean(rend')';
                end
                rendimentiS_PA=logreturn(giorno_iniziale:(giorno_iniziale+giorniA-1),1);
                clear logreturn
                rendimentiattesiS_PA=mean(rendimentiS_PA);
                xi=rendimentiS_PA-rendimentiattesiS_PA;
                varS_PA=(std(rendimentiS_PA)).^2;
                vi=varS_PA*bb;
                covarianzaAS_PA=xx'*xi/(giorniA-1);
                ci=covarianzaAS_PA*bb;
                
                betaA=ci/vi;
                alphaA=rendimentiattesiA'-covarianzaAS_PA*rendimentiattesiS_PA/varS_PA;
                deltaA=zeros(azioni,1);
                for i=1:azioni
                    deltaA(i)=var(rendimentiA(:,i)-alphaA(i)-betaA(i)*rendimentiS_PA)*bb;
                end
                Target=vi*betaA*betaA'+diag(deltaA);
        end
        %intensita' dello Shrinkage: alpha
        switch algoritmo
            case {'ConstantCorrelation','DiagonalCommonVariance','SIS_Pw','SIS_Pu'}
                %do nothing
            otherwise
                wm=covarianzaA*(giorniA-1)/giorniA;    
                for i=1:azioni
                    for j=1:azioni
                        w(:,i,j)=xx(:,i).*xx(:,j);  
                        ww(:,i,j)=(w(:,i,j)-wm(i,j));
                    end
                end
                switch algoritmo
                    case 'ShrinkageConstantCorrelation'
                        varcovarianzaA=zeros(azioni);
                        covarcovarianzaA=[];
                         for i=1:azioni
                            for j=i:azioni
                                for l=1:azioni
                                    for m=l:azioni
                                        covarcovarianzaA(i,j,l,m)=ww(:,i,j)'*ww(:,l,m)*giorniA/(giorniA-1)^3;
                                        covarcovarianzaA(i,j,m,l)=covarcovarianzaA(i,j,l,m);
                                        covarcovarianzaA(j,i,m,l)=covarcovarianzaA(i,j,l,m);
                                        covarcovarianzaA(j,i,l,m)=covarcovarianzaA(i,j,l,m);
                                    end
                                end
                            end
                         end
                         for i=1:azioni
                            for j=i:azioni
                                varcovarianzaA(i,j)=covarcovarianzaA(i,j,i,j);
                                varcovarianzaA(j,i)=varcovarianzaA(i,j);
                            end
                         end
                         save('covarcovarianzaA','covarcovarianzaA')
                         save('varcovarianzaA','varcovarianzaA')
                    case {'ShrinkagePerfectPositiveCorrelation','ShrinkageCommonCoVariance','ShrinkageDiagonalCommonVariance','ShrinkageDiagonalUnequalVariance','ShrinkageSIS_Pw','ShrinkageSIS_Pu'}
%                        load 'varcovarianzaA'           
                         varcovarianzaA=zeros(azioni);
                         for i=1:azioni
                             for j=i:azioni 
                                 varcovarianzaA(i,j)=ww(:,i,j)'*ww(:,i,j)*giorniA/(giorniA-1)^3;
                                 varcovarianzaA(j,i)=varcovarianzaA(i,j);
                             end
                         end         
                end
        end
        switch algoritmo
            case 'ShrinkageConstantCorrelation'
                f=zeros(azioni);
                for i=1:azioni
                    for j=i+1:azioni
                        f(i,j)=(1/2)*(sqrt(covarianzaA(j,j)/covarianzaA(i,i))*covarcovarianzaA(i,i,i,j)+sqrt(covarianzaA(i,i)/covarianzaA(j,j))*covarcovarianzaA(j,j,i,j));
                        f(j,i)=f(i,j);
                    end
                    f(i,i)=1;
                end
                aa=1;
                clear covarcovarianzaA
                save('f','f')
            case {'ConstantCorrelation','DiagonalCommonVariance','SIS_Pw','SIS_Pu'}
                %do nothing
            case 'ShrinkagePerfectPositiveCorrelation'
                load 'f'
                aa=1;
            case {'ShrinkageDiagonalCommonVariance','ShrinkageCommonCoVariance'}
                f=0;
                aa=0;
            case 'ShrinkageDiagonalUnequalVariance'
                f=0;
                aa=1;
                r=0;
            case {'ShrinkageSIS_Pw','ShrinkageSIS_Pu'}
                f=zeros(azioni);
                for i=1:azioni
                    for j=i+1:azioni
                        f(i,j)=(sum(xi.*w(:,i,j).*(betaA(i)*xx(:,j)+betaA(j)*xx(:,i)-betaA(i)*betaA(j)*xi))-Target(i,j)*wm(i,j))/giorniA^2;
                        f(j,i)=f(i,j);
                    end
                    f(i,i)=varcovarianzaA(i,i)*bb^3;
                end
                aa=0;
                r=1;
                clear xi ci
        end
        switch algoritmo
            case {'DiagonalCommonVariance','ConstantCorrelation','SIS_Pw','SIS_Pu'}
                alpha=1;
            otherwise
                NN=sum(sum(varcovarianzaA*bb^3-r*f))-aa*sum(diag(varcovarianzaA*bb^3-r*f));
                DD=sum(sum((covarianzaA*bb-Target).^2));
                alpha=max(0,min(1,NN/DD));
                clear varcovarianzaA
                clear xx ww wm w
                clear f 
        end
        %filtraggio matrice covarianzaA
        covarianzafiltrataA=(alpha*Target+(1-alpha)*covarianzaA*bb)/bb;
        clear Target bb r
end

%OTTIMIZZAZIONE
switch algoritmo
    case 'Naive'
        pesi=1/azioni*ones(1,azioni);
        rischiopredetto=sqrt(pesi(1,:)*covarianzafiltrataA*pesi(1,:)');
        minimoautovalore=min(min(eig(covarianzafiltrataA)));
    otherwise

switch shortselling
    case 'N'
            rendimentiattesiA=mean(rendimentiA);
            optimization_algorithm='quadprog';
            [rischiopredetto,pesi,minimoautovalore]=minvarport(covarianzafiltrataA);
    case 'S'
        rendimentiattesiA=mean(rendimentiA);
        switch algoritmo
            case 'PseudoinverseTola'
                optimization_algorithm='tola';
            case 'PseudoinverseMarkowitz'
                optimization_algorithm='markowitz';
            otherwise
                optimization_algorithm='quadprog';
        end
        [rischiopredetto,pesi,minimoautovalore]=minvarssport(covarianzafiltrataA,optimization_algorithm);
        
    otherwise
        error('wrong argument: %s',shortselling)
end 
end %switch algoritmo
%C A L C O L O   RISCHIO REALIZZATO E RISCHIO PERCENTUALE
covarianzaB=cov(rendimentiB);

rischiorealizzato=sqrt(pesi*covarianzaB*pesi');
rischiopercentuale=(rischiorealizzato-rischiopredetto)/rischiopredetto;
end


