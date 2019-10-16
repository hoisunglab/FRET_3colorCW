% If isfastbinding == true, stateid == 1 (DA1A2 + DA1) is analyzed
% together by default. Labeling efficiency of A2 is a fitting parameter.
% By setting stateid = [1 3], state 3 (D-onyl + DA2) is also analyzed.
%
% If is fastbinding == false, stateid == 1 (DA1A2) is analyzed by default.
% According to stateid setting, state 2 (DA1) or state 3 (DA2) can be
% analyzed together.
%
% In the analyses of DA1, A1 and A2 are not distinguished but
% treated as a single kind. In the analyses of DA1, A1 photons are treated
% as D.
%
% fixparam:
% {[effb1 effb2], effs2, effs3, effb12}
% [effb1 effb2 effb12] are donor and accetor leak into other channels
% (effb1 = A1/(A2+A1+D), effb2 = A2/(A2+A1+D), effb12 = A2/(A1+A2)
% effs2 and effs3 are FRET efficiencies in state 2
% (A2 dark, E = (A1+A2)/(D+A1+A2)) and 3 (A1 dark, E = A2/(D+A1+A2)).
% These values are fitting parameters in some conditions

function respconv=mlhrateeiglingen3c1ExABIC(initparams,LUbounds,fixparam,frburstdata,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateid,chans)

    function logmlh=mlhratesub(params)
        pconv=diff(LUbounds,1,2)./(1+params.^2)+LUbounds(:,1);
        darkparams=pconv(end-3:end);
        pconv(end-3:end)=[];
        
        ratesum=pconv(end-2*nstate+3:end-nstate+1);
        frn=pconv(end-nstate+2:end);
        pfactor=cumprod([1; 1-frn(1:end-1)]);
        peq0=frn.*pfactor;
        peq0=[peq0; 1-sum(peq0)];
        
        effb1=fixparam{1}(1);  % E1 of A1,A2 dark state (= A1/(A2 + A1 + D))
        effb2=fixparam{1}(2);  % E2 of A1,A2 dark state (= A2/(A2 + A1 + D)) (depends on the A2-protein concentration)
        ktobrt=darkparams(1:2);  % rate from dark to bright state of A1 and A2
        frnb0=darkparams(3:4);    % A1 and A2 bright population at photon count rate of 100 ms-1
        
        effs2=fixparam{2};  % Pre-determined E_DA1 (= (A1+A2)/(D+A1+A2)) will be used unless modified below using fitting parameters
        effs3=fixparam{3};  % Pre-determined E_DA2 (= A2/(D+A1+A2)) will be used unless modified below using fitting parameters
        effb12=fixparam{4};  % E12 when no A2: A1 leak into A2 channel

        eff1s1=pconv(1:nstate);             % 3 color E1app
        eff2s1=pconv(nstate+1:2*nstate);    % 3 color E2app
        if isfastbinding,
            effs2=pconv(2*nstate+1:3*nstate-1);   % E_DA1 except unbound state
            eff12ub=eff2s1(ubid)/(eff1s1(ubid)+eff2s1(ubid));       % E12 of unbound state
            if all(stateid == [1 3]),
                effs3=pconv(3*nstate:4*nstate-1);   % E_DA2
                labelF=pconv(4*nstate);       % labeling efficiency of A2
                if isA2bleach == 1, A2brate=pconv(4*nstate+1); end     % Rate for A2 bleaching or blinking
                if isA2bleach == 2,     % For the case of no detailed balance of A2 bleaching and blinking
                    labelFapp=pconv(4*nstate+1);
                    A2brate=pconv(4*nstate+2);
                end
                effb2=effs3(ubid);  % E of A2 dark state without A1 is the same as E_DA2 of unbound state
            else
                labelF=pconv(3*nstate);       % labeling efficiency of A2
                if isA2bleach == 1, A2brate=pconv(3*nstate+1); end     % Rate for A2 bleaching or blinking
                if isA2bleach == 2,     % For the case of no detailed balance of A2 bleaching and blinking
                    labelFapp=pconv(3*nstate+1);
                    A2brate=pconv(3*nstate+2);
                end
            end
        else
            if any(stateid == 2), effs2=pconv(2*nstate+1:3*nstate); end     % E_DA1
            if any(stateid == 3), effs3=pconv(3*nstate+1:4*nstate); end     % E_DA2
        end
        
        if isfastbinding,
                fixparamnew=[effb1; effb2; effs2; effs3; effb12];
                initparamtemp=[pconv; darkparams];
                initparamtemp(3*nstate:4*nstate-1)=[];
                LUboundsdummy=[initparamtemp*0.8 initparamtemp*1.2];
            if isA2bleach < 2,
                logmlhtemp=mlhrateeiglingen3c1ExBindABI_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata1,cumindex1,indexone1,cntrate1,fixparamnew,ubid);
            else
                logmlhtemp=mlhrateeiglingen3c1ExBindnoDBABI_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata1,cumindex1,indexone1,cntrate1,fixparamnew,ubid);
            end
                logmlh=logmlhtemp(1);
                
                cntrate=zeros(length(indexone2),1);
                for k=1:length(indexone2)
                    oneburst=frburstdata2(cumindex2(indexone2(k))+1:cumindex2(indexone2(k)+1),:);
                    photoninterval=diff(oneburst(:,end-2))*1e-4; % 1 ms timeunit
                    cntrate(k)=1/mean(photoninterval);
                end
                
            ratemat0=zeros(2*nstate-1);
            for jj=1:nstate-1;
                peqtemp=peq0(jj:jj+1)/sum(peq0(jj:jj+1));
                ratemat0(jj:jj+1,jj:jj+1)=ratemat0(jj:jj+1,jj:jj+1)+ratesum(jj)*[-peqtemp(2) peqtemp(1);peqtemp(2) -peqtemp(1)];
            end
            ratemat1=ratemat0(1:nstate,1:nstate);
            ratemat0(:,ubid)=ratemat0(:,ubid)*labelF;
            ratemat1(:,ubid)=ratemat1(:,ubid)*(1-labelF);
            A2darkStid=nstate+1:2*nstate;
            A2darkStid(ubid+1:end)=A2darkStid(ubid+1:end)-1;
            A2darkStid(ubid)=ubid;
            ratemat0(A2darkStid,A2darkStid)=ratemat0(A2darkStid,A2darkStid)+ratemat1;
            peqfb=[peq0*labelF; peq0*(1-labelF)];
            peqfb(ubid)=peqfb(ubid)+peqfb(ubid+nstate);
            peqfb(nstate+ubid)=[];
            if isA2bleach,
                brateaddid=1:nstate;
                brateaddid(ubid)=[];
                for jj=brateaddid;
                    if isA2bleach < 2,
                        ratemat0([jj A2darkStid(jj)],[jj A2darkStid(jj)])=ratemat0([jj A2darkStid(jj)],[jj A2darkStid(jj)])+A2brate*[-(1-labelF) labelF; 1-labelF -labelF];
                    else
                        ratemat0([jj A2darkStid(jj)],[jj A2darkStid(jj)])=ratemat0([jj A2darkStid(jj)],[jj A2darkStid(jj)])+A2brate*[-(1-labelFapp) labelFapp; 1-labelFapp -labelFapp];
                    end
                end
            end
                
            effs3New=[effs3; repmat(effs3(ubid),size(effs3,1)-1,1)];
            
            pinput=[peqfb peqfb ones(size(peqfb))];     % First column is not used. Second is equilibrium population and the third is summation over states in the likelihood calculation

            inputparam=[effs3New; effb2; ktobrt(2); frnb0(2)];
            LUboundsdummy=[min(inputparam,[],2)*0.8 max(inputparam,[],2)*1.2]; % dummy LUbounds
            logmlhtemp=mlhTPgenABIcal([inputparam inputparam],LUboundsdummy,frburstdata2,cumindex2,indexone2,cntrate,ratemat0,pinput);
        
                logmlh=logmlh+logmlhtemp(1);                
        else
            % 3-color segment likelihood
            fixparamnew=[effb1; effb2; effs2; effs3; effb12];
            initparamtemp=[pconv(1:2*nstate); ratesum; frn; darkparams];
            LUboundsdummy=[initparamtemp*0.8 initparamtemp*1.2];
            logmlhtemp=mlhrateeiglingen3c1ExABI_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata1,cumindex1,indexone1,cntrate1,fixparamnew);
            logmlh=logmlhtemp(1);

            if any(stateid == 2),       % 2-color (DA1) segment likelihood
                cntrate=zeros(length(indexone2),1);
                for k=1:length(indexone2)
                    oneburst=frburstdata2(cumindex2(indexone2(k))+1:cumindex2(indexone2(k)+1),:);
                    photoninterval=diff(oneburst(:,end-2))*1e-4; % 1 ms timeunit
                    cntrate(k)=1/mean(photoninterval);
                end
                % For the case of the definition, E2 = A2/(D+A2)
%                 initparamtemp=[effs2; 1-(1-effb1)*(1-effb2)/(1-effb1*effb2); ratesum; frn; ktobrt(1); frnb0(1)];
                % For the case of the definition, E2 = A2/(D+A1+A2)
                initparamtemp=[effs2; effb1; ratesum; frn; ktobrt(1); frnb0(1)];                
                LUboundsdummy=[initparamtemp*0.8 initparamtemp*1.2]; % dummy LUbounds
                logmlhtemp=mlhrateeiglingenABI_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata2,cumindex2,indexone2,cntrate);
                logmlh=logmlh+logmlhtemp(1);
            end
            if any(stateid == 3),       % 2-color (DA2) segment likelihood
                cntrate=zeros(length(indexone3),1);
                for k=1:length(indexone3)
                    oneburst=frburstdata3(cumindex3(indexone3(k))+1:cumindex3(indexone3(k)+1),:);
                    photoninterval=diff(oneburst(:,end-2))*1e-4; % 1 ms timeunit
                    cntrate(k)=1/mean(photoninterval);
                end
                initparamtemp=[effs3; effb2; ratesum; frn; ktobrt(2); frnb0(2)];
                LUboundsdummy=[initparamtemp*0.8 initparamtemp*1.2]; % dummy LUbounds
                logmlhtemp=mlhrateeiglingenABI_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata3,cumindex3,indexone3,cntrate);
                logmlh=logmlh+logmlhtemp(1);
            end
        end
    end

    donchan=chans(3);
    ac1chan=chans(2);
    ac2chan=chans(1);

    frburstdata1=frburstdata{1};
    cumindex1=cumindex{1};
    indexone1=indexone{1};
    donid=find(frburstdata1(:,end) == donchan);
    ac1id=find(frburstdata1(:,end) == ac1chan);
    frburstdata1(:,end)=1;
    frburstdata1(ac1id,end)=2;
    frburstdata1(donid,end)=3;  % Convert colors for 3-color analysis
    
    fixparam1=[fixparam{1}; fixparam{2}; fixparam{3}; fixparam{4}];
    cntrate1=zeros(length(indexone1),1);
    for kk=1:length(indexone1);
        oneburst1=frburstdata1(cumindex1(indexone1(kk))+1:cumindex1(indexone1(kk)+1),:);
        photoninterval1=diff(oneburst1(:,end-2))*1e-4;        % in ms
        cntrate1(kk)=1/mean(photoninterval1);
    end

    if isfastbinding,
        if all (stateid == [1 3]),
            nstate=(length(initparams)-2-isA2bleach)/6;
            frburstdata2=frburstdata{2};
            cumindex2=cumindex{2};
            indexone2=indexone{2};
            ac2id=find(frburstdata2(:,end) == ac2chan);
            frburstdata2(:,end)=2;
            frburstdata2(ac2id,end)=1;  % Convert colors for 2-color analysis
        else
            nstate=(length(initparams)-2-isA2bleach)/5;
            if isA2bleach < 2,
                respconv=mlhrateeiglingen3c1ExBindABI_MT(initparams,LUbounds,frburstdata1,cumindex1,indexone1,cntrate1,fixparam1,ubid);
            else
                respconv=mlhrateeiglingen3c1ExBindnoDBABI_MT(initparams,LUbounds,frburstdata1,cumindex1,indexone1,cntrate1,fixparam1,ubid);
            end
        end
    else
        nstate=(length(initparams)-2)/(length(stateid)+3);
        if all(stateid == 1),
            respconv=mlhrateeiglingen3c1ExABI_MT(initparams,LUbounds,frburstdata1,cumindex1,indexone1,cntrate1,fixparam1);
        end
        if any(stateid == 2),
            frburstdata2=frburstdata{2};
            cumindex2=cumindex{2};
            indexone2=indexone{2};
            donid=find(frburstdata2(:,end) == donchan);
            frburstdata2(:,end)=1;
            frburstdata2(donid,end)=2;  % Convert colors for 2-color analysis
        end            
        if any(stateid == 3),
            frburstdata3=frburstdata{3};
            cumindex3=cumindex{3};
            indexone3=indexone{3};
            ac2id=find(frburstdata3(:,end) == ac2chan);
            frburstdata3(:,end)=2;
            frburstdata3(ac2id,end)=1;  % Convert colors for 2-color analysis
        end        
    end

    if length(stateid) > 1,
        invinitparams=sqrt(diff(LUbounds,1,2)./(initparams-LUbounds(:,1))-1);
        options=optimset('MaxFunEval', 40000, 'MaxIter', 2000);
        resparams=fminsearch(@mlhratesub,invinitparams,options);
        respconv=diff(LUbounds,1,2)./(1+resparams.^2)+LUbounds(:,1);
    end

end