% If isfastbinding == true, stateid == 1 (DA1A2 + DA1) is analyzed
% together by default. Labeling efficiency of A2 is a fitting parameter.
% By setting stateid = [1 3], state 3 (D-onyl + DA2) is also analyzed.
%
% If fastbinding == false, stateid == 1 (DA1A2) is analyzed by default.
% According to stateid setting, state 2 (DA1) or state 3 (DA2) can be
% analyzed together.
%
% In the analyses of DA1 or DA2, A1 and A2 are not distinguished but
% treated as a single kind. If only A1 (or A2) photons are need to be
% analyzed, A2 (or A1) photons should be removed prior to a function call.

function respconv=mlhrateeiglingen3c1ExC(initparams,LUbounds,frburstdata,cumindex,indexone,isfastbinding,isA2bleach,ubid,stateid,chans)

    function logmlh=mlhratesub(params)
        pconv=diff(LUbounds,1,2)./(1+params.^2)+LUbounds(:,1);
        
        ratesum=pconv(end-2*nstate+3:end-nstate+1);
        frn=pconv(end-nstate+2:end);
        pfactor=cumprod([1; 1-frn(1:end-1)]);
        peq=frn.*pfactor;
        peq=[peq; 1-sum(peq)];
        
%        eff1s1=pconv(1:nstate);             % 3 color E1app
%        eff2s1=pconv(nstate+1:2*nstate);    % 3 color E2app
        if isfastbinding,
%            effs2=pconv(2*nstate+1:3*nstate-1);   % E_DA1 except unbound state
%            eff12ub=eff2s1(ubid)/(eff1s1(ubid)+eff2s1(ubid));       % E12 of unbound state
            effs3=pconv(3*nstate:4*nstate-1);   % E_DA2
            labelF=pconv(4*nstate);       % labeling efficiency of A2
            if isA2bleach == 1, A2brate=pconv(4*nstate+1); end     % Rate for A2 bleaching or blinking
            if isA2bleach == 2,     % For the case of no detailed balance of A2 bleaching and blinking
                labelFapp=pconv(4*nstate+1);
                A2brate=pconv(4*nstate+2);
            end
        else
            if any(stateid == 2), effs2=pconv(2*nstate+1:3*nstate); end     % E_DA1
            if any(stateid == 3), effs3=pconv(3*nstate+1:4*nstate); end     % E_DA2
        end
        
        if isfastbinding,
                initparamtemp=pconv;
                initparamtemp(3*nstate:4*nstate-1)=[];
                LUboundsdummy=[initparamtemp*0.8 initparamtemp*1.2];
            if isA2bleach < 2,
                logmlhtemp=mlhrateeiglingen3c1ExBind_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata1,cumindex1,indexone1,ubid);
            else
                logmlhtemp=mlhrateeiglingen3c1ExBindnoDB_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata1,cumindex1,indexone1,ubid);
            end                
                logmlh=logmlhtemp(1);
                
            ratemat0=zeros(2*nstate-1);
            for jj=1:nstate-1;
                peqtemp=peq(jj:jj+1)/sum(peq(jj:jj+1));
                ratemat0(jj:jj+1,jj:jj+1)=ratemat0(jj:jj+1,jj:jj+1)+ratesum(jj)*[-peqtemp(2) peqtemp(1);peqtemp(2) -peqtemp(1)];
            end
            ratemat1=ratemat0(1:nstate,1:nstate);
            ratemat0(:,ubid)=ratemat0(:,ubid)*labelF;
            ratemat1(:,ubid)=ratemat1(:,ubid)*(1-labelF);
            A2darkStid=nstate+1:2*nstate;
            A2darkStid(ubid+1:end)=A2darkStid(ubid+1:end)-1;
            A2darkStid(ubid)=ubid;
            ratemat0(A2darkStid,A2darkStid)=ratemat0(A2darkStid,A2darkStid)+ratemat1;
            peqfb=[peq*labelF; peq*(1-labelF)];
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

            inputparam=effs3New;
            LUboundsdummy=[min(inputparam,[],2)*0.8 max(inputparam,[],2)*1.2]; % dummy LUbounds
            logmlhtemp=mlhTPgencal([inputparam inputparam],LUboundsdummy,frburstdata2,cumindex2,indexone2,ratemat0,pinput);
            
                logmlh=logmlh+logmlhtemp(1);
        else %if isfastbinding~=1
            initparamtemp=[pconv(1:2*nstate); ratesum; frn];
            LUboundsdummy=[initparamtemp*0.8 initparamtemp*1.2];
            logmlhtemp=mlhrateeiglingen3c1Ex_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata1,cumindex1,indexone1);
            logmlh=logmlhtemp(1);
            
            if any(stateid == 2),       % 2-color (DA1) segment likelihood
                initparamtemp=[effs2; ratesum; frn];
                LUboundsdummy=[initparamtemp*0.8 initparamtemp*1.2]; % dummy LUbounds
                logmlhtemp=mlhrateeiglingen_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata2,cumindex2,indexone2);
                logmlh=logmlh+logmlhtemp(1);
            end
            if any(stateid == 3),       % 2-color (DA2) segment likelihood
                initparamtemp=[effs3; ratesum; frn];
                LUboundsdummy=[initparamtemp*0.8 initparamtemp*1.2]; % dummy LUbounds
                logmlhtemp=mlhrateeiglingen_MT([initparamtemp initparamtemp],LUboundsdummy,frburstdata3,cumindex3,indexone3);
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
    
    if isfastbinding,
        if all (stateid == [1 3]),
            nstate=(length(initparams)+2-isA2bleach)/6;
            frburstdata2=frburstdata{2};
            cumindex2=cumindex{2};
            indexone2=indexone{2};
            ac2id=find(frburstdata2(:,end) == ac2chan);
            frburstdata2(:,end)=2;
            frburstdata2(ac2id,end)=1;  % Convert colors for 2-color analysis
        else
            nstate=(length(initparams)+2-isA2bleach)/5;
            if isA2bleach < 2,
                respconv=mlhrateeiglingen3c1ExBind_MT(initparams,LUbounds,frburstdata1,cumindex1,indexone1,ubid);
            else
                respconv=mlhrateeiglingen3c1ExBindnoDB_MT(initparams,LUbounds,frburstdata1,cumindex1,indexone1,ubid);
            end
        end
    else
        nstate=(length(initparams)+2)/(length(stateid)+3);
        if all(stateid == 1),
            respconv=mlhrateeiglingen3c1Ex_MT(initparams,LUbounds,frburstdata1,cumindex1,indexone1);
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