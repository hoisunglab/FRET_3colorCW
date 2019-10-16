function [errorparams, logmlh0, bic]=mlhrateeiglingen3cAlexABICerror(initparams,LUbounds,frburstdata,cumindex,indexone,fixed, stateid,exid,chans,denatConc,drvintv,drvpnt)


minValue=1E-15; %replace 0 to minValue

if isempty(denatConc) %no denaturant conc. information
    denatConc=indexone*0+1;  %set to 1 (no dependency on denaturant conc.)
end

initparamsMat=initparams;%cell2mat(initparams');
LUboundsMat=cell2mat(LUbounds');
fixedMat=cell2mat(fixed');
%fixedMat=fixedMat*0+1;
%fixedMat(1:3)=0;
%assert(size(initparamsMat,1)==size(LUboundsMat,1) && size(LUboundsMat,1)==size(fixedMat,1)); %Row lengths conincide

donchan=chans(2);
ac1chan=chans(1);
ac2chan=chans(3);

donid=find(frburstdata(:,end) == donchan);
ac1id=find(frburstdata(:,end) == ac1chan);
frburstdata(:,end)=1;
frburstdata(ac1id,end)=2;
frburstdata(donid,end)=3;  % Convert colors for 3-color analysis

donBackId=find(frburstdata(:,end)==3 & exid==2);
frburstdata(donBackId,:)=[];
exid(donBackId)=[];
%we need 2 countrate for each excitation
cntrate1=zeros(length(indexone),2); %acceptor 1,2 countrate for exc1
cntrate2=zeros(length(indexone),2); %acceptor 1,2 countrate for exc2

for iione=1:length(indexone)
    oneburst=frburstdata(cumindex(indexone(iione))+1:cumindex(indexone(iione)+1),:);
    exidone=exid(cumindex(indexone(iione))+1:cumindex(indexone(iione)+1),:);
    oneburst1=oneburst(exidone==1,:);
    oneburst2=oneburst(exidone==2,:);
    photoninterval1=diff(oneburst1(:,end-2))*1e-4; % 1 ms timeunit
    photoninterval2=diff(oneburst2(:,end-2))*1e-4; % 1 ms timeunit
    totalCntrate1=1/mean(photoninterval1);
    totalCntrate2=1/mean(photoninterval2);
    
    tot1Exc=size(oneburst1,1);
    num1Acc1=numel(find(oneburst1(:,4)==2));
    num1Acc2=numel(find(oneburst1(:,4)==1));
    
    tot2Exc=size(oneburst2,1);    
    num2Acc1=numel(find(oneburst2(:,4)==2));
    num2Acc2=numel(find(oneburst2(:,4)==1));
    
    cntrate1(iione,1)=totalCntrate1/tot1Exc*num1Acc1; %exc1 acc1
    cntrate1(iione,2)=totalCntrate1/tot1Exc*num1Acc2; %exc1 acc2
    
    cntrate2(iione,1)=totalCntrate2/tot2Exc*num2Acc1; %exc2 acc1
    cntrate2(iione,2)=totalCntrate2/tot2Exc*num2Acc2; %exc2 acc2
    
    
end
cntrate1(isnan(cntrate1))=0;
cntrate2(isnan(cntrate2))=0;

cntrate=mean([sum(cntrate1,2) sum(cntrate2,2)],2); %average acceptor countrate

initparamsNew=initparams;
initparamsNew(initparamsNew<=0)=minValue;

LUboundsNew=cell2mat(LUbounds');
fixedNew=fixedMat;%cell2mat(fixed');

%fixedNew=fixedNew*0;
initparams=initparamsNew;
drvintvs=1+drvintv*(-drvpnt:drvpnt);
initparammat=[];
for kk=1:length(initparams)
    if fixedNew(kk)==1
        continue;
    end
    drvintvsOne=drvintvs;
    %shirinkg drvintvs
    %while drvintvsOne(1)*initparams(kk) <= LUboundsNew(kk,1) || drvintvsOne(end)*initparams(kk) >= LUboundsNew(kk,2)
    %    drvintvsOne=1+(drvintvsOne-1)/2;
    %end
    if drvintvsOne(1)*initparams(kk) <= LUboundsNew(kk,1) || drvintvsOne(end)*initparams(kk) >= LUboundsNew(kk,2)
        fixedNew(kk)=1; %out of the boundary
        continue;
    end
    for nn=1:length(drvintvsOne)
        initparamone=initparams;
        initparamone(kk)=initparams(kk)*drvintvsOne(nn);
        initparammat=[initparammat initparamone];
    end
end

for kk=1:length(initparams)
    if fixedNew(kk)==1
        continue;
    end
    for mm=kk+1:length(initparams)
        if fixedNew(mm)==1
            continue;
        end
        drvintvsOne=drvintvs;
        %shirink drvintvs
        %while drvintvsOne(1)*initparams(kk) <= LUboundsNew(kk,1) || drvintvsOne(end)*initparams(kk) >= LUboundsNew(kk,2) ...
        %        || drvintvsOne(1)*initparams(mm) <= LUboundsNew(mm,1) || drvintvsOne(end)*initparams(mm) >= LUboundsNew(mm,2)
        %    drvintvsOne=1+(drvintvsOne-1)/2;
        %end
        for nn=1:length(drvintvsOne)
            initparamone=initparams;
            initparamone(kk)=initparams(kk)*drvintvsOne(nn);
            initparamone(mm)=initparams(mm)*drvintvsOne(nn);
            initparammat=[initparammat initparamone];
        end
    end
end

LUboundsNew=[min(initparammat,[],2)*0.9 max(initparammat,[],2)*1.1]; %just for a case, put a broader LUbounds
if numel(initparammat)<1
    fprintf(2,'Nothing to evaluate in mn2ABIerrorC\n');
    return;
end
logmlh=mlhrateeiglingen3cAlexABI_MT(initparammat,LUboundsNew,frburstdata,cumindex,indexone,cntrate,fixedMat*0,stateid,exid,denatConc);
Hessianmat=zeros(length(initparams)-sum(fixedNew),length(initparams)-sum(fixedNew));
logmlhdiag=zeros(1,length(drvintvs));

mlhIndex=1;
for kk=1:length(initparams)
    if fixedNew(kk)==1
        continue;
    end
    logmlhdiag=logmlhdiag*0;
    for nn=1:length(drvintvs)
        logmlhdiag(nn)=logmlh(mlhIndex);
        mlhIndex=mlhIndex+1;
    end
    drvintvsOne=drvintvs;
    %while drvintvsOne(1)*initparams(kk) <= LUboundsNew(kk,1) || drvintvsOne(end)*initparams(kk) >= LUboundsNew(kk,2)
    %    drvintvsOne=1+(drvintvsOne-1)/2;
    %end
    diagcoeff=polyfit(drvintvsOne,logmlhdiag,2);
    Hessianmat(kk-sum(fixedNew(1:kk)),kk-sum(fixedNew(1:kk)))=diagcoeff(1)/initparams(kk)^2;
end

logmlhoffdiag=zeros(1,length(drvintvs));
for kk=1:length(initparams)
    if fixedNew(kk)==1
        continue;
    end
    for mm=kk+1:length(initparams)
        if fixedNew(mm)==1
            continue;
        end
        logmlhoffdiag=logmlhoffdiag*0;
        
        drvintvsOne=drvintvs;
        %while drvintvsOne(1)*initparams(kk) <= LUboundsNew(kk,1) || drvintvsOne(end)*initparams(kk) >= LUboundsNew(kk,2) ...
        %        || drvintvsOne(1)*initparams(mm) <= LUboundsNew(mm,1) || drvintvsOne(end)*initparams(mm) >= LUboundsNew(mm,2)
        %    drvintvsOne=1+(drvintvsOne-1)/2;
        %end
        for nn=1:length(drvintvs)
            logmlhoffdiag(nn)=logmlh(mlhIndex);
            mlhIndex=mlhIndex+1;
        end
        offdiagcoeff=polyfit(drvintvsOne,logmlhoffdiag,2);
        
        kkNew=kk-sum(fixedNew(1:kk));
        mmNew=mm-sum(fixedNew(1:mm));
        Hessianmat(kkNew,mmNew)=(offdiagcoeff(1)-Hessianmat(kkNew,kkNew)*initparams(kk)^2-Hessianmat(mmNew,mmNew)*initparams(mm)^2)/2/initparams(kk)/initparams(mm);
        Hessianmat(mmNew,kkNew)=Hessianmat(kkNew,mmNew);
    end
end

covmat=inv(2*Hessianmat);
errorparamsNew=sqrt(diag(covmat));
errorparams=zeros(numel(initparams),1);
errorparams(logical(1-fixedNew))=errorparamsNew;

% Change parameter
% errorparams(nState+i)=
logmlh0=-logmlh(drvpnt+1); %logmlh for the found parameter
numphotons=sum(cumindex(indexone+1)-cumindex(indexone));
bic=-2*logmlh0 + length(errorparamsNew)*log(numphotons);

end