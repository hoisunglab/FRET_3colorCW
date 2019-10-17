%initparams{1}=[eff1s1 ... eff1sn ... eff2s1 ... eff2sn] for DA1A2
%initparams{2}=[eff1s1 ... eff1sn leakage] for DA1
%initparams{3}=[eff2s1 ... eff2sn] for DA2
%initparams{4}=[eff12s1 ... eff12sn] for A1A2 (acceptor excitation)
%initparams{5}=[eff1D eff2D] for D-only
%initparams{6}=[eff2D gamma] for A2 only (A1 excitation).
%gamma=A2onlyCountrate/A1A2countrate
%initparams{7}=[k(1) ~ k(n-1)]  k(1) is the sum of foward and backward rates between state 1 and 2;
%initparams{8}=[f(1) ~ f(n-1)] p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
%Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
%initparams{9}=[kb1 pb1 kb2 pb2] blinking description for A1 and A2.

function respconv=mlhrateeiglingen3cAlexABIC(initparams,LUbounds,frburstdata,cumindex,indexone,fixed, stateid,exid,chans,denatConc)

if ~exist('denatConc','var') %no denaturant conc. information
    denatConc=indexone*0+1;  %set to 1 (no dependency on denaturant conc.)
end


initparamsMat=cell2mat(initparams');
LUboundsMat=cell2mat(LUbounds');
fixedMat=cell2mat(fixed');
assert(size(initparamsMat,1)==size(LUboundsMat,1) && size(LUboundsMat,1)==size(fixedMat,1)); %Row lengths conincide

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

numSb=(numel(initparamsMat)-7)/7;
if any(denatConc~=1) && numSb~=2
    fprintf(2,'Concentration dependecy is for 2-state only\n');
    assert(numSb==2);
    return;
end

respconv=mlhrateeiglingen3cAlexABI_MT(initparamsMat,LUboundsMat,frburstdata,cumindex,indexone,cntrate,fixedMat,stateid,exid,denatConc);
end