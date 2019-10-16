stateid3=find(stateid==1);
stateid2=find(stateid==2);
stateid4=find(stateid==3);

burstbint3r_foo=[];
cumindex_foo=[0];
for ii=1:numel(stateid3)
    burstbint3r_foo=[burstbint3r_foo;burstbint3r(cumindex(stateid3(ii))+1:cumindex(stateid3(ii)+1),:)];
    cumindex_foo=[cumindex_foo;cumindex(stateid3(ii)+1)-cumindex(stateid3(ii))];
end
cumindex_foo=cumsum(cumindex_foo);
save('D:\GitHub\FRET3cCW_GUI\Folding\trj_DA12.mat','burstbint3r_foo');
save('D:\GitHub\FRET3cCW_GUI\Folding\index_DA12.mat','cumindex_foo');

burstbint3r_foo=[];
cumindex_foo=[0];
stateid3=stateid2;
for ii=1:numel(stateid3)
    burstbint3r_foo=[burstbint3r_foo;burstbint3r(cumindex(stateid3(ii))+1:cumindex(stateid3(ii)+1),:)];
    cumindex_foo=[cumindex_foo;cumindex(stateid3(ii)+1)-cumindex(stateid3(ii))];
end
cumindex_foo=cumsum(cumindex_foo);
save('D:\GitHub\FRET3cCW_GUI\Folding\trj_DA1.mat','burstbint3r_foo');
save('D:\GitHub\FRET3cCW_GUI\Folding\index_DA1.mat','cumindex_foo');

burstbint3r_foo=[];
cumindex_foo=[0];
stateid3=stateid4;
for ii=1:numel(stateid3)
    burstbint3r_foo=[burstbint3r_foo;burstbint3r(cumindex(stateid3(ii))+1:cumindex(stateid3(ii)+1),:)];
    cumindex_foo=[cumindex_foo;cumindex(stateid3(ii)+1)-cumindex(stateid3(ii))];
end
cumindex_foo=cumsum(cumindex_foo);
save('D:\GitHub\FRET3cCW_GUI\Folding\trj_DA2.mat','burstbint3r_foo');
save('D:\GitHub\FRET3cCW_GUI\Folding\index_DA2.mat','cumindex_foo');