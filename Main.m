load data
d1 = data.quark7;
d1(17:end,:)=[];
d{1} = d1(:,1);
d{2} = d1(:,2);
d{3} = d1(:,2);ones(length(d1(:,1)),1);

pt = linspace(0,1,20).';
[param,sp,err]=MainFit(d,pt)

save lhc