%% Analyze roaming, dwelling, and quiescence from activity data captured with WorMotel
% Matt Churgin 
% Part 1a: Fit individual worms to extract roaming and dwelling information
clear all
close all
load ExampleData

ActVal=ActVal(:,1:59000);
onel=5;
ActVal=conv2(double(ActVal(:,1:59000)),ones(1,onel)/onel,'same');
nbins=[0:25:2000];
for i=1:48 %number of worms
    temp=ActVal(i,:);
    temp(temp==0)=[]; % remove quiescent frames from fit analysis
    [ah ax]=hist(temp,nbins); % generate histogram
    [af{i} gof1{i}]=createMattFit2(ax,ah); % create fit
    disp(['fit ' num2str(i)])
end

save ExampleFits
%% Part 1b: Plot individual worms fit versus actual data
close all

i=1; % example region of interest to plot

temp=ActVal(i,:);
temp(temp==0)=[];
[ah ax]=hist(temp,nbins);

x=ax;
% Fit parameters
a=af{i}.a;
b=af{i}.b;
c=af{i}.c;
d=af{i}.d;
e=af{i}.e;
f=af{i}.f;
g=af{i}.g;

plot(ax,ah,'k')
hold on
plot(ax,a*exp(-b*x)+c*exp(-d*(x-e).^2)+f*exp(-g*x),'ro')
legend('Data','Fit')
legend boxoff
box off
xlabel('Activity (A.U.)')
ylabel('Number of Frames')
set(gca,'FontSize',15)

%% Part 2: Calculate Quiescence
clear all
close all
load ExampleData

% designate genotype groups based on region of interest
a=[1:4 (1:4)+8];
b=a+16;
c=a+32;
d=a+4;
e=b+4;
f=c+4;

onel=5;
ActSmooth=conv2(double(ActVal(:,1:59000)),ones(1,onel)/onel,'same');
q=sum(ActSmooth<1,2)/3600;
qa=mean(q(a)); qas=std(q(a))/sqrt(length(a));
qb=mean(q(b)); qbs=std(q(b))/sqrt(length(b));
qc=mean(q(c)); qcs=std(q(c))/sqrt(length(c));
qd=mean(q(d)); qds=std(q(d))/sqrt(length(d));
qe=mean(q(e)); qes=std(q(e))/sqrt(length(e));
qf=mean(q(f)); qfs=std(q(f))/sqrt(length(f));


hold on
plot(1+(rand(length(a),1)-0.5)*0.1,q(a),'k.','MarkerSize',10)
plot(2+(rand(length(b),1)-0.5)*0.1,q(b),'k.','MarkerSize',10)
plot(3+(rand(length(c),1)-0.5)*0.1,q(c),'k.','MarkerSize',10)
plot(4+(rand(length(d),1)-0.5)*0.1,q(d),'k.','MarkerSize',10)
plot(5+(rand(length(e),1)-0.5)*0.1,q(e),'k.','MarkerSize',10)
plot(6+(rand(length(f),1)-0.5)*0.1,q(f),'k.','MarkerSize',10)

ylabel('Quiescence (Hours)')
xlabel('Genotype')
set(gca,'FontSize',15)
xlim([0 7])

qa=q(a);
qb=q(b);
qc=q(c);
qd=q(d);
qe=q(e);
qf=q(f);

save ExampleQuiescence qa qb qc qd qe qf a b c d e f
%% Part 3: Calculate roaming, dwelling, and quiescence fractions and plot
clear all
close all

load ExampleFits
load ExampleQuiescence

G1=[];
G2=[];
G3=[];
G4=[];
G5=[];
G6=[];
for i=1:length(a)
    G1=[G1; af{a(i)}.a af{a(i)}.b af{a(i)}.c af{a(i)}.d af{a(i)}.e af{a(i)}.f af{a(i)}.g qa(i)];
    G2=[G2; af{d(i)}.a af{d(i)}.b af{d(i)}.c af{d(i)}.d af{d(i)}.e af{d(i)}.f af{d(i)}.g qd(i)];
    G3=[G3; af{b(i)}.a af{b(i)}.b af{b(i)}.c af{b(i)}.d af{b(i)}.e af{b(i)}.f af{b(i)}.g qb(i)];
    G4=[G4; af{e(i)}.a af{e(i)}.b af{e(i)}.c af{e(i)}.d af{e(i)}.e af{e(i)}.f af{e(i)}.g qe(i)];
    G5=[G5; af{c(i)}.a af{c(i)}.b af{c(i)}.c af{c(i)}.d af{c(i)}.e af{c(i)}.f af{c(i)}.g qc(i)];
    G6=[G6; af{f(i)}.a af{f(i)}.b af{f(i)}.c af{f(i)}.d af{f(i)}.e af{f(i)}.f af{f(i)}.g qf(i)];
end
% convert data to area under curve
divis=1;
x=[0:25:2000];
for i=1:length(G1)
    G1E(i)=G1(i,1).*sum(exp(-G1(i,2).*x))+G1(i,6).*sum(exp(-G1(i,7).*x));
    G1G(i)=G1(i,3).*sum(exp(-G1(i,4).*(x-G1(i,5)).^2))/divis;
    
    tot=1-G1(i,8)/16;
    tot2=G1E(i)+G1G(i);
    G1E(i)=G1E(i)/tot2*tot;
    G1G(i)=G1G(i)/tot2*tot;
end

for i=1:length(G2)
    G2E(i)=G2(i,1).*sum(exp(-G2(i,2).*x))+G2(i,6).*sum(exp(-G2(i,7).*x));
    G2G(i)=G2(i,3).*sum(exp(-G2(i,4).*(x-G2(i,5)).^2))/divis;
    
    tot=1-G2(i,8)/16;
    tot2=G2E(i)+G2G(i);
    G2E(i)=G2E(i)/tot2*tot;
    G2G(i)=G2G(i)/tot2*tot;
end

for i=1:length(G3)
    G3E(i)=G3(i,1).*sum(exp(-G3(i,2).*x))+G3(i,6).*sum(exp(-G3(i,7).*x));
    G3G(i)=G3(i,3).*sum(exp(-G3(i,4).*(x-G3(i,5)).^2))/divis;
    
    tot=1-G3(i,8)/16;
    tot2=G3E(i)+G3G(i);
    G3E(i)=G3E(i)/tot2*tot;
    G3G(i)=G3G(i)/tot2*tot;
end

for i=1:length(G4)
    G4E(i)=G4(i,1).*sum(exp(-G4(i,2).*x))+G4(i,6).*sum(exp(-G4(i,7).*x));
    G4G(i)=G4(i,3).*sum(exp(-G4(i,4).*(x-G4(i,5)).^2))/divis;
    
    tot=1-G4(i,8)/16;
    tot2=G4E(i)+G4G(i);
    G4E(i)=G4E(i)/tot2*tot;
    G4G(i)=G4G(i)/tot2*tot;
end

for i=1:length(G5)
    G5E(i)=G5(i,1).*sum(exp(-G5(i,2).*x))+G5(i,6).*sum(exp(-G5(i,7).*x));
    G5G(i)=G5(i,3).*sum(exp(-G5(i,4).*(x-G5(i,5)).^2))/divis;
    
    tot=1-G5(i,8)/16;
    tot2=G5E(i)+G5G(i);
    G5E(i)=G5E(i)/tot2*tot;
    G5G(i)=G5G(i)/tot2*tot;
end
for i=1:length(G6)
    G6E(i)=G6(i,1).*sum(exp(-G6(i,2).*x))+G6(i,6).*sum(exp(-G6(i,7).*x));
    G6G(i)=G6(i,3).*sum(exp(-G6(i,4).*(x-G6(i,5)).^2))/divis;
    
    tot=1-G6(i,8)/16;
    tot2=G6E(i)+G6G(i);
    G6E(i)=G6E(i)/tot2*tot;
    G6G(i)=G6G(i)/tot2*tot;
end
rois{1}=[G1E; G1G];
rois{2}=[G2E; G2G];
rois{3}=[G3E; G3G];
rois{4}=[G4E; G4G];
rois{5}=[G5E; G5G];
rois{6}=[G6E; G6G];


C(:,1)=[0 0 0];
C(:,2)=[0.1 0.8 0.1];
C(:,3)=[0.8 0.1 0.2];
C(:,4)=[0.1 0.3 0.8];
C(:,5)=[0.8 0.3 0.8];
C(:,6)=[0.35 0.65 0.85];
C(:,7)=[0.65 0.65 0.65];
C(:,8)=[0.85 0.85 0.85];
figure
hold on
et=0:0.25:2*pi;

for z=1:6
    
    [myPC{z} myScore{z}]=princomp([rois{z}(1,:)' rois{z}(2,:)']);
    theta=atan(myPC{z}(1,1)/myPC{z}(2,1));
   
    % use standard deviation or standard error for ellipse size
    a=std(myScore{z}(:,2))/sqrt(length(a));
    b=std(myScore{z}(:,1))/sqrt(length(a));
    
    ex{z}=mean(rois{z}(1,:))-a*cos(et)*cos(theta)+b*sin(et)*sin(theta);
    ey{z}=mean(rois{z}(2,:))+a*cos(et)*sin(theta)+b*sin(et)*cos(theta);

        plot(ex{z},ey{z},'Color',C(:,z),'LineWidth',3)

    
end
x=0:16;
xlabel('Dwelling Fraction')
ylabel('Roaming Fraction')
legend('Condition 1','Condition 2','Condition 3','Condition 4','Condition 5','Condition 6')
legend boxoff
axis([0 6 0 16])
axis([0 0.75 0.0 .95])
set(gca,'FontSize',15)
set(gca,'YTick',[0:.2:1])