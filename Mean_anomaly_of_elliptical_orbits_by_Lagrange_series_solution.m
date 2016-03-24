%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg


% Mean anomaly of elliptical orbits by Lagrange series solution
close all; clear all; clc;
%% inputs
Me=linspace(0,2*pi,10);
e=0.65;
N=[3,10];
%% calculations
% exact solution
Me2=Me-e*sin(Me);
% Lagrange series solution
for m=1:length(N)
    A{m}=zeros(N(m),length(Me));
    E1{m}=zeros(N(m),length(Me));
    for n=1:N(m)
        for k=0:floor(n/2)
            A{m}(n,:)=A{m}(n,:)+(-1)^k/factorial(n-k)/factorial(k)*(n-2*k)^(n-1)*sin((n-2*k)*Me);
        end
        a(n,:)=(A{m}(n,:)).*0.5.^(n-1);
        E1{m}(n,:)=E1{m}(n,:)+a(n,:)*e^n;
    end
end
% [cc,bb]=size(E1);
for mm=1:length(N)
    [cc(mm),bb(mm)]=size(E1{mm});
    EE{mm}=zeros(1,bb(mm));
    for kk=1:cc(mm)
        EE{mm}=EE{mm}+E1{mm}(kk,1:length(Me));
    end
    EEE{mm}=EE{mm}+Me;
end
%% plotting
figure(1);
set(gcf,'color','w');
hold all;
% exact solution
plot(Me,Me2,'linewidth',2);
legendInfo1{1} = ['Exact soluion'];
% Lagrange series solution
for nn1=1:length(N)
    plot(EEE{nn1},Me,'linewidth',2);
    legendInfo1{nn1+1} = ['N = ' num2str(N(nn1))];
end
xlabel('E','fontsize',18);
ylabel('M_e','fontsize',18);
title('Mean anomaly vs True anomaly with Lagrange series solution','fontsize',18);
legend(legendInfo1,'location','northwest');
grid on;
xlim([0,2*pi]);
ylim([0,2*pi]);
%---------------------------------------------------------------------------------------------------------------------------------------------------------
% error
figure(2);
set(gcf,'color','w');
hold all;
for b1=1:length(N)
    plot(Me,(Me-EEE{b1}),'linewidth',2);
    legendInfo1{b1} = ['Error of N = ' num2str(N(b1))];
end
xlabel('Me','fontsize',18);
ylabel('Exact-Bessel solution','fontsize',18);
title('Error in Eccentric anomaly  E','fontsize',18);
legend(legendInfo1,'location','northwest');
grid on;
