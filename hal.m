clc; clear all;

global L1 L2 L3 L4 m1 m2 m3 T F P g q1_0 dq1_0

m1=1; m2=2; m3=3; L1=0.9; L2=2; L3=1.1; L4=2; g=9.81; T=1; F=2; P=1.5;
q1_0=0; dq1_0=-1;
option = odeset('maxstep',0.001);
z0 = init()

% ode
for j=1:1
tic
[t,z] = ode45(@Augmented,[0:0.01:5],z0,option);
toc
q1_Aug = z(:,1);  q2_Aug = z(:,2);  q3_Aug = z(:,3); 
dq1_Aug = z(:,4); dq2_Aug = z(:,5); dq3_Aug = z(:,6);
E_error_Aug = check(t,z);

tic
[t,z] = ode45(@Elimination,[0:0.01:5],z0,option);
toc
q1_Eli = z(:,1);  q2_Eli = z(:,2);  q3_Eli = z(:,3); 
dq1_Eli = z(:,4); dq2_Eli = z(:,5); dq3_Eli = z(:,6);
E_error_Eli = check(t,z);

tic
[t,z] = ode45(@Greenwood,[0:0.01:5],z0,option);
toc
q1_Gre = z(:,1);  q2_Gre = z(:,2);  q3_Gre = z(:,3); 
dq1_Gre = z(:,4); dq2_Gre = z(:,5); dq3_Gre = z(:,6);
E_error_Gre = check(t,z);

tic
[t,z] = ode45(@Integrated_Multiplier,[0:0.01:5],[z0,0,0],option);
toc
q1_IMM = z(:,1);  q2_IMM = z(:,2);  q3_IMM = z(:,3); 
dq1_IMM = z(:,4); dq2_IMM = z(:,5); dq3_IMM = z(:,6);
E_error_IMM = check(t,z);

tic
[t,z] = ode45(@Embedding,[0:0.01:5],[z0(1:3),z0(4)],option);
toc
q1_Emb = z(:,1);  q2_Emb = z(:,2);  q3_Emb = z(:,3); 
u_Emb = z(:,4);
size = length(t);
for i=1:size
    W =[                                    1
     -(L1*sin(q1_Emb(i) - q3_Emb(i)))/(L2*sin(q2_Emb(i) - q3_Emb(i)))
     -(L1*sin(q1_Emb(i) - q2_Emb(i)))/(L3*sin(q2_Emb(i) - q3_Emb(i)))];
 dq1_Emb(i) = W(1)*u_Emb(i);
 dq2_Emb(i) = W(2)*u_Emb(i);
 dq3_Emb(i) = W(3)*u_Emb(i);
end
E_error_Emb = check(t,[z(:,1:3) dq1_Emb' dq2_Emb' dq3_Emb']);
tic
[t,z] = ode45(@Modified_Lagrange,[0:0.01:5],[z0(1:3),-z0(4)+z0(5)],option);
toc
q1_ML = z(:,1);  q2_ML = z(:,2);  q3_ML = z(:,3); 
u_ML = z(:,4);
size = length(t);
for i=1:size
    W =[     -(L2*sin(q2_ML(i) - q3_ML(i)))/(L1*sin(q1_ML(i) - q3_ML(i)) + L2*sin(q2_ML(i) - q3_ML(i)))
              (L1*sin(q1_ML(i) - q3_ML(i)))/(L1*sin(q1_ML(i) - q3_ML(i)) + L2*sin(q2_ML(i) - q3_ML(i)))
     (L1*L2*sin(q1_ML(i) - q2_ML(i)))/(L1*L3*sin(q1_ML(i) - q3_ML(i)) + L2*L3*sin(q2_ML(i) - q3_ML(i)))];

 dq1_ML(i) = W(1)*u_ML(i);
 dq2_ML(i) = W(2)*u_ML(i);
 dq3_ML(i) = W(3)*u_ML(i);
end
E_error_ML = check(t,[z(:,1:3) dq1_ML' dq2_ML' dq3_ML']);
end

% q(j) plots 
for i=1:1
figure
hold on
grid on
plot(t,q1_Aug,'r','linewidth',10)
plot(t,q1_Eli,'g--','linewidth',10)
plot(t,q1_Gre,'b','linewidth',6)
plot(t,q1_IMM,'k--','linewidth',6)
plot(t,q1_Emb,'y','linewidth',2)
plot(t,q1_ML,'c--','linewidth',2)
ylabel('q1','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('Augmented','Elimination','Greenwood','Integrated Multiplier','Embedding','Modified Lagrange') 

figure
hold on
grid on
plot(t,q2_Aug,'r','linewidth',10)
plot(t,q2_Eli,'g--','linewidth',10)
plot(t,q2_Gre,'b','linewidth',6)
plot(t,q2_IMM,'k--','linewidth',6)
plot(t,q2_Emb,'y','linewidth',2)
plot(t,q2_ML,'c--','linewidth',2)
ylabel('q2','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('Augmented','Elimination','Greenwood','Integrated Multiplier','Embedding','Modified Lagrange') 

figure
hold on
grid on
plot(t,q3_Aug,'r','linewidth',10)
plot(t,q3_Eli,'g--','linewidth',10)
plot(t,q3_Gre,'b','linewidth',6)
plot(t,q3_IMM,'k--','linewidth',6)
plot(t,q3_Emb,'y','linewidth',2)
plot(t,q3_ML,'c--','linewidth',2)
ylabel('q3','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('Augmented','Elimination','Greenwood','Integrated Multiplier','Embedding','Modified Lagrange') 
end

% constant errors
for i=1:length(z)
   Cx_Aug(i) = (L1*sin(q1_Aug(i)) - L4 + L2*sin(q2_Aug(i)) - L3*sin(q3_Aug(i)));
   Cy_Aug(i) = (L3*cos(q3_Aug(i)) - L2*cos(q2_Aug(i)) - L1*cos(q1_Aug(i)));
   
   Cx_Eli(i) = (1*sin(q1_Eli(i)) - L4 + L2*sin(q2_Eli(i)) - L3*sin(q3_Eli(i)));
   Cy_Eli(i) = (L3*cos(q3_Eli(i)) - L2*cos(q2_Eli(i)) - L1*cos(q1_Eli(i)));
   
   Cx_Gre(i) = (L1*sin(q1_Gre(i)) - L4 + L2*sin(q2_Gre(i)) - L3*sin(q3_Gre(i)));
   Cy_Gre(i) = (L3*cos(q3_Gre(i)) - L2*cos(q2_Gre(i)) - L1*cos(q1_Gre(i)));
   
   Cx_IMM(i) = (L1*sin(q1_IMM(i)) - L4 + L2*sin(q2_IMM(i)) - L3*sin(q3_IMM(i)));
   Cy_IMM(i) = (L3*cos(q3_IMM(i)) - L2*cos(q2_IMM(i)) - L1*cos(q1_IMM(i)));
   
   Cx_Emb(i) = (L1*sin(q1_Emb(i)) - L4 + L2*sin(q2_Emb(i)) - L3*sin(q3_Emb(i)));
   Cy_Emb(i) = (L3*cos(q3_Emb(i)) - L2*cos(q2_Emb(i)) - L1*cos(q1_Emb(i)));
   
   Cx_ML(i) = (L1*sin(q1_ML(i)) - L4 + L2*sin(q2_ML(i)) - L3*sin(q3_ML(i)));
   Cy_ML(i) = (L3*cos(q3_ML(i)) - L2*cos(q2_ML(i)) - L1*cos(q1_ML(i)));
end

% error plots
for i=1:1
figure
hold on
grid on
plot(t,Cx_Aug,'r','linewidth',10)
plot(t,Cx_Eli,'g--','linewidth',10)
plot(t,Cx_Gre,'b','linewidth',6)
plot(t,Cx_IMM,'k--','linewidth',6)
plot(t,Cx_Emb,'y','linewidth',2)
plot(t,Cx_ML,'c--','linewidth',2)
ylabel('x error','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('Augmented','Elimination','Greenwood','Integrated Multiplier','Embedding','Modified Lagrange') 

figure
hold on
grid on
plot(t,Cy_Aug,'r','linewidth',10)
plot(t,Cy_Eli,'g--','linewidth',10)
plot(t,Cy_Gre,'b','linewidth',6)
plot(t,Cy_IMM,'k--','linewidth',6)
plot(t,Cy_Emb,'y','linewidth',2)
plot(t,Cy_ML,'c--','linewidth',2)
ylabel('y error','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('Augmented','Elimination','Greenwood','Integrated Multiplier','Embedding','Modified Lagrange') 

figure
hold on
grid on
plot(t,E_error_Aug,'r','linewidth',10)
plot(t,E_error_Eli,'g--','linewidth',10)
plot(t,E_error_Gre,'b','linewidth',6)
plot(t,E_error_IMM,'k--','linewidth',6)
plot(t,E_error_Emb,'y','linewidth',2)
plot(t,E_error_ML,'c--','linewidth',2)
ylabel('E error','fontsize',18)
xlabel('Time (s)','fontsize',18)
legend('Augmented','Elimination','Greenwood','Integrated Multiplier','Embedding','Modified Lagrange') 
end
