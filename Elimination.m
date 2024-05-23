function dz = Elimination (t,z,t_init)

global L1 L2 L3 L4 m1 m2 m3 T F P g
q1 = z(1);  q2 = z(2);  q3 = z(3);
dq1= z(4);  dq2= z(5);  dq3= z(6);
q=[q1;q2;q3];
dq=[dq1;dq2;dq3];

M =[      (L1^2*(m1 + 3*m2))/3, (L1*L2*m2*cos(q1 - q2))/2,           0
     (L1*L2*m2*cos(q1 - q2))/2,               (L2^2*m2)/3,           0
                             0,                         0, (L3^2*m3)/3];
B =[T - F*L1*cos(q1 - q2) - (L1*g*m1*sin(q1))/2 - L1*g*m2*sin(q1) - (L1*L2*dq2^2*m2*sin(q1 - q2))/2
                                                   -(L2*m2*(- L1*sin(q1 - q2)*dq1^2 + g*sin(q2)))/2
                                                                        (L3*(2*P - g*m3*sin(q3)))/2];a =[ L1*cos(q1), L2*cos(q2), -L3*cos(q3)
     L1*sin(q1), L2*sin(q2), -L3*sin(q3)];
% a_dot =[ -L1*dq1*sin(q1), -L2*dq2*sin(q2),  L3*dq3*sin(q3)
%           L1*dq1*cos(q1),  L2*dq2*cos(q2), -L3*dq3*cos(q3)];
a_dot_dq =[- L1*sin(q1)*dq1^2 - L2*sin(q2)*dq2^2 + L3*sin(q3)*dq3^2
             L1*cos(q1)*dq1^2 + L2*cos(q2)*dq2^2 - L3*cos(q3)*dq3^2];
         
MM =[M(3,:)-a(:,3)'*(a(:,1:2)'\M(1:2,:))
       -a];
BB =[B(3)-a(:,3)'*(a(:,1:2)'\B(1:2))
    a_dot_dq];
X=MM\BB;
dz=[dq;X];
