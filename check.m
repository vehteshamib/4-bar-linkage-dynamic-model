function E_error = check(t,z)

global L1 L2 L3 L4 m1 m2 m3 T F P g q1_0 dq1_0

q1 = z(:,1);  q2 = z(:,2);  q3 = z(:,3);  
dq1 = z(:,4);  dq2 = z(:,5);  dq3 = z(:,6);  

size = length(q1);
i=1;
E(i) = (L1^2*dq1(i)^2*m1)/6 + (L1^2*dq1(i)^2*m2)/2 + (L2^2*dq2(i)^2*m2)/6 + (L3^2*dq3(i)^2*m3)/6 - (L1*g*m1*cos(q1(i)))/2 - L1*g*m2*cos(q1(i)) - (L2*g*m2*cos(q2(i)))/2 - (L3*g*m3*cos(q3(i)))/2 + (L1*L2*dq1(i)*dq2(i)*m2*cos(q1(i) - q2(i)))/2;
for i = 2 : size
    E(i) = (L1^2*dq1(i)^2*m1)/6 + (L1^2*dq1(i)^2*m2)/2 + (L2^2*dq2(i)^2*m2)/6 + (L3^2*dq3(i)^2*m3)/6 - (L1*g*m1*cos(q1(i)))/2 - L1*g*m2*cos(q1(i)) - (L2*g*m2*cos(q2(i)))/2 - (L3*g*m3*cos(q3(i)))/2 + (L1*L2*dq1(i)*dq2(i)*m2*cos(q1(i) - q2(i)))/2;
    W(i) =T*dq1(i) + L3*P*dq3(i) - F*L1*dq1(i)*cos(q1(i) - q2(i));

    dE(i)= (E(i) - E(i-1))/(t(i) - t(i-1));
    error(i) = (dE(i) - W(i))*100;
end
error(1) = 0;
E_error = error;
