function z_init = init ()

global q1_0 dq1_0
fun1 = @root2d1;
x0 = [0,0,0,0];
x = fsolve(fun1,x0);
z_init =[q1_0 x(1:2) dq1_0 x(3:4)];