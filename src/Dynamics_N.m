
function [G4] = Dynamics_N(q,dq,ddq)

%Recursive Newton-Euler for standard Denavit-Hartenberg notation

th1=q(1);
th2=q(2);
th3=q(3);
th4=q(4);
th5=q(5);
th6=q(6);
th7=q(7);

% standard DH Parameters of the robot
d = [0.2703 0 0.3644 0 0.3743 0 0.2295];
a = [0.069 0 0.069 0 0.01 0 0];
th = [th1 th2 th3 th4 th5 th6 th7]; 
al = [-pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2 0];


% Generalized transformation matrix for a link CS
% Tr = [cos(th) -sin(th)*cos(al)  sin(th)*sin(al) a*cos(th); ...
%       sin(th)  cos(th)*cos(al) -cos(th)*sin(al) a*sin(th); ...
%             0          sin(al)          cos(al)         d; ...
%             0               0                0          1]

%Forward kinematic transformation  matrices for individual link,compute the link rotation matrices


A1 = [cos(th(1))  -sin(th(1))*cos(al(1))  sin(th(1))*sin(al(1)) a(1)*cos(th(1)); ...
      sin(th(1))   cos(th(1))*cos(al(1)) -cos(th(1))*sin(al(1)) a(1)*sin(th(1)); ...
               0              sin(al(1))             cos(al(1))            d(1); ...
               0                       0                      0              1];
A2 = [cos(th(2))  -sin(th(2))*cos(al(2))  sin(th(2))*sin(al(2)) a(2)*cos(th(2)); ...
      sin(th(2))   cos(th(2))*cos(al(2)) -cos(th(2))*sin(al(2)) a(2)*sin(th(2)); ...
               0              sin(al(2))             cos(al(2))            d(2); ...
               0                       0                      0              1];
A3 = [cos(th(3))  -sin(th(3))*cos(al(3))  sin(th(3))*sin(al(3)) a(3)*cos(th(3)); ...
      sin(th(3))   cos(th(3))*cos(al(3)) -cos(th(3))*sin(al(3)) a(3)*sin(th(3)); ...
               0              sin(al(3))             cos(al(3))            d(3); ...
               0                       0                      0              1];
A4 = [cos(th(4))  -sin(th(4))*cos(al(4))  sin(th(4))*sin(al(4)) a(4)*cos(th(4)); ...
      sin(th(4))   cos(th(4))*cos(al(4)) -cos(th(4))*sin(al(4)) a(4)*sin(th(4)); ...
               0              sin(al(4))             cos(al(4))            d(4); ...
               0                       0                      0              1];
A5 = [cos(th(5))  -sin(th(5))*cos(al(5))  sin(th(5))*sin(al(5)) a(5)*cos(th(5)); ...
      sin(th(5))   cos(th(5))*cos(al(5)) -cos(th(5))*sin(al(5)) a(5)*sin(th(5)); ...
               0              sin(al(5))             cos(al(5))            d(5); ...
               0                       0                      0              1];
A6 = [cos(th(6))  -sin(th(6))*cos(al(6))  sin(th(6))*sin(al(6)) a(6)*cos(th(6)); ...
      sin(th(6))   cos(th(6))*cos(al(6)) -cos(th(6))*sin(al(6)) a(6)*sin(th(6)); ...
               0              sin(al(6))             cos(al(6))            d(6); ...
               0                       0                      0              1];
           
A7 = [cos(th(7))  -sin(th(7))*cos(al(7))  sin(th(7))*sin(al(7)) a(7)*cos(th(7)); ...
      sin(th(7))   cos(th(7))*cos(al(7)) -cos(th(7))*sin(al(7)) a(7)*sin(th(7)); ...
               0              sin(al(7))             cos(al(7))            d(7); ...
               0                       0                      0              1];

R(:,:,1)=A1(1:3,1:3);
R(:,:,2)=A2(1:3,1:3);
R(:,:,3)=A3(1:3,1:3);
R(:,:,4)=A4(1:3,1:3);
R(:,:,5)=A5(1:3,1:3);
R(:,:,6)=A6(1:3,1:3);
R(:,:,7)=A7(1:3,1:3);

P(:,1)=A1(1:3,4);
P(:,2)=A2(1:3,4);
P(:,3)=A3(1:3,4);
P(:,4)=A4(1:3,4);
P(:,5)=A5(1:3,4);
P(:,6)=A6(1:3,4);
P(:,7)=A7(1:3,4);
%init some variables
g0=[0 0 -9.80665]';
z0=[0 0 1]';
m1=5.70044;m2=3.22698;m3=4.31272;m4=2.07206;m5=2.24665;m6=1.60979;m7=0.54218;
m=[m1 m2 m3 m4 m5 m6 m7];


x(:,1)=[-0.05117 0.07908 0.00086]';
x(:,2)=[0.00269 -0.00529 0.06845]';
x(:,3)=[-0.07176 0.08149 0.00132]';
x(:,4)=[0.00159 -0.01117 0.02618]';
x(:,5)=[-0.01168 0.13111 0.00460]';
x(:,6)=[0.00697  0.00600 0.06048]';
x(:,7)=[0.005137 0.0009572 -0.06682]';

% fprintf('x1: '); disp(x(:,1))
%惯量
Ixx1=0.0470910226;     Iyy1=0.035959884;    Izz1=0.0376697645;
Ixx2=0.027885975;      Iyy2=0.020787492;    Izz2=0.0117520941;
Ixx3=0.0266173355;     Iyy3=0.012480083;    Izz3=0.0284435520;
Ixx4=0.0131822787;     Iyy4=0.009268520;    Izz4=0.0071158268;
Ixx5=0.0166774282;     Iyy5=0.003746311;    Izz5=0.0167545726;
Ixx6=0.0070053791;     Iyy6=0.005527552;    Izz6=0.0038760715;
Ixx7=0.0008162135;     Iyy7=0.0008735012;   Izz7=0.0005494148;
%惯量
Ixy1=-0.0061487003;    Ixz1= 0.0001278755 ;  Iyz1=-0.0007808689;
Ixy2=-0.0001882199;    Ixz2=-0.00030096397;  Iyz2= 0.0020767576;
Ixy3=-0.0039218988;    Ixz3= 0.0002927063 ;  Iyz3=-0.001083893 ;
Ixy4=-0.0001966341;    Ixz4= 0.0003603617 ;  Iyz4= 0.000745949 ; 
Ixy5=-0.0001865762;    Ixz5= 0.0001840370 ;  Iyz5= 0.0006473235;
Ixy6= 0.0001534806;    Ixz6=-0.0004438478 ;  Iyz6=-0.0002111503;
Ixy7= 0.000128440 ;    Ixz7= 0.00018969891;  Iyz7= 0.0001057726;

%inertia tensor of augmented link

I(:,:,1)=[ Ixx1 -Ixy1 -Ixz1
          -Ixy1  Iyy1 -Iyz1
          -Ixz1 -Iyz1  Izz1];                            
                                  
I(:,:,2)=[ Ixx2 -Ixy2 -Ixz2
          -Ixy2  Iyy2 -Iyz2
          -Ixz2 -Iyz2  Izz2]; 
    
I(:,:,3)=[ Ixx3 -Ixy3 -Ixz3
          -Ixy3  Iyy3 -Iyz3
          -Ixz3 -Iyz3  Izz3];  
    
I(:,:,4)=[ Ixx4 -Ixy4 -Ixz4
          -Ixy4  Iyy4 -Iyz4
          -Ixz4 -Iyz4  Izz4];      
    
I(:,:,5)=[ Ixx5 -Ixy5 -Ixz5
          -Ixy5  Iyy5 -Iyz5
          -Ixz5 -Iyz5  Izz5];     
    
I(:,:,6)=[ Ixx6 -Ixy6 -Ixz6
          -Ixy6  Iyy6 -Iyz6
          -Ixz6 -Iyz6  Izz6];

I(:,:,7)=[ Ixx7 -Ixy7 -Ixz7
          -Ixy7  Iyy7 -Iyz7
          -Ixz7 -Iyz7  Izz7];  
          
              
%the forward recursion 1-7
w=zeros(3,1);%angular velocity of link
dw=zeros(3,1);%angular acceleration of link
ddp=-g0;%linear acceleration of origin of Frame i
ddpc=zeros(3,7);
nn=zeros(3,7);
f=zeros(3,8);
u=zeros(3,8);
tou=zeros(7,1);
%statement order is important here
%revolute axis

for i=1:7
    dw=transpose(R(:,:,i))*(dw+z0*ddq(i)+cross(w,z0*dq(i)));
    w=transpose(R(:,:,i))*(w+z0*dq(i));
%     fprintf('w: '); disp(w)
%     fprintf('dw: '); disp(dw)
    ddp=transpose(R(:,:,i))*ddp+cross(dw,transpose(R(:,:,i))*P(:,i))+cross(w,cross(w,transpose(R(:,:,i))*P(:,i)));
    ddpc(:,i)=ddp+cross(dw,x(:,i))+cross(w,cross(w,x(:,i)));  
    nn(:,i)=I(:,:,i)*dw+cross(w,I(:,:,i)*w);
end


%order of these statements is important, since both
%u and f are functions of previous f.
for j=7:-1:1
    if j==7
        T=eye(3,3);
    else
        T=R(:,:,j+1);
    end
    f(:,j)=T*f(:,j+1)+m(j)*ddpc(:,j);
    u(:,j)=nn(:,j)+T*u(:,j+1)+cross(T*f(:,j+1),x(:,j))+cross(x(:,j)+transpose(R(:,:,j))*P(:,j),f(:,j));
    
    Q=R(:,:,j);
    tou(j,1)=u(:,j)'*(transpose(Q)*z0);   
end
G4=tou;




































































































