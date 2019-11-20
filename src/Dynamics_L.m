
%Recursive Lagrange for standard Denavit-Hartenberg notation

function [G3] = Dynamics_L(q,dq,ddq)


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
th = [th1 th2 th3 th4 th5 th6,th7]; 
al = [-pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2 0];


% Generalized transformation matrix for a link CS
% Tr = [cos(th) -sin(th)*cos(al)  sin(th)*sin(al) a*cos(th); ...
%       sin(th)  cos(th)*cos(al) -cos(th)*sin(al) a*sin(th); ...
%             0  sin(al)                  cos(al)         d; ...
%             0        0                       0          1]

%Forward kinematic transformation  matrices for individual link

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

           
g=[0 0 -9.80665 0];
Q= [0 -1 0 0
    1  0 0 0
    0  0 0 0
    0  0 0 0];
m1=5.70044;m2=3.22698;m3=4.31272;m4=2.07206;m5=2.24665;m6=1.60979;m7=0.54218;
x1=[-0.05117 0.07908 0.00086 1];%Column vector
x2=[0.00269 -0.00529 0.06845 1];
x3=[-0.07176 0.08149 0.00132 1];
x4=[0.00159 -0.01117 0.02618 1];
x5=[-0.01168 0.13111 0.00460 1];
x6=[0.00697  0.00600 0.06048 1];
x7=[0.005137 0.0009572 -0.06682 1];

%link Inertia tensors (all units kg·m2)
% link      Ixx            Iyy           Izz
% 1     0.0470910226   0.035959884   0.0376697645
% 2     0.027885975    0.020787492   0.0117520941
% 3     0.0266173355   0.012480083   0.0284435520
% 4     0.0131822787   0.009268520   0.0071158268
% 5     0.0166774282   0.003746311   0.0167545726
% 6     0.0070053791   0.005527552   0.0038760715
% 7     0.0008162135   0.0008735012  0.0005494148

%link        Ixy            Iyz             Ixz
% 1     -0.0061487003   -0.0007808689    0.0001278755
% 2     -0.0001882199    0.0020767576   -0.00030096397
% 3     -0.0039218988   -0.001083893     0.0002927063
% 4     -0.0001966341    0.000745949     0.0003603617
% 5     -0.0001865762    0.0006473235    0.0001840370
% 6      0.0001534806   -0.0002111503   -0.0004438478
% 7      0.000128440     0.0001057726    0.00018969891

%calculate Pseudo-inertia matrix of each link
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

%(xi,yi,zi) is center of mass
%I=[(-Ixx+Iyy+Izz)/2          Ixy             Ixz         mixi
%         Ixy           (Ixx-Iyy+Izz)/2       Iyz         miyi
%         Ixz                 Iyz        (Ixx+Iyy-Izz)/2  mizi
%         mixi                miyi           mizi         mi  ];

I1=[(-Ixx1+Iyy1+Izz1)/2       Ixy1             Ixz1           m1*x1(1)
          Ixy1         (Ixx1-Iyy1+Izz1)/2      Iyz1           m1*x1(2)
          Ixz1                Iyz1        (Ixx1+Iyy1-Izz1)/2  m1*x1(3)
        m1*x1(1)            m1*x1(2)          m1*x1(3)             m1];                            
                                  
I2=[(-Ixx2+Iyy2+Izz2)/2       Ixy2             Ixz2           m2*x2(1)
          Ixy2         (Ixx2-Iyy2+Izz2)/2      Iyz2           m2*x2(2)
          Ixz2                Iyz2        (Ixx2+Iyy2-Izz2)/2  m2*x2(3)
        m2*x2(1)            m2*x2(2)          m2*x2(3)             m2]; 
    
I3=[(-Ixx3+Iyy3+Izz3)/2       Ixy3             Ixz3           m3*x3(1)
          Ixy3         (Ixx3-Iyy3+Izz3)/2      Iyz3           m3*x3(2)
          Ixz3                Iyz3        (Ixx3+Iyy3-Izz3)/2  m3*x3(3)
        m3*x3(1)            m3*x3(2)           m3*x3(3)            m3];  
    
I4=[(-Ixx4+Iyy4+Izz4)/2       Ixy4             Ixz4           m4*x4(1)
          Ixy4         (Ixx4-Iyy4+Izz4)/2      Iyz4           m4*x4(2)
          Ixz4                Iyz4        (Ixx4+Iyy4-Izz4)/2  m4*x4(3)
        m4*x4(1)            m4*x4(2)           m4*x4(3)            m4];      
    
I5=[(-Ixx5+Iyy5+Izz5)/2       Ixy5             Ixz5           m5*x5(1)
          Ixy5         (Ixx5-Iyy5+Izz5)/2      Iyz5           m5*x5(2)
          Ixz5                Iyz5        (Ixx5+Iyy5-Izz5)/2  m5*x5(3)
        m5*x5(1)            m5*x5(2)           m5*x5(3)            m5];     
    
I6=[(-Ixx6+Iyy6+Izz6)/2       Ixy6             Ixz6           m6*x6(1)
          Ixy6         (Ixx6-Iyy6+Izz6)/2      Iyz6           m6*x6(2)
          Ixz6                Iyz6        (Ixx6+Iyy6-Izz6)/2  m6*x6(3)
        m6*x6(1)            m6*x6(2)           m6*x6(3)            m6];     
    
I7=[(-Ixx7+Iyy7+Izz7)/2       Ixy7             Ixz7           m7*x7(1)
          Ixy7         (Ixx7-Iyy7+Izz7)/2      Iyz7           m7*x7(2)
          Ixz7                Iyz7        (Ixx7+Iyy7-Izz7)/2  m7*x7(3)
        m7*x7(1)            m7*x7(2)           m7*x7(3)            m7];      
    
    

%calculate gravity term
T10=A1;T20=A1*A2;T30=A1*A2*A3;T40=A1*A2*A3*A4;T50=A1*A2*A3*A4*A5;T60=A1*A2*A3*A4*A5*A6;T70=A1*A2*A3*A4*A5*A6*A7;
T21=A2;T31=A2*A3;T41=A2*A3*A4;T51=A2*A3*A4*A5;T61=A2*A3*A4*A5*A6;T71=A2*A3*A4*A5*A6*A7;
T32=A3;T42=A3*A4;T52=A3*A4*A5;T62=A3*A4*A5*A6;T72=A3*A4*A5*A6*A7;
T43=A4;T53=A4*A5;T63=A4*A5*A6;T73=A4*A5*A6*A7;
T54=A5;T64=A5*A6;T74=A5*A6*A7;
T65=A6;T75=A6*A7;
T76=A7;

u11=Q*T10;      u21=Q*T20;      u31=Q*T30;     u41=Q*T40;      u51=Q*T50;     u61=Q*T60;     u71=Q*T60;
u22=T10*Q*T21;  u32=T10*Q*T31;  u42=T10*Q*T41; u52=T10*Q*T51;  u62=T10*Q*T61; u72=T10*Q*T71;
u33=T20*Q*T32;  u43=T20*Q*T42;  u53=T20*Q*T52; u63=T20*Q*T62;  u73=T20*Q*T72;
u44=T30*Q*T43;  u54=T30*Q*T53;  u64=T30*Q*T63; u74=T30*Q*T73;
u55=T40*Q*T54;  u65=T40*Q*T64;  u75=T40*Q*T74;
u66=T50*Q*T65;  u76=T50*Q*T75;
u77=T60*Q*T76;

g1=-(m1*g*u11*x1'+m2*g*u21*x2'+m3*g*u31*x3'+m4*g*u41*x4'+m5*g*u51*x5'+m6*g*u61*x6'+m7*g*u71*x7');
g2=-(m2*g*u22*x2'+m3*g*u32*x3'+m4*g*u42*x4'+m5*g*u52*x5'+m6*g*u62*x6'+m7*g*u72*x7');
g3=-(m3*g*u33*x3'+m4*g*u43*x4'+m5*g*u53*x5'+m6*g*u63*x6'+m7*g*u73*x7');
g4=-(m4*g*u44*x4'+m5*g*u54*x5'+m6*g*u64*x6'+m7*g*u74*x7');
g5=-(m5*g*u55*x5'+m6*g*u65*x6'+m7*g*u75*x7');
g6=-(m6*g*u66*x6'+m7*g*u76*x7');
g7=-(m7*g*u77*x7');

%与加速度相关的惯量矩阵项
%calculate inertia term D(q),symmetric, bounded, positive definite
% D=[D11 D12 D13 D14 D15 D16 D17
%    D21 D22 D23 D24 D25 D26 D27
%    D31 D32 D33 D34 D35 D36 D37
%    D41 D42 D43 D44 D45 D46 D47 
%    D51 D52 D53 D54 D55 D56 D57
%    D61 D62 D63 D64 D65 D66 D67 
%    D71 D72 D73 D74 D75 D76 D77];

D11=trace(u11*I1*transpose(u11))+trace(u21*I2*transpose(u21))+trace(u31*I3*transpose(u31))+...
    trace(u41*I4*transpose(u41))+trace(u51*I5*transpose(u51))+trace(u61*I6*transpose(u61))+trace(u71*I7*transpose(u71));
D12=trace(u22*I2*transpose(u21))+trace(u32*I3*transpose(u31))+trace(u42*I4*transpose(u41))+...
    trace(u52*I5*transpose(u51))+trace(u62*I6*transpose(u61))+trace(u72*I7*transpose(u71));
D13=trace(u33*I3*transpose(u31))+trace(u43*I4*transpose(u41))+trace(u53*I5*transpose(u51))+...
    trace(u63*I6*transpose(u61))+trace(u73*I7*transpose(u71));      
D14=trace(u44*I4*transpose(u41))+trace(u54*I5*transpose(u51))+trace(u64*I6*transpose(u61))+...
    trace(u74*I7*transpose(u71));   
D15=trace(u55*I5*transpose(u51))+trace(u65*I6*transpose(u61))+trace(u75*I7*transpose(u71));
D16=trace(u66*I6*transpose(u61))+trace(u76*I7*transpose(u71));
D17=trace(u77*I7*transpose(u71));

D22=trace(u22*I2*transpose(u22))+trace(u32*I3*transpose(u32))+trace(u42*I4*transpose(u42))+...
    trace(u52*I5*transpose(u52))+trace(u62*I6*transpose(u62))+trace(u72*I7*transpose(u72));
D23=trace(u33*I3*transpose(u32))+trace(u43*I4*transpose(u42))+trace(u53*I5*transpose(u52))+...
    trace(u63*I6*transpose(u62))+trace(u73*I7*transpose(u72));
D24=trace(u44*I4*transpose(u42))+trace(u54*I5*transpose(u52))+trace(u64*I6*transpose(u62))+...
    trace(u74*I7*transpose(u72));
D25=trace(u55*I5*transpose(u52))+trace(u65*I6*transpose(u62))+trace(u75*I7*transpose(u72));
D26=trace(u66*I6*transpose(u62))+trace(u76*I7*transpose(u72));
D27=trace(u77*I7*transpose(u72));

D33=trace(u33*I3*transpose(u33))+trace(u43*I4*transpose(u43))+trace(u53*I5*transpose(u53))+...
    trace(u63*I6*transpose(u63))+trace(u73*I7*transpose(u73));
D34=trace(u44*I4*transpose(u43))+trace(u54*I5*transpose(u53))+...
    trace(u64*I6*transpose(u63))+trace(u74*I7*transpose(u73));
D35=trace(u55*I5*transpose(u53))+trace(u65*I6*transpose(u63))+trace(u75*I7*transpose(u73));
D36=trace(u66*I6*transpose(u63))+trace(u76*I7*transpose(u73));
D37=trace(u77*I7*transpose(u73));

D44=trace(u44*I4*transpose(u44))+trace(u54*I5*transpose(u54))+...
    trace(u64*I6*transpose(u64))+trace(u74*I7*transpose(u74));
D45=trace(u55*I5*transpose(u54))+trace(u65*I6*transpose(u64))+...
    trace(u75*I7*transpose(u74));
D46=trace(u66*I6*transpose(u64))+trace(u76*I7*transpose(u74));
D47=trace(u77*I7*transpose(u74));

D55=trace(u55*I5*transpose(u55))+trace(u65*I6*transpose(u65))+...
    trace(u75*I7*transpose(u75));
D56=trace(u66*I6*transpose(u65))+trace(u76*I7*transpose(u75));
D57=trace(u77*I7*transpose(u75));

D66=trace(u66*I6*transpose(u66))+trace(u76*I7*transpose(u76));
D67=trace(u77*I7*transpose(u76));

D77=trace(u77*I7*transpose(u77));


D21=D12;
D31=D13;D32=D23;
D41=D14;D42=D24;D43=D34;
D51=D15;D52=D25;D53=D35;D54=D45;
D61=D16;D62=D26;D63=D36;D64=D46;D65=D56;
D71=D17;D72=D27;D73=D37;D74=D47;D75=D57;D76=D67;

D1=D11*ddq(1)+D12*ddq(2)+D13*ddq(3)+D14*ddq(4)+D15*ddq(5)+D16*ddq(6)+D17*ddq(7);
D2=D21*ddq(1)+D22*ddq(2)+D23*ddq(3)+D24*ddq(4)+D25*ddq(5)+D26*ddq(6)+D27*ddq(7);
D3=D31*ddq(1)+D32*ddq(2)+D33*ddq(3)+D34*ddq(4)+D35*ddq(5)+D36*ddq(6)+D37*ddq(7);
D4=D41*ddq(1)+D42*ddq(2)+D43*ddq(3)+D44*ddq(4)+D45*ddq(5)+D46*ddq(6)+D47*ddq(7);
D5=D51*ddq(1)+D52*ddq(2)+D53*ddq(3)+D54*ddq(4)+D55*ddq(5)+D56*ddq(6)+D57*ddq(7);
D6=D61*ddq(1)+D62*ddq(2)+D63*ddq(3)+D64*ddq(4)+D65*ddq(5)+D66*ddq(6)+D67*ddq(7);
D7=D71*ddq(1)+D72*ddq(2)+D73*ddq(3)+D74*ddq(4)+D75*ddq(5)+D76*ddq(6)+D77*ddq(7);

%calculate Coriolis and Centrifugal term h(q,dot(q))
%下面的矩阵是方阵,并且是对称矩�?,�?以只�?要推导一半数量即�?

u111=Q^2*T10;

u211=Q^2*T20;    u212=Q*T10*Q*T21;
u221=Q*T10*Q*T21;u222=T10*Q^2*T21;

u311=Q^2*T30;    u312=Q*T10*Q*T31;    u313=Q*T20*Q*T32;
u321=Q*T10*Q*T31;u322=T10*Q^2*T31;    u323=T10*Q*T21*Q*T32;
u331=Q*T20*Q*T32;u332=T10*Q*T21*Q*T32;u333=T20*Q^2*T32;

u411=Q^2*T40;     u412=Q*T10*Q*T41;     u413=Q*T20*Q*T42;    u414=Q*T30*Q*T43;
u421=Q*T10*Q*T41; u422=T10*Q^2*T41;     u423=T10*Q*T21*Q*T42;u424=T10*Q*T31*Q*T43;
u431=Q*T20*Q*T42; u432=T10*Q*T21*Q*T42; u433=T20*Q^2*T42;    u434=T20*Q*T32*Q*T43;
u441=Q*T30*Q*T43; u442=T10*Q*T31*Q*T43; u443=T20*Q*T32*Q*T43;u444=T30*Q^2*T43;

u511=Q^2*T50;     u512=Q*T10*Q*T51;     u513=Q*T20*Q*T52;    u514=Q*T30*Q*T53;    u515=Q*T40*Q*T54;
u521=Q*T10*Q*T51; u522=T10*Q^2*T51;     u523=T10*Q*T21*Q*T52;u524=T10*Q*T31*Q*T53;u525=T10*Q*T41*Q*T54;
u531=Q*T20*Q*T52; u532=T10*Q*T21*Q*T52; u533=T20*Q^2*T52;    u534=T20*Q*T32*Q*T53;u535=T20*Q*T42*Q*T54;
u541=Q*T30*Q*T53; u542=T10*Q*T31*Q*T53; u543=T20*Q*T32*Q*T53;u544=T30*Q^2*T53;    u545=T30*Q*T43*Q*T54;
u551=Q*T40*Q*T54; u552=T10*Q*T41*Q*T54; u553=T20*Q*T42*Q*T54;u554=T30*Q*T43*Q*T54;u555=T40*Q^2*T54;

u611=Q^2*T60;     u612=Q*T10*Q*T61;     u613=Q*T20*Q*T62;    u614=Q*T30*Q*T63;    u615=Q*T40*Q*T64;    u616=Q*T50*Q*T65;
u621=Q*T10*Q*T61; u622=T10*Q^2*T61;     u623=T10*Q*T21*Q*T62;u624=T10*Q*T31*Q*T63;u625=T10*Q*T41*Q*T64;u626=T10*Q*T51*Q*T65;
u631=Q*T20*Q*T62; u632=T10*Q*T21*Q*T62; u633=T20*Q^2*T62;    u634=T20*Q*T32*Q*T63;u635=T20*Q*T42*Q*T64;u636=T20*Q*T52*Q*T65;
u641=Q*T30*Q*T63; u642=T10*Q*T31*Q*T63; u643=T20*Q*T32*Q*T63;u644=T30*Q^2*T63;    u645=T30*Q*T43*Q*T64;u646=T30*Q*T53*Q*T65;
u651=Q*T40*Q*T64; u652=T10*Q*T41*Q*T64; u653=T20*Q*T42*Q*T64;u654=T30*Q*T43*Q*T64;u655=T40*Q^2*T64;    u656=T40*Q*T54*Q*T65;
u661=Q*T50*Q*T65; u662=T10*Q*T51*Q*T65; u663=T20*Q*T52*Q*T65;u664=T30*Q*T53*Q*T65;u665=T40*Q*T54*Q*T65;u666=T50*Q^2*T65;

u711=Q^2*T70;     u712=Q*T10*Q*T71;     u713=Q*T20*Q*T72;    u714=Q*T30*Q*T73;    u715=Q*T40*Q*T74;    u716=Q*T50*Q*T75;    u717=Q*T60*Q*T76;
u721=Q*T10*Q*T71; u722=T10*Q^2*T71;     u723=T10*Q*T21*Q*T72;u724=T10*Q*T31*Q*T73;u725=T10*Q*T41*Q*T74;u726=T10*Q*T51*Q*T75;u727=T10*Q*T61*Q*T76;
u731=Q*T20*Q*T72; u732=T10*Q*T21*Q*T72; u733=T20*Q^2*T72;    u734=T20*Q*T32*Q*T73;u735=T20*Q*T42*Q*T74;u736=T20*Q*T52*Q*T75;u737=T20*Q*T62*Q*T76;
u741=Q*T30*Q*T73; u742=T10*Q*T31*Q*T73; u743=T20*Q*T32*Q*T73;u744=T30*Q^2*T73;    u745=T30*Q*T43*Q*T74;u746=T30*Q*T53*Q*T75;u747=T30*Q*T63*Q*T76;
u751=Q*T40*Q*T74; u752=T10*Q*T41*Q*T74; u753=T20*Q*T42*Q*T74;u754=T30*Q*T43*Q*T74;u755=T40*Q^2*T74;    u756=T40*Q*T54*Q*T75;u757=T40*Q*T64*Q*T76;
u761=Q*T50*Q*T75; u762=T10*Q*T51*Q*T75; u763=T20*Q*T52*Q*T75;u764=T30*Q*T53*Q*T75;u765=T40*Q*T54*Q*T75;u766=T50*Q^2*T75;    u767=T50*Q*T65*Q*T76;
u771=Q*T60*Q*T76; u772=T10*Q*T61*Q*T76; u773=T20*Q*T62*Q*T76;u774=T30*Q*T63*Q*T76;u775=T40*Q*T64*Q*T76;u776=T50*Q*T65*Q*T76;u777=T60*Q^2*T76;

h111=trace(u111*I1*transpose(u11))+trace(u211*I2*transpose(u21))+trace(u311*I3*transpose(u31))+trace(u411*I4*transpose(u41))+trace(u511*I5*transpose(u51))+trace(u611*I6*transpose(u61))+trace(u711*I7*transpose(u71));
h112=trace(u212*I2*transpose(u21))+trace(u312*I3*transpose(u31))+trace(u412*I4*transpose(u41))+trace(u512*I5*transpose(u51))+trace(u612*I6*transpose(u61))+trace(u712*I7*transpose(u71));
h113=trace(u313*I3*transpose(u31))+trace(u413*I4*transpose(u41))+trace(u513*I5*transpose(u51))+trace(u613*I6*transpose(u61))+trace(u713*I7*transpose(u71));
h114=trace(u414*I4*transpose(u41))+trace(u514*I5*transpose(u51))+trace(u614*I6*transpose(u61))+trace(u714*I7*transpose(u71));
h115=trace(u515*I5*transpose(u51))+trace(u615*I6*transpose(u61))+trace(u715*I7*transpose(u71));
h116=trace(u616*I6*transpose(u61))+trace(u716*I7*transpose(u71));
h117=trace(u717*I7*transpose(u71));

h121=trace(u221*I2*transpose(u21))+trace(u321*I3*transpose(u31))+trace(u421*I4*transpose(u41))+trace(u521*I5*transpose(u51))+trace(u621*I6*transpose(u61))+trace(u721*I7*transpose(u71));
h122=trace(u222*I2*transpose(u21))+trace(u322*I3*transpose(u31))+trace(u422*I4*transpose(u41))+trace(u522*I5*transpose(u51))+trace(u622*I6*transpose(u61))+trace(u722*I7*transpose(u71));
h123=trace(u323*I3*transpose(u31))+trace(u423*I4*transpose(u41))+trace(u523*I5*transpose(u51))+trace(u623*I6*transpose(u61))+trace(u723*I7*transpose(u71));
h124=trace(u424*I4*transpose(u41))+trace(u524*I5*transpose(u51))+trace(u624*I6*transpose(u61))+trace(u724*I7*transpose(u71));
h125=trace(u525*I5*transpose(u51))+trace(u625*I6*transpose(u61))+trace(u725*I7*transpose(u71));
h126=trace(u626*I6*transpose(u61))+trace(u726*I7*transpose(u71));
h127=trace(u727*I7*transpose(u71));

h131=trace(u331*I3*transpose(u31))+trace(u431*I4*transpose(u41))+trace(u531*I5*transpose(u51))+trace(u631*I6*transpose(u61))+trace(u731*I7*transpose(u71));
h132=trace(u332*I3*transpose(u31))+trace(u432*I4*transpose(u41))+trace(u532*I5*transpose(u51))+trace(u632*I6*transpose(u61))+trace(u732*I7*transpose(u71));
h133=trace(u333*I3*transpose(u31))+trace(u433*I4*transpose(u41))+trace(u533*I5*transpose(u51))+trace(u633*I6*transpose(u61))+trace(u733*I7*transpose(u71));
h134=trace(u434*I4*transpose(u41))+trace(u534*I5*transpose(u51))+trace(u634*I6*transpose(u61))+trace(u734*I7*transpose(u71));
h135=trace(u535*I5*transpose(u51))+trace(u635*I6*transpose(u61))+trace(u735*I7*transpose(u71));
h136=trace(u636*I6*transpose(u61))+trace(u736*I7*transpose(u71));
h137=trace(u737*I7*transpose(u71));

h141=trace(u441*I4*transpose(u41))+trace(u541*I5*transpose(u51))+trace(u641*I6*transpose(u61))+trace(u741*I7*transpose(u71));
h142=trace(u442*I4*transpose(u41))+trace(u542*I5*transpose(u51))+trace(u642*I6*transpose(u61))+trace(u742*I7*transpose(u71));
h143=trace(u443*I4*transpose(u41))+trace(u543*I5*transpose(u51))+trace(u643*I6*transpose(u61))+trace(u743*I7*transpose(u71));
h144=trace(u444*I4*transpose(u41))+trace(u544*I5*transpose(u51))+trace(u644*I6*transpose(u61))+trace(u744*I7*transpose(u71));
h145=trace(u545*I5*transpose(u51))+trace(u645*I6*transpose(u61))+trace(u745*I7*transpose(u71));
h146=trace(u646*I6*transpose(u61))+trace(u746*I7*transpose(u61));
h147=trace(u747*I7*transpose(u71));

h151=trace(u551*I5*transpose(u51))+trace(u651*I6*transpose(u61))+trace(u751*I7*transpose(u71));
h152=trace(u552*I5*transpose(u51))+trace(u652*I6*transpose(u61))+trace(u752*I7*transpose(u71));
h153=trace(u553*I5*transpose(u51))+trace(u653*I6*transpose(u61))+trace(u753*I7*transpose(u71));
h154=trace(u554*I5*transpose(u51))+trace(u654*I6*transpose(u61))+trace(u754*I7*transpose(u71));
h155=trace(u555*I5*transpose(u51))+trace(u655*I6*transpose(u61))+trace(u755*I7*transpose(u71));
h156=trace(u656*I6*transpose(u61))+trace(u756*I7*transpose(u71));
h157=trace(u757*I7*transpose(u71));

h161=trace(u661*I6*transpose(u61))+trace(u761*I7*transpose(u71));
h162=trace(u662*I6*transpose(u61))+trace(u762*I7*transpose(u71));
h163=trace(u663*I6*transpose(u61))+trace(u763*I7*transpose(u71));
h164=trace(u664*I6*transpose(u61))+trace(u764*I7*transpose(u71));
h165=trace(u665*I6*transpose(u61))+trace(u765*I7*transpose(u71));
h166=trace(u666*I6*transpose(u61))+trace(u766*I7*transpose(u71));
h167=trace(u767*I7*transpose(u71));

h171=trace(u771*I7*transpose(u71));
h172=trace(u772*I7*transpose(u71));
h173=trace(u773*I7*transpose(u71));
h174=trace(u774*I7*transpose(u71));
h175=trace(u775*I7*transpose(u71));
h176=trace(u776*I7*transpose(u71));
h177=trace(u777*I7*transpose(u71));

h1=h111*dq(1)*dq(1)+h112*dq(1)*dq(2)+h113*dq(1)*dq(3)+h114*dq(1)*dq(4)+h115*dq(1)*dq(5)+h116*dq(1)*dq(6)+h117*dq(1)*dq(7)+...
   h121*dq(2)*dq(1)+h122*dq(2)*dq(2)+h123*dq(2)*dq(3)+h124*dq(2)*dq(4)+h125*dq(2)*dq(5)+h126*dq(2)*dq(6)+h127*dq(2)*dq(7)+...
   h131*dq(3)*dq(1)+h132*dq(3)*dq(2)+h133*dq(3)*dq(3)+h134*dq(3)*dq(4)+h135*dq(3)*dq(5)+h136*dq(3)*dq(6)+h137*dq(3)*dq(7)+...
   h141*dq(4)*dq(1)+h142*dq(4)*dq(2)+h143*dq(4)*dq(3)+h144*dq(4)*dq(4)+h145*dq(4)*dq(5)+h146*dq(4)*dq(6)+h147*dq(4)*dq(7)+...
   h151*dq(5)*dq(1)+h152*dq(5)*dq(2)+h153*dq(5)*dq(3)+h154*dq(5)*dq(4)+h155*dq(5)*dq(5)+h156*dq(5)*dq(6)+h157*dq(5)*dq(7)+...
   h161*dq(6)*dq(1)+h162*dq(6)*dq(2)+h163*dq(6)*dq(3)+h164*dq(6)*dq(4)+h165*dq(6)*dq(5)+h166*dq(6)*dq(6)+h167*dq(6)*dq(7)+...
   h171*dq(7)*dq(1)+h172*dq(7)*dq(2)+h173*dq(7)*dq(3)+h174*dq(7)*dq(4)+h175*dq(7)*dq(5)+h176*dq(7)*dq(6)+h177*dq(7)*dq(7);

h211=trace(u211*I2*transpose(u22))+trace(u311*I3*transpose(u32))+trace(u411*I4*transpose(u42))+trace(u511*I5*transpose(u52))+trace(u611*I6*transpose(u62))+trace(u711*I7*transpose(u72));
h212=trace(u212*I2*transpose(u22))+trace(u312*I3*transpose(u32))+trace(u412*I4*transpose(u42))+trace(u512*I5*transpose(u52))+trace(u612*I6*transpose(u62))+trace(u712*I7*transpose(u72));
h213=trace(u313*I3*transpose(u32))+trace(u413*I4*transpose(u42))+trace(u513*I5*transpose(u52))+trace(u613*I6*transpose(u62))+trace(u713*I7*transpose(u72));
h214=trace(u414*I4*transpose(u42))+trace(u514*I5*transpose(u52))+trace(u614*I6*transpose(u62))+trace(u714*I7*transpose(u72));
h215=trace(u515*I5*transpose(u52))+trace(u615*I6*transpose(u62))+trace(u715*I7*transpose(u72));
h216=trace(u616*I6*transpose(u62))+trace(u716*I7*transpose(u72));
h217=trace(u717*I7*transpose(u72));

h221=trace(u221*I2*transpose(u22))+trace(u321*I3*transpose(u32))+trace(u421*I4*transpose(u42))+trace(u521*I5*transpose(u52))+trace(u621*I6*transpose(u62))+trace(u721*I7*transpose(u72));
h222=trace(u222*I2*transpose(u22))+trace(u322*I3*transpose(u32))+trace(u422*I4*transpose(u42))+trace(u522*I5*transpose(u52))+trace(u622*I6*transpose(u62))+trace(u722*I7*transpose(u72));
h223=trace(u323*I3*transpose(u32))+trace(u423*I4*transpose(u42))+trace(u523*I5*transpose(u52))+trace(u623*I6*transpose(u62))+trace(u723*I7*transpose(u72));
h224=trace(u424*I4*transpose(u42))+trace(u524*I5*transpose(u52))+trace(u624*I6*transpose(u62))+trace(u724*I7*transpose(u72));
h225=trace(u525*I5*transpose(u52))+trace(u625*I6*transpose(u62))+trace(u725*I7*transpose(u72));
h226=trace(u626*I6*transpose(u62))+trace(u726*I7*transpose(u72));
h227=trace(u727*I7*transpose(u72));

h231=trace(u331*I3*transpose(u32))+trace(u431*I4*transpose(u42))+trace(u531*I5*transpose(u52))+trace(u631*I6*transpose(u62))+trace(u731*I7*transpose(u72));
h232=trace(u332*I3*transpose(u32))+trace(u432*I4*transpose(u42))+trace(u532*I5*transpose(u52))+trace(u632*I6*transpose(u62))+trace(u732*I7*transpose(u72));
h233=trace(u333*I3*transpose(u32))+trace(u433*I4*transpose(u42))+trace(u533*I5*transpose(u52))+trace(u633*I6*transpose(u62))+trace(u733*I7*transpose(u72));
h234=trace(u434*I4*transpose(u42))+trace(u534*I5*transpose(u52))+trace(u634*I6*transpose(u62))+trace(u734*I7*transpose(u72));
h235=trace(u535*I5*transpose(u52))+trace(u635*I6*transpose(u62))+trace(u735*I7*transpose(u72));
h236=trace(u636*I6*transpose(u62))+trace(u736*I7*transpose(u72));
h237=trace(u737*I7*transpose(u72));

h241=trace(u441*I4*transpose(u42))+trace(u541*I5*transpose(u52))+trace(u641*I6*transpose(u62))+trace(u741*I7*transpose(u72));
h242=trace(u442*I4*transpose(u42))+trace(u542*I5*transpose(u52))+trace(u642*I6*transpose(u62))+trace(u742*I7*transpose(u72));
h243=trace(u443*I4*transpose(u42))+trace(u543*I5*transpose(u52))+trace(u643*I6*transpose(u62))+trace(u743*I7*transpose(u72));
h244=trace(u444*I4*transpose(u42))+trace(u544*I5*transpose(u52))+trace(u644*I6*transpose(u62))+trace(u744*I7*transpose(u72));
h245=trace(u545*I5*transpose(u52))+trace(u645*I6*transpose(u62))+trace(u745*I7*transpose(u72));
h246=trace(u646*I6*transpose(u62))+trace(u746*I7*transpose(u72));
h247=trace(u747*I7*transpose(u72));

h251=trace(u551*I5*transpose(u52))+trace(u651*I6*transpose(u62))+trace(u751*I7*transpose(u72));
h252=trace(u552*I5*transpose(u52))+trace(u652*I6*transpose(u62))+trace(u752*I7*transpose(u72));
h253=trace(u553*I5*transpose(u52))+trace(u653*I6*transpose(u62))+trace(u753*I7*transpose(u72));
h254=trace(u554*I5*transpose(u52))+trace(u654*I6*transpose(u62))+trace(u754*I7*transpose(u72));
h255=trace(u555*I5*transpose(u52))+trace(u655*I6*transpose(u62))+trace(u755*I7*transpose(u72));
h256=trace(u656*I6*transpose(u62))+trace(u756*I7*transpose(u72));
h257=trace(u757*I7*transpose(u72));

h261=trace(u661*I6*transpose(u62))+trace(u761*I7*transpose(u72));
h262=trace(u662*I6*transpose(u62))+trace(u762*I7*transpose(u72));
h263=trace(u663*I6*transpose(u62))+trace(u763*I7*transpose(u72));
h264=trace(u664*I6*transpose(u62))+trace(u764*I7*transpose(u72));
h265=trace(u665*I6*transpose(u62))+trace(u765*I7*transpose(u72));
h266=trace(u666*I6*transpose(u62))+trace(u766*I7*transpose(u72));
h267=trace(u767*I7*transpose(u72));

h271=trace(u771*I7*transpose(u72));
h272=trace(u772*I7*transpose(u72));
h273=trace(u773*I7*transpose(u72));
h274=trace(u774*I7*transpose(u72));
h275=trace(u775*I7*transpose(u72));
h276=trace(u776*I7*transpose(u72));
h277=trace(u777*I7*transpose(u72));

h2=h211*dq(1)*dq(1)+h212*dq(1)*dq(2)+h213*dq(1)*dq(3)+h214*dq(1)*dq(4)+h215*dq(1)*dq(5)+h216*dq(1)*dq(6)+h217*dq(1)*dq(7)+...
   h221*dq(2)*dq(1)+h222*dq(2)*dq(2)+h223*dq(2)*dq(3)+h224*dq(2)*dq(4)+h225*dq(2)*dq(5)+h226*dq(2)*dq(6)+h227*dq(2)*dq(7)+...
   h231*dq(3)*dq(1)+h232*dq(3)*dq(2)+h233*dq(3)*dq(3)+h234*dq(3)*dq(4)+h235*dq(3)*dq(5)+h236*dq(3)*dq(6)+h237*dq(3)*dq(7)+...
   h241*dq(4)*dq(1)+h242*dq(4)*dq(2)+h243*dq(4)*dq(3)+h244*dq(4)*dq(4)+h245*dq(4)*dq(5)+h246*dq(4)*dq(6)+h247*dq(4)*dq(7)+...
   h251*dq(5)*dq(1)+h252*dq(5)*dq(2)+h253*dq(5)*dq(3)+h254*dq(5)*dq(4)+h255*dq(5)*dq(5)+h256*dq(5)*dq(6)+h257*dq(5)*dq(7)+...
   h261*dq(6)*dq(1)+h262*dq(6)*dq(2)+h263*dq(6)*dq(3)+h264*dq(6)*dq(4)+h265*dq(6)*dq(5)+h266*dq(6)*dq(6)+h267*dq(6)*dq(7)+...
   h271*dq(7)*dq(1)+h272*dq(7)*dq(2)+h273*dq(7)*dq(3)+h274*dq(7)*dq(4)+h275*dq(7)*dq(5)+h276*dq(7)*dq(6)+h277*dq(7)*dq(7);

h311=trace(u311*I3*transpose(u33))+trace(u411*I4*transpose(u43))+trace(u511*I5*transpose(u53))+trace(u611*I6*transpose(u63))+trace(u711*I7*transpose(u73));
h312=trace(u312*I3*transpose(u33))+trace(u412*I4*transpose(u43))+trace(u512*I5*transpose(u53))+trace(u612*I6*transpose(u63))+trace(u712*I7*transpose(u73));
h313=trace(u313*I3*transpose(u33))+trace(u413*I4*transpose(u43))+trace(u513*I5*transpose(u53))+trace(u613*I6*transpose(u63))+trace(u713*I7*transpose(u73));
h314=trace(u414*I4*transpose(u43))+trace(u514*I5*transpose(u53))+trace(u614*I6*transpose(u63))+trace(u714*I7*transpose(u73));
h315=trace(u515*I5*transpose(u53))+trace(u615*I6*transpose(u63))+trace(u715*I7*transpose(u73));
h316=trace(u616*I6*transpose(u63))+trace(u716*I7*transpose(u73));
h317=trace(u717*I7*transpose(u73));

h321=trace(u321*I3*transpose(u33))+trace(u421*I4*transpose(u43))+trace(u521*I5*transpose(u53))+trace(u621*I6*transpose(u63))+trace(u721*I7*transpose(u73));
h322=trace(u322*I3*transpose(u33))+trace(u422*I4*transpose(u43))+trace(u522*I5*transpose(u53))+trace(u622*I6*transpose(u63))+trace(u722*I7*transpose(u73));
h323=trace(u323*I3*transpose(u33))+trace(u423*I4*transpose(u43))+trace(u523*I5*transpose(u53))+trace(u623*I6*transpose(u63))+trace(u723*I7*transpose(u73));
h324=trace(u424*I4*transpose(u43))+trace(u524*I5*transpose(u53))+trace(u624*I6*transpose(u63))+trace(u724*I7*transpose(u73));
h325=trace(u525*I5*transpose(u53))+trace(u625*I6*transpose(u63))+trace(u725*I7*transpose(u73));
h326=trace(u626*I6*transpose(u63))+trace(u726*I7*transpose(u73));
h327=trace(u727*I7*transpose(u73));

h331=trace(u331*I3*transpose(u33))+trace(u431*I4*transpose(u43))+trace(u531*I5*transpose(u53))+trace(u631*I6*transpose(u63))+trace(u731*I7*transpose(u73));
h332=trace(u332*I3*transpose(u33))+trace(u432*I4*transpose(u43))+trace(u532*I5*transpose(u53))+trace(u632*I6*transpose(u63))+trace(u732*I7*transpose(u73));
h333=trace(u333*I3*transpose(u33))+trace(u433*I4*transpose(u43))+trace(u533*I5*transpose(u53))+trace(u633*I6*transpose(u63))+trace(u733*I7*transpose(u73));
h334=trace(u434*I4*transpose(u43))+trace(u534*I5*transpose(u53))+trace(u634*I6*transpose(u63))+trace(u734*I7*transpose(u73));
h335=trace(u535*I5*transpose(u53))+trace(u635*I6*transpose(u63))+trace(u735*I7*transpose(u73));
h336=trace(u636*I6*transpose(u63))+trace(u736*I7*transpose(u73));
h337=trace(u737*I7*transpose(u73));

h341=trace(u441*I4*transpose(u43))+trace(u541*I5*transpose(u53))+trace(u641*I6*transpose(u63))+trace(u741*I7*transpose(u73));
h342=trace(u442*I4*transpose(u43))+trace(u542*I5*transpose(u53))+trace(u642*I6*transpose(u63))+trace(u742*I7*transpose(u73));
h343=trace(u443*I4*transpose(u43))+trace(u543*I5*transpose(u53))+trace(u643*I6*transpose(u63))+trace(u743*I7*transpose(u73));
h344=trace(u444*I4*transpose(u43))+trace(u544*I5*transpose(u53))+trace(u644*I6*transpose(u63))+trace(u744*I7*transpose(u73));
h345=trace(u545*I5*transpose(u53))+trace(u645*I6*transpose(u63))+trace(u745*I7*transpose(u73));
h346=trace(u646*I6*transpose(u63))+trace(u746*I7*transpose(u73));
h347=trace(u747*I7*transpose(u73));

h351=trace(u551*I5*transpose(u53))+trace(u651*I6*transpose(u63))+trace(u751*I7*transpose(u73));
h352=trace(u552*I5*transpose(u53))+trace(u652*I6*transpose(u63))+trace(u752*I7*transpose(u73));
h353=trace(u553*I5*transpose(u53))+trace(u653*I6*transpose(u63))+trace(u753*I7*transpose(u73));
h354=trace(u554*I5*transpose(u53))+trace(u654*I6*transpose(u63))+trace(u754*I7*transpose(u73));
h355=trace(u555*I5*transpose(u53))+trace(u655*I6*transpose(u63))+trace(u755*I7*transpose(u73));
h356=trace(u656*I6*transpose(u63))+trace(u756*I7*transpose(u73));
h357=trace(u757*I7*transpose(u73));

h361=trace(u661*I6*transpose(u63))+trace(u761*I7*transpose(u73));
h362=trace(u662*I6*transpose(u63))+trace(u762*I7*transpose(u73));
h363=trace(u663*I6*transpose(u63))+trace(u763*I7*transpose(u73));
h364=trace(u664*I6*transpose(u63))+trace(u764*I7*transpose(u73));
h365=trace(u665*I6*transpose(u63))+trace(u765*I7*transpose(u73));
h366=trace(u666*I6*transpose(u63))+trace(u766*I7*transpose(u73));
h367=trace(u767*I7*transpose(u73));

h371=trace(u771*I7*transpose(u73));
h372=trace(u772*I7*transpose(u73));
h373=trace(u773*I7*transpose(u73));
h374=trace(u774*I7*transpose(u73));
h375=trace(u775*I7*transpose(u73));
h376=trace(u776*I7*transpose(u73));
h377=trace(u777*I7*transpose(u73));

h3=h311*dq(1)*dq(1)+h312*dq(1)*dq(2)+h313*dq(1)*dq(2)+h314*dq(1)*dq(4)+h315*dq(1)*dq(5)+h316*dq(1)*dq(6)+h317*dq(1)*dq(7)+...
   h321*dq(2)*dq(1)+h322*dq(2)*dq(2)+h323*dq(2)*dq(3)+h324*dq(2)*dq(4)+h325*dq(2)*dq(5)+h326*dq(2)*dq(6)+h327*dq(2)*dq(7)+...
   h331*dq(3)*dq(1)+h332*dq(3)*dq(2)+h333*dq(3)*dq(3)+h334*dq(3)*dq(4)+h335*dq(3)*dq(5)+h336*dq(3)*dq(6)+h337*dq(3)*dq(7)+...
   h341*dq(4)*dq(1)+h342*dq(4)*dq(2)+h343*dq(4)*dq(3)+h344*dq(4)*dq(4)+h345*dq(4)*dq(5)+h346*dq(4)*dq(6)+h347*dq(4)*dq(7)+...
   h351*dq(5)*dq(1)+h352*dq(5)*dq(2)+h353*dq(5)*dq(3)+h354*dq(5)*dq(4)+h355*dq(5)*dq(5)+h356*dq(5)*dq(6)+h357*dq(5)*dq(7)+...
   h361*dq(6)*dq(1)+h362*dq(6)*dq(2)+h363*dq(6)*dq(3)+h364*dq(6)*dq(4)+h365*dq(6)*dq(5)+h366*dq(6)*dq(6)+h367*dq(6)*dq(7)+...
   h371*dq(7)*dq(1)+h372*dq(7)*dq(2)+h373*dq(7)*dq(3)+h374*dq(7)*dq(4)+h375*dq(7)*dq(5)+h376*dq(7)*dq(6)+h377*dq(7)*dq(7);

h411=trace(u411*I4*transpose(u44))+trace(u511*I5*transpose(u54))+trace(u611*I6*transpose(u64))+trace(u711*I7*transpose(u74));
h412=trace(u412*I4*transpose(u44))+trace(u512*I5*transpose(u54))+trace(u612*I6*transpose(u64))+trace(u712*I7*transpose(u74));
h413=trace(u413*I4*transpose(u44))+trace(u513*I5*transpose(u54))+trace(u613*I6*transpose(u64))+trace(u713*I7*transpose(u74));
h414=trace(u414*I4*transpose(u44))+trace(u514*I5*transpose(u54))+trace(u614*I6*transpose(u64))+trace(u714*I7*transpose(u74));
h415=trace(u515*I5*transpose(u54))+trace(u615*I6*transpose(u64))+trace(u715*I7*transpose(u74));
h416=trace(u616*I6*transpose(u64))+trace(u716*I7*transpose(u74));
h417=trace(u717*I7*transpose(u74));

h421=trace(u421*I4*transpose(u44))+trace(u521*I5*transpose(u54))+trace(u621*I6*transpose(u64))+trace(u721*I7*transpose(u74));
h422=trace(u422*I4*transpose(u44))+trace(u522*I5*transpose(u54))+trace(u622*I6*transpose(u64))+trace(u722*I7*transpose(u74));
h423=trace(u423*I4*transpose(u44))+trace(u523*I5*transpose(u54))+trace(u623*I6*transpose(u64))+trace(u723*I7*transpose(u74));
h424=trace(u424*I4*transpose(u44))+trace(u524*I5*transpose(u54))+trace(u624*I6*transpose(u64))+trace(u724*I7*transpose(u74));
h425=trace(u525*I5*transpose(u54))+trace(u625*I6*transpose(u64))+trace(u725*I7*transpose(u74));
h426=trace(u626*I6*transpose(u64))+trace(u726*I7*transpose(u74));
h427=trace(u727*I7*transpose(u74));

h431=trace(u431*I4*transpose(u44))+trace(u531*I5*transpose(u54))+trace(u631*I6*transpose(u64))+trace(u731*I7*transpose(u74));
h432=trace(u432*I4*transpose(u44))+trace(u532*I5*transpose(u54))+trace(u632*I6*transpose(u64))+trace(u732*I7*transpose(u74));
h433=trace(u433*I4*transpose(u44))+trace(u533*I5*transpose(u54))+trace(u633*I6*transpose(u64))+trace(u733*I7*transpose(u74));
h434=trace(u434*I4*transpose(u44))+trace(u534*I5*transpose(u54))+trace(u634*I6*transpose(u64))+trace(u734*I7*transpose(u74));
h435=trace(u535*I5*transpose(u54))+trace(u635*I6*transpose(u64))+trace(u735*I7*transpose(u74));
h436=trace(u636*I6*transpose(u64))+trace(u736*I7*transpose(u74));
h437=trace(u737*I7*transpose(u74));

h441=trace(u441*I4*transpose(u44))+trace(u541*I5*transpose(u54))+trace(u641*I6*transpose(u64))+trace(u741*I7*transpose(u74));
h442=trace(u442*I4*transpose(u44))+trace(u542*I5*transpose(u54))+trace(u642*I6*transpose(u64))+trace(u742*I7*transpose(u74));
h443=trace(u443*I4*transpose(u44))+trace(u543*I5*transpose(u54))+trace(u643*I6*transpose(u64))+trace(u743*I7*transpose(u74));
h444=trace(u444*I4*transpose(u44))+trace(u544*I5*transpose(u54))+trace(u644*I6*transpose(u64))+trace(u744*I7*transpose(u74));
h445=trace(u545*I5*transpose(u54))+trace(u645*I6*transpose(u64))+trace(u745*I7*transpose(u74));
h446=trace(u646*I6*transpose(u64))+trace(u746*I7*transpose(u74));
h447=trace(u747*I7*transpose(u74));

h451=trace(u551*I5*transpose(u54))+trace(u651*I6*transpose(u64))+trace(u751*I7*transpose(u74));
h452=trace(u552*I5*transpose(u54))+trace(u652*I6*transpose(u64))+trace(u752*I7*transpose(u74));
h453=trace(u553*I5*transpose(u54))+trace(u653*I6*transpose(u64))+trace(u753*I7*transpose(u74));
h454=trace(u554*I5*transpose(u54))+trace(u654*I6*transpose(u64))+trace(u754*I7*transpose(u74));
h455=trace(u555*I5*transpose(u54))+trace(u655*I6*transpose(u64))+trace(u755*I7*transpose(u74));
h456=trace(u656*I6*transpose(u64))+trace(u756*I7*transpose(u74));
h457=trace(u757*I7*transpose(u74));

h461=trace(u661*I6*transpose(u64))+trace(u761*I7*transpose(u74));
h462=trace(u662*I6*transpose(u64))+trace(u762*I7*transpose(u74));
h463=trace(u663*I6*transpose(u64))+trace(u763*I7*transpose(u74));
h464=trace(u664*I6*transpose(u64))+trace(u764*I7*transpose(u74));
h465=trace(u665*I6*transpose(u64))+trace(u765*I7*transpose(u74));
h466=trace(u666*I6*transpose(u64))+trace(u766*I7*transpose(u74));
h467=trace(u767*I7*transpose(u74));

h471=trace(u771*I7*transpose(u74));
h472=trace(u772*I7*transpose(u74));
h473=trace(u773*I7*transpose(u74));
h474=trace(u774*I7*transpose(u74));
h475=trace(u775*I7*transpose(u74));
h476=trace(u776*I7*transpose(u74));
h477=trace(u777*I7*transpose(u74));

h4=h411*dq(1)*dq(1)+h412*dq(1)*dq(2)+h413*dq(1)*dq(3)+h414*dq(1)*dq(4)+h415*dq(1)*dq(5)+h416*dq(1)*dq(6)+h417*dq(1)*dq(7)+...
   h421*dq(2)*dq(1)+h422*dq(2)*dq(2)+h423*dq(2)*dq(3)+h424*dq(2)*dq(4)+h425*dq(2)*dq(5)+h426*dq(2)*dq(6)+h427*dq(2)*dq(7)+...
   h431*dq(3)*dq(1)+h432*dq(3)*dq(2)+h433*dq(3)*dq(3)+h434*dq(3)*dq(4)+h435*dq(3)*dq(5)+h436*dq(3)*dq(6)+h437*dq(3)*dq(7)+...
   h441*dq(4)*dq(1)+h442*dq(4)*dq(2)+h443*dq(4)*dq(3)+h444*dq(4)*dq(4)+h445*dq(4)*dq(5)+h446*dq(4)*dq(6)+h447*dq(4)*dq(7)+...
   h451*dq(5)*dq(1)+h452*dq(5)*dq(2)+h453*dq(5)*dq(3)+h454*dq(5)*dq(4)+h455*dq(5)*dq(5)+h456*dq(5)*dq(6)+h457*dq(5)*dq(7)+...
   h461*dq(6)*dq(1)+h462*dq(6)*dq(2)+h463*dq(6)*dq(3)+h464*dq(6)*dq(4)+h465*dq(6)*dq(5)+h466*dq(6)*dq(6)+h467*dq(6)*dq(7)+...
   h471*dq(7)*dq(1)+h472*dq(7)*dq(2)+h473*dq(7)*dq(3)+h474*dq(7)*dq(4)+h475*dq(7)*dq(5)+h476*dq(7)*dq(6)+h477*dq(7)*dq(7);

h511=trace(u511*I5*transpose(u55))+trace(u611*I6*transpose(u65))+trace(u711*I7*transpose(u75));
h512=trace(u512*I5*transpose(u55))+trace(u612*I6*transpose(u65))+trace(u712*I7*transpose(u75));
h513=trace(u513*I5*transpose(u55))+trace(u613*I6*transpose(u65))+trace(u713*I7*transpose(u75));
h514=trace(u514*I5*transpose(u55))+trace(u614*I6*transpose(u65))+trace(u714*I7*transpose(u75));
h515=trace(u515*I5*transpose(u55))+trace(u615*I6*transpose(u65))+trace(u715*I7*transpose(u75));
h516=trace(u616*I6*transpose(u65))+trace(u716*I7*transpose(u75));
h517=trace(u717*I7*transpose(u75));

h521=trace(u521*I5*transpose(u55))+trace(u621*I6*transpose(u65))+trace(u721*I7*transpose(u75));
h522=trace(u522*I5*transpose(u55))+trace(u622*I6*transpose(u65))+trace(u722*I7*transpose(u75));
h523=trace(u523*I5*transpose(u55))+trace(u623*I6*transpose(u65))+trace(u723*I7*transpose(u75));
h524=trace(u524*I5*transpose(u55))+trace(u624*I6*transpose(u65))+trace(u724*I7*transpose(u75));
h525=trace(u525*I5*transpose(u55))+trace(u625*I6*transpose(u65))+trace(u725*I7*transpose(u75));
h526=trace(u626*I6*transpose(u65))+trace(u726*I7*transpose(u75));
h527=trace(u727*I7*transpose(u75));

h531=trace(u531*I5*transpose(u55))+trace(u631*I6*transpose(u65))+trace(u731*I7*transpose(u75));
h532=trace(u532*I5*transpose(u55))+trace(u632*I6*transpose(u65))+trace(u732*I7*transpose(u75));
h533=trace(u533*I5*transpose(u55))+trace(u633*I6*transpose(u65))+trace(u733*I7*transpose(u75));
h534=trace(u534*I5*transpose(u55))+trace(u634*I6*transpose(u65))+trace(u734*I7*transpose(u75));
h535=trace(u535*I5*transpose(u55))+trace(u635*I6*transpose(u65))+trace(u735*I7*transpose(u75));
h536=trace(u636*I6*transpose(u65))+trace(u736*I7*transpose(u75));
h537=trace(u737*I7*transpose(u75));

h541=trace(u541*I5*transpose(u55))+trace(u641*I6*transpose(u65))+trace(u741*I7*transpose(u75));
h542=trace(u542*I5*transpose(u55))+trace(u642*I6*transpose(u65))+trace(u742*I7*transpose(u75));
h543=trace(u543*I5*transpose(u55))+trace(u643*I6*transpose(u65))+trace(u743*I7*transpose(u75));
h544=trace(u544*I5*transpose(u55))+trace(u644*I6*transpose(u65))+trace(u744*I7*transpose(u75));
h545=trace(u545*I5*transpose(u55))+trace(u645*I6*transpose(u65))+trace(u745*I7*transpose(u75));
h546=trace(u646*I6*transpose(u65))+trace(u746*I7*transpose(u75));
h547=trace(u746*I7*transpose(u75));

h551=trace(u551*I5*transpose(u55))+trace(u651*I6*transpose(u65))+trace(u751*I7*transpose(u75));
h552=trace(u552*I5*transpose(u55))+trace(u652*I6*transpose(u65))+trace(u752*I7*transpose(u75));
h553=trace(u553*I5*transpose(u55))+trace(u653*I6*transpose(u65))+trace(u753*I7*transpose(u75));
h554=trace(u554*I5*transpose(u55))+trace(u654*I6*transpose(u65))+trace(u754*I7*transpose(u75));
h555=trace(u555*I5*transpose(u55))+trace(u655*I6*transpose(u65))+trace(u755*I7*transpose(u75));
h556=trace(u656*I6*transpose(u65))+trace(u756*I7*transpose(u75));
h557=trace(u757*I7*transpose(u75));

h561=trace(u661*I6*transpose(u65))+trace(u761*I7*transpose(u75));
h562=trace(u662*I6*transpose(u65))+trace(u762*I7*transpose(u75));
h563=trace(u663*I6*transpose(u65))+trace(u763*I7*transpose(u75));
h564=trace(u664*I6*transpose(u65))+trace(u764*I7*transpose(u75));
h565=trace(u665*I6*transpose(u65))+trace(u765*I7*transpose(u75));
h566=trace(u666*I6*transpose(u65))+trace(u766*I7*transpose(u75));
h567=trace(u767*I7*transpose(u75));

h571=trace(u771*I7*transpose(u75));
h572=trace(u772*I7*transpose(u75));
h573=trace(u773*I7*transpose(u75));
h574=trace(u774*I7*transpose(u75));
h575=trace(u775*I7*transpose(u75));
h576=trace(u776*I7*transpose(u75));
h577=trace(u777*I7*transpose(u75));

h5=h511*dq(1)*dq(1)+h512*dq(1)*dq(2)+h513*dq(1)*dq(3)+h514*dq(1)*dq(4)+h515*dq(1)*dq(5)+h516*dq(1)*dq(6)+h517*dq(1)*dq(7)+...
   h521*dq(2)*dq(1)+h522*dq(2)*dq(2)+h523*dq(2)*dq(3)+h524*dq(2)*dq(4)+h525*dq(2)*dq(5)+h526*dq(2)*dq(6)+h527*dq(2)*dq(7)+...
   h531*dq(3)*dq(1)+h532*dq(3)*dq(2)+h533*dq(3)*dq(3)+h534*dq(3)*dq(4)+h535*dq(3)*dq(5)+h536*dq(3)*dq(6)+h537*dq(3)*dq(7)+...
   h541*dq(4)*dq(1)+h542*dq(4)*dq(2)+h543*dq(4)*dq(3)+h544*dq(4)*dq(4)+h545*dq(4)*dq(5)+h546*dq(4)*dq(6)+h547*dq(4)*dq(7)+...
   h551*dq(5)*dq(1)+h552*dq(5)*dq(2)+h553*dq(5)*dq(3)+h554*dq(5)*dq(4)+h555*dq(5)*dq(5)+h556*dq(5)*dq(6)+h557*dq(5)*dq(7)+...
   h561*dq(6)*dq(1)+h562*dq(6)*dq(2)+h563*dq(6)*dq(3)+h564*dq(6)*dq(4)+h565*dq(6)*dq(5)+h566*dq(6)*dq(6)+h567*dq(6)*dq(7)+...
   h571*dq(7)*dq(1)+h572*dq(7)*dq(2)+h573*dq(7)*dq(3)+h574*dq(7)*dq(4)+h575*dq(7)*dq(5)+h576*dq(7)*dq(6)+h577*dq(7)*dq(7);

h611=trace(u611*I6*transpose(u66))+trace(u711*I7*transpose(u76));
h612=trace(u612*I6*transpose(u66))+trace(u712*I7*transpose(u76));
h613=trace(u613*I6*transpose(u66))+trace(u713*I7*transpose(u76));
h614=trace(u614*I6*transpose(u66))+trace(u714*I7*transpose(u76));
h615=trace(u615*I6*transpose(u66))+trace(u715*I7*transpose(u76));
h616=trace(u616*I6*transpose(u66))+trace(u716*I7*transpose(u76));
h617=trace(u717*I7*transpose(u76));

h621=trace(u621*I6*transpose(u66))+trace(u721*I7*transpose(u76));
h622=trace(u622*I6*transpose(u66))+trace(u722*I7*transpose(u76));
h623=trace(u623*I6*transpose(u66))+trace(u723*I7*transpose(u76));
h624=trace(u624*I6*transpose(u66))+trace(u724*I7*transpose(u76));
h625=trace(u625*I6*transpose(u66))+trace(u725*I7*transpose(u76));
h626=trace(u626*I6*transpose(u66))+trace(u726*I7*transpose(u76));
h627=trace(u727*I7*transpose(u76));

h631=trace(u631*I6*transpose(u66))+trace(u731*I7*transpose(u76));
h632=trace(u632*I6*transpose(u66))+trace(u732*I7*transpose(u76));
h633=trace(u633*I6*transpose(u66))+trace(u733*I7*transpose(u76));
h634=trace(u634*I6*transpose(u66))+trace(u734*I7*transpose(u76));
h635=trace(u635*I6*transpose(u66))+trace(u735*I7*transpose(u76));
h636=trace(u636*I6*transpose(u66))+trace(u736*I7*transpose(u76));
h637=trace(u737*I7*transpose(u76));

h641=trace(u641*I6*transpose(u66))+trace(u741*I7*transpose(u76));
h642=trace(u642*I6*transpose(u66))+trace(u742*I7*transpose(u76));
h643=trace(u643*I6*transpose(u66))+trace(u743*I7*transpose(u76));
h644=trace(u644*I6*transpose(u66))+trace(u744*I7*transpose(u76));
h645=trace(u645*I6*transpose(u66))+trace(u745*I7*transpose(u76));
h646=trace(u646*I6*transpose(u66))+trace(u746*I7*transpose(u76));
h647=trace(u747*I7*transpose(u76));

h651=trace(u651*I6*transpose(u66))+trace(u751*I7*transpose(u76));
h652=trace(u652*I6*transpose(u66))+trace(u752*I7*transpose(u76));
h653=trace(u653*I6*transpose(u66))+trace(u753*I7*transpose(u76));
h654=trace(u654*I6*transpose(u66))+trace(u754*I7*transpose(u76));
h655=trace(u655*I6*transpose(u66))+trace(u755*I7*transpose(u76));
h656=trace(u656*I6*transpose(u66))+trace(u756*I7*transpose(u76));
h657=trace(u757*I7*transpose(u76));

h661=trace(u661*I6*transpose(u66))+trace(u761*I7*transpose(u76));
h662=trace(u662*I6*transpose(u66))+trace(u762*I7*transpose(u76));
h663=trace(u663*I6*transpose(u66))+trace(u763*I7*transpose(u76));
h664=trace(u664*I6*transpose(u66))+trace(u764*I7*transpose(u76));
h665=trace(u665*I6*transpose(u66))+trace(u765*I7*transpose(u76));
h666=trace(u666*I6*transpose(u66))+trace(u766*I7*transpose(u76));
h667=trace(u767*I7*transpose(u76));

h671=trace(u771*I7*transpose(u76));
h672=trace(u772*I7*transpose(u76));
h673=trace(u773*I7*transpose(u76));
h674=trace(u774*I7*transpose(u76));
h675=trace(u775*I7*transpose(u76));
h676=trace(u776*I7*transpose(u76));
h677=trace(u777*I7*transpose(u76));

h6=h611*dq(1)*dq(1)+h612*dq(1)*dq(2)+h613*dq(1)*dq(3)+h614*dq(1)*dq(4)+h615*dq(1)*dq(5)+h616*dq(1)*dq(6)+h617*dq(1)*dq(7)+...
   h621*dq(2)*dq(1)+h622*dq(2)*dq(2)+h623*dq(2)*dq(3)+h624*dq(2)*dq(4)+h625*dq(2)*dq(5)+h626*dq(2)*dq(6)+h627*dq(2)*dq(7)+...
   h631*dq(3)*dq(1)+h632*dq(3)*dq(2)+h633*dq(3)*dq(3)+h634*dq(3)*dq(4)+h635*dq(3)*dq(5)+h636*dq(3)*dq(6)+h637*dq(3)*dq(7)+...
   h641*dq(4)*dq(1)+h642*dq(4)*dq(2)+h643*dq(4)*dq(3)+h644*dq(4)*dq(4)+h645*dq(4)*dq(5)+h646*dq(4)*dq(6)+h647*dq(4)*dq(7)+...
   h651*dq(5)*dq(1)+h652*dq(5)*dq(2)+h653*dq(5)*dq(3)+h654*dq(5)*dq(4)+h655*dq(5)*dq(5)+h656*dq(5)*dq(6)+h657*dq(5)*dq(7)+...
   h661*dq(6)*dq(1)+h662*dq(6)*dq(2)+h663*dq(6)*dq(3)+h664*dq(6)*dq(4)+h665*dq(6)*dq(5)+h666*dq(6)*dq(6)+h667*dq(6)*dq(7)+...
   h671*dq(7)*dq(1)+h672*dq(7)*dq(2)+h673*dq(7)*dq(3)+h674*dq(7)*dq(4)+h675*dq(7)*dq(5)+h676*dq(7)*dq(6)+h677*dq(7)*dq(7);


h711=trace(u711*I7*transpose(u77));
h712=trace(u712*I7*transpose(u77));
h713=trace(u713*I7*transpose(u77));
h714=trace(u714*I7*transpose(u77));
h715=trace(u715*I7*transpose(u77));
h716=trace(u716*I7*transpose(u77));
h717=trace(u717*I7*transpose(u77));

h721=trace(u721*I7*transpose(u77));
h722=trace(u722*I7*transpose(u77));
h723=trace(u723*I7*transpose(u77));
h724=trace(u724*I7*transpose(u77));
h725=trace(u725*I7*transpose(u77));
h726=trace(u726*I7*transpose(u77));
h727=trace(u727*I7*transpose(u77));

h731=trace(u731*I7*transpose(u77));
h732=trace(u732*I7*transpose(u77));
h733=trace(u733*I7*transpose(u77));
h734=trace(u734*I7*transpose(u77));
h735=trace(u735*I7*transpose(u77));
h736=trace(u736*I7*transpose(u77));
h737=trace(u737*I7*transpose(u77));

h741=trace(u741*I7*transpose(u77));
h742=trace(u742*I7*transpose(u77));
h743=trace(u743*I7*transpose(u77));
h744=trace(u744*I7*transpose(u77));
h745=trace(u745*I7*transpose(u77));
h746=trace(u746*I7*transpose(u77));
h747=trace(u747*I7*transpose(u77));

h751=trace(u751*I7*transpose(u77));
h752=trace(u752*I7*transpose(u77));
h753=trace(u753*I7*transpose(u77));
h754=trace(u754*I7*transpose(u77));
h755=trace(u755*I7*transpose(u77));
h756=trace(u756*I7*transpose(u77));
h757=trace(u757*I7*transpose(u77));

h761=trace(u761*I7*transpose(u77));
h762=trace(u762*I7*transpose(u77));
h763=trace(u763*I7*transpose(u77));
h764=trace(u764*I7*transpose(u77));
h765=trace(u765*I7*transpose(u77));
h766=trace(u766*I7*transpose(u77));
h767=trace(u767*I7*transpose(u77));

h771=trace(u771*I7*transpose(u77));
h772=trace(u772*I7*transpose(u77));
h773=trace(u773*I7*transpose(u77));
h774=trace(u774*I7*transpose(u77));
h775=trace(u775*I7*transpose(u77));
h776=trace(u776*I7*transpose(u77));
h777=trace(u777*I7*transpose(u77));

h7=h711*dq(1)*dq(1)+h712*dq(1)*dq(2)+h713*dq(1)*dq(3)+h714*dq(1)*dq(4)+h715*dq(1)*dq(5)+h716*dq(1)*dq(6)+h717*dq(1)*dq(7)+...
   h721*dq(2)*dq(1)+h722*dq(2)*dq(2)+h723*dq(2)*dq(3)+h724*dq(2)*dq(4)+h725*dq(2)*dq(5)+h726*dq(2)*dq(6)+h727*dq(2)*dq(7)+...
   h731*dq(3)*dq(1)+h732*dq(3)*dq(2)+h733*dq(3)*dq(3)+h734*dq(3)*dq(4)+h735*dq(3)*dq(5)+h736*dq(3)*dq(6)+h737*dq(3)*dq(7)+...
   h741*dq(4)*dq(1)+h742*dq(4)*dq(2)+h743*dq(4)*dq(3)+h744*dq(4)*dq(4)+h745*dq(4)*dq(5)+h746*dq(4)*dq(6)+h747*dq(4)*dq(7)+...
   h751*dq(5)*dq(1)+h752*dq(5)*dq(2)+h753*dq(5)*dq(3)+h754*dq(5)*dq(4)+h755*dq(5)*dq(5)+h756*dq(5)*dq(6)+h757*dq(5)*dq(7)+...
   h761*dq(6)*dq(1)+h762*dq(6)*dq(2)+h763*dq(6)*dq(3)+h764*dq(6)*dq(4)+h765*dq(6)*dq(5)+h766*dq(6)*dq(6)+h767*dq(6)*dq(7)+...
   h771*dq(7)*dq(1)+h772*dq(7)*dq(2)+h773*dq(7)*dq(3)+h774*dq(7)*dq(4)+h775*dq(7)*dq(5)+h776*dq(7)*dq(6)+h777*dq(7)*dq(7);

t1=D1+h1+g1;
t2=D2+h2+g2;
t3=D3+h3+g3;
t4=D4+h4+g4;
t5=D5+h5+g5;
t6=D6+h6+g6;
t7=D7+h7+g7;

G3=[t1;t2;t3;t4;t5;t6;t7];









