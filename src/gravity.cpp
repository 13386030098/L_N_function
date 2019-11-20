#include "ros/ros.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include </usr/include/eigen3/Eigen/Eigen>
#include <sensor_msgs/JointState.h>
#include <baxter_core_msgs/JointCommand.h>


#define pi M_PI

using namespace std;
using namespace Eigen;

double q[7]={0};
double dq[7]={0};
double ddq[7]={0};
VectorXd torque(7);
//torque << 0,0,0,0,0,0,0;

VectorXd dynamic(double *q,double *dq,double *ddq)
{
    // standard DH Parameters of the robot
    double d[7]={0.2703,0,0.3644,0,0.3743,0,0.2295};//280,0.2295
    double a[7]={0.069,0,0.069,0,0.01,0,0};
    double th[7]={q[0],q[1]+pi/2,q[2],q[3],q[4],q[5],q[6]};
    double al[7]={-pi/2,pi/2,-pi/2,pi/2,-pi/2,pi/2,0};//-pi/2,pi/2,-pi/2,pi/2,-pi/2,pi/2,0

    Matrix4d A1,A2,A3,A4,A5,A6,A7;

    A1 << cos(th[0]), -sin(th[0])*cos(al[0]), sin(th[0])*sin(al[0]),a[0]*cos(th[0]),
          sin(th[0]),  cos(th[0])*cos(al[0]),-cos(th[0])*sin(al[0]),a[0]*sin(th[0]),
                   0,             sin(al[0]),            cos(al[0]),           d[0],
                   0,                      0,                     0,             1;

   //cout <<A1<< endl<<endl;

    A2 << cos(th[1]), -sin(th[1])*cos(al[1]), sin(th[1])*sin(al[1]),a[1]*cos(th[1]),
          sin(th[1]),  cos(th[1])*cos(al[1]),-cos(th[1])*sin(al[1]),a[1]*sin(th[1]),
                   0,             sin(al[1]),            cos(al[1]),           d[1],
                   0,                      0,                     0,             1;
   // cout <<A2<< endl<<endl;

    A3 << cos(th[2]), -sin(th[2])*cos(al[2]), sin(th[2])*sin(al[2]),a[2]*cos(th[2]),
          sin(th[2]),  cos(th[2])*cos(al[2]),-cos(th[2])*sin(al[2]),a[2]*sin(th[2]),
                   0,             sin(al[2]),            cos(al[2]),           d[2],
                   0,                      0,                     0,             1;

    A4 << cos(th[3]), -sin(th[3])*cos(al[3]), sin(th[3])*sin(al[3]),a[3]*cos(th[3]),
          sin(th[3]),  cos(th[3])*cos(al[3]),-cos(th[3])*sin(al[3]),a[3]*sin(th[3]),
                   0,             sin(al[3]),            cos(al[3]),           d[3],
                   0,                      0,                     0,             1;

    A5 << cos(th[4]), -sin(th[4])*cos(al[4]), sin(th[4])*sin(al[4]),a[4]*cos(th[4]),
          sin(th[4]),  cos(th[4])*cos(al[4]),-cos(th[4])*sin(al[4]),a[4]*sin(th[4]),
                   0,             sin(al[4]),            cos(al[4]),           d[4],
                   0,                      0,                     0,             1;

    A6 << cos(th[5]), -sin(th[5])*cos(al[5]), sin(th[5])*sin(al[5]),a[5]*cos(th[5]),
          sin(th[5]),  cos(th[5])*cos(al[5]),-cos(th[5])*sin(al[5]),a[5]*sin(th[5]),
                   0,             sin(al[5]),            cos(al[5]),           d[5],
                   0,                      0,                     0,             1;

    A7 << cos(th[6]), -sin(th[6])*cos(al[6]), sin(th[6])*sin(al[6]),a[6]*cos(th[6]),
          sin(th[6]),  cos(th[6])*cos(al[6]),-cos(th[6])*sin(al[6]),a[6]*sin(th[6]),
                   0,             sin(al[6]),            cos(al[6]),           d[6],
                   0,                      0,                     0,             1;

    MatrixXd R(21,21);
    R.block(0,0,3,3)=A1.block(0,0,3,3);
    R.block(3,3,3,3)=A2.block(0,0,3,3);
    R.block(6,6,3,3)=A3.block(0,0,3,3);
    R.block(9,9,3,3)=A4.block(0,0,3,3);
    R.block(12,12,3,3)=A5.block(0,0,3,3);
    R.block(15,15,3,3)=A6.block(0,0,3,3);
    R.block(18,18,3,3)=A7.block(0,0,3,3);

    MatrixXd P(3,7);
    P.col(0)=A1.block(0,3,3,1);
    P.col(1)=A2.block(0,3,3,1);
    P.col(2)=A3.block(0,3,3,1);
    P.col(3)=A4.block(0,3,3,1);
    P.col(4)=A5.block(0,3,3,1);
    P.col(5)=A6.block(0,3,3,1);
    P.col(6)=A7.block(0,3,3,1);

    //init some variables
    Vector3d g0(0,0,-9.80665);
    Vector3d z0(0,0,1);
    VectorXd m(7);
    m << 5.70044,3.22698,4.31272,2.07206,2.24665,1.60979,0.54218;
    //cout <<m<< endl<<endl;
   //center of mass
    MatrixXd x(3,7);
    x << -0.05117, 0.00269,-0.07176, 0.00159,-0.01168,0.00697,0.005137,
          0.07908,-0.00529, 0.08149,-0.01117, 0.13111,0.00600,0.0009572,
          0.00086, 0.06845, 0.00132, 0.02618, 0.00460,0.06048,-0.06682;
    //cout <<A1<< endl<<endl;

    //惯量
    double Ixx1=0.0470910226,Iyy1=0.035959884,Izz1=0.0376697645;
    double Ixx2=0.027885975, Iyy2=0.020787492,Izz2=0.0117520941;
    double Ixx3=0.0266173355,Iyy3=0.012480083,Izz3=0.0284435520;
    double Ixx4=0.0131822787,Iyy4=0.009268520,Izz4=0.0071158268;
    double Ixx5=0.0166774282,Iyy5=0.003746311,Izz5=0.0167545726;
    double Ixx6=0.0070053791,Iyy6=0.005527552,Izz6=0.0038760715;
    double Ixx7=0.0008162135,Iyy7=0.0008735012,Izz7=0.0005494148;
    //惯量
    double Ixy1=-0.0061487003,Ixz1= 0.0001278755, Iyz1=-0.0007808689;
    double Ixy2=-0.0001882199,Ixz2=-0.00030096397,Iyz2= 0.0020767576;
    double Ixy3=-0.0039218988,Ixz3= 0.0002927063, Iyz3=-0.001083893 ;
    double Ixy4=-0.0001966341,Ixz4= 0.0003603617, Iyz4= 0.000745949 ;
    double Ixy5=-0.0001865762,Ixz5= 0.0001840370, Iyz5= 0.0006473235;
    double Ixy6= 0.0001534806,Ixz6=-0.0004438478, Iyz6=-0.0002111503;
    double Ixy7= 0.000128440, Ixz7= 0.00018969891,Iyz7= 0.0001057726;

    Matrix3d I1,I2,I3,I4,I5,I6,I7;
    MatrixXd I(21,21);

    I1 <<  Ixx1,-Ixy1,-Ixz1,
          -Ixy1, Iyy1,-Iyz1,
          -Ixz1,-Iyz1, Izz1;
    //cout <<I1<< endl<<endl;


    I2 <<  Ixx2,-Ixy2,-Ixz2,
          -Ixy2, Iyy2,-Iyz2,
          -Ixz2,-Iyz2, Izz2;

    I3 <<  Ixx3,-Ixy3,-Ixz3,
          -Ixy3, Iyy3,-Iyz3,
          -Ixz3,-Iyz3, Izz3;

    I4 <<  Ixx4,-Ixy4,-Ixz4,
          -Ixy4, Iyy4,-Iyz4,
          -Ixz4,-Iyz4, Izz4;

    I5 <<  Ixx5,-Ixy5,-Ixz5,
          -Ixy5, Iyy5,-Iyz5,
          -Ixz5,-Iyz5, Izz5;

    I6 <<  Ixx6,-Ixy6,-Ixz6,
          -Ixy6, Iyy6,-Iyz6,
          -Ixz6,-Iyz6, Izz6;

    I7 <<  Ixx7,-Ixy7,-Ixz7,
          -Ixy7, Iyy7,-Iyz7,
          -Ixz7,-Iyz7, Izz7;

    I.block(0,0,3,3)=I1;
    I.block(3,3,3,3)=I2;
    I.block(6,6,3,3)=I3;
    I.block(9,9,3,3)=I4;
    I.block(12,12,3,3)=I5;
    I.block(15,15,3,3)=I6;
    I.block(18,18,3,3)=I7;

    Vector3d w(0,0,0);
    Vector3d dw(0,0,0);
    Vector3d ddp=-g0;

    MatrixXd ddpc(3,7);
    MatrixXd nn(3,7);
    MatrixXd f(3,8);
    MatrixXd u(3,8);
    MatrixXd T(3,3);
    MatrixXd Q(3,3);
    VectorXd tou(7);
    Vector3d t1;
    Vector3d t2;
    Vector3d t3;
    Vector3d t4;
    Vector3d t5;
    Vector3d t6;
    Vector3d t7;
    tou=VectorXd::Zero(7);


    ddpc=MatrixXd::Zero(3,7);
    nn=MatrixXd::Zero(3,7);
    f=MatrixXd::Zero(3,8);
    u=MatrixXd::Zero(3,8);

    //the forward recursion 1-7

    int i=0,j=0;

    for(i=0;i<7;i++)
    {
        dw=(R.block(3*i,3*i,3,3).transpose())*(dw+z0*ddq[i]+w.cross(z0*dq[i]));
        w=(R.block(3*i,3*i,3,3).transpose())*(w+z0*dq[i]);
        t1=R.block(3*i,3*i,3,3).transpose()*P.col(i);
        t2=x.col(i);
        t3=I.block(3*i,3*i,3,3)*w;
        ddp=(R.block(3*i,3*i,3,3).transpose())*ddp+dw.cross(t1)+w.cross(w.cross(t1));
        ddpc.col(i)=ddp+dw.cross(t2)+w.cross(w.cross(t2));
        nn.col(i)=I.block(3*i,3*i,3,3)*dw+w.cross(t3);
    }

    for(j=6;j>-1;j--)
    {
        if(j==6)
        {
             T=MatrixXd::Identity(3,3);
        }
        else
        {
            T=R.block(3*(j+1),3*(j+1),3,3);
        }

        f.col(j)=T*f.col(j+1)+m(j)*ddpc.col(j);

        t4=T*f.col(j+1);
        t5=x.col(j);
        t6=(x.col(j)+(R.block(3*j,3*j,3,3).transpose())*P.col(j));
        t7=f.col(j);
        u.col(j)=nn.col(j)+T*u.col(j+1)+t4.cross(t5)+t6.cross(t7);
        Q=R.block(3*j,3*j,3,3);
        tou(j)=(u.col(j).transpose())*(Q.transpose()*z0);
    }
return tou;
}
void chatterCallback(const sensor_msgs::JointState::ConstPtr& msg)
{
    q[0]=msg->position[5];
    q[1]=msg->position[6];
    q[2]=msg->position[3];
    q[3]=msg->position[4];
    q[4]=msg->position[7];
    q[5]=msg->position[8];
    q[6]=msg->position[9];
   // cout<<"ok"<<endl;

    torque=dynamic(q,dq,ddq);
    torque=torque;
    cout<<torque<<endl<<endl;
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "gravity");
    ros::NodeHandle node;
    ros::Subscriber position_sub=node.subscribe("/robot/joint_states",1000,chatterCallback);
    ros::Publisher torque_pub=node.advertise<baxter_core_msgs::JointCommand>("/robot/limb/left/joint_command",1);
    baxter_core_msgs::JointCommand left_cmd;
    torque<<0,0,0,0,0,0,0;
    string name[7]={"left_s0","left_s1","left_e0","left_e1","left_w0","left_w1","left_w2"};
    int mode=3;

    ros::Rate loop_rate(20);
     while(ros::ok())
    {
         ros::spinOnce();
         left_cmd.mode=mode;//TORQUE_MODE=3

         left_cmd.names.resize(7);
         left_cmd.names[0]=name[0];
         left_cmd.names[1]=name[1];
         left_cmd.names[2]=name[2];
         left_cmd.names[3]=name[3];
         left_cmd.names[4]=name[4];
         left_cmd.names[5]=name[5];
         left_cmd.names[6]=name[6];

         left_cmd.command.resize(7);
         left_cmd.command[0]=0;
         left_cmd.command[1]=0;
         left_cmd.command[2]=0;
         left_cmd.command[3]=0;
         left_cmd.command[4]=0;
         left_cmd.command[5]=0;
         left_cmd.command[6]=0;
         torque_pub.publish(left_cmd);
       //  cout<<"ok"<<endl;
        loop_rate.sleep();
    }
    return 0;
}















