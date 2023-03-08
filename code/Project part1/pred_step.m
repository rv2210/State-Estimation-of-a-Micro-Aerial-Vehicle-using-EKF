 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time
% Process Model
%% Extracting the variable data from uPrev for defining X and X_dot and For At calculations
x= uPrev(1);
y= uPrev(2);
z= uPrev(3);
roll= uPrev(4);
pitch= uPrev(5);
yaw= uPrev(6);
vx= uPrev(7);
vy= uPrev(8);
vz= uPrev(9);
bg1= uPrev(10);
bg2= uPrev(11);
bg3= uPrev(12);
ba1= uPrev(13);
ba2= uPrev(14);
ba3= uPrev(15);
wx= angVel(1);
wy= angVel(2);
wz= angVel(3);
ax= acc(1);
ay= acc(2);
az= acc(3);
g= [0 0 -9.81]; %Defining gravity matrix and since gravity acts only in negative z direction remaining elements are zero and this is in world frame
X= uPrev; % Defining X as uPrev 
bg=[bg1;bg2;bg3];   %Defining gyroscopic bias in x, y, z axis as bg1 bg2 bg3
ba=[ba1;ba2;ba3];   %Defining accelerometer bias in x, y, z axis as ba1 ba2 ba3
am=[ax;ay;az];      %defining accelerometer matrix
wm=[wx;wy;wz];      % defining angular velocity matrix
t=wm-bg;
u=am-ba;

%% Defining elements in X_dot matrix 
X1_dot=[vx;vy;vz];                      %First element in X_Dot matrix
X2_dot=[wx-bg1-(cos(roll)*sin(pitch)*(bg3-wz))/cos(pitch)-(sin(pitch)*sin(roll)*(bg2-wy))/cos(pitch);
        sin(roll)*(bg3-wz)-cos(roll)*(bg2 - wy); 
        -(cos(roll)*(bg3-wz))/cos(pitch)-(sin(roll)*(bg2-wy))/cos(pitch)];              %2nd element in X_dot matrix
X3_dot= [(az-ba3)*(sin(roll)*sin(yaw)+cos(roll)*cos(yaw)*sin(pitch))-(ay-ba2)*(cos(roll)*sin(yaw)-cos(yaw)*sin(pitch)*sin(roll))+cos(pitch)*cos(yaw)*(ax-ba1); 
        (ay-ba2)*(cos(roll)*cos(yaw)+sin(pitch)*sin(roll)*sin(yaw))-(az-ba3)*(cos(yaw)*sin(roll)-cos(roll)*sin(pitch)*sin(yaw))+cos(pitch)*sin(yaw)*(ax-ba1); 
        cos(pitch)*cos(roll)*(az-ba3)-sin(pitch)*(ax-ba1)+cos(pitch)*sin(roll)*(ay-ba2)-(981/100)];  %3rd element in X_dot Matrix

X_dot= [X1_dot; X2_dot; X3_dot;0;0;0;0;0;0];

%Defining G Matrix "2nd element in X_dot matrix"
G= [ cos(pitch)*cos(yaw) -sin(yaw) 0;
    cos(pitch)*sin(yaw) cos(yaw) 0 ;
    -sin(pitch) 0 1;];
%G inverse matrix calculated in background script file attached with the files

G_inverse=  [cos(yaw)/cos(pitch), sin(yaw)/cos(pitch), 0;
             -sin(yaw), cos(yaw), 0; 
             (cos(yaw)*sin(pitch))/cos(pitch), (sin(pitch)*sin(yaw))/cos(pitch), 1];

%Finding rotation matrix ZYX orientation
R=rotz(yaw)*roty(pitch)*rotx(roll);
%% Defining various elements in At Matrix these elements were calculated in background live script file that is attached with the folder
PDGinR= [0, (cos(roll)*sin(pitch))/cos(pitch), -(sin(pitch)*sin(roll))/cos(pitch);
         0, -sin(roll), -cos(roll);
         0, cos(roll)/cos(pitch), -sin(roll)/cos(pitch)];
PD_Gin_R= PDGinR*t;

PDGinP= [0, -sin(roll)/(sin(pitch)^2 - 1), cos(roll)/cos(pitch)^2;
         0, 0, 0;
         0, -(sin(pitch)*sin(roll))/(sin(pitch)^2 - 1), (cos(roll)*sin(pitch))/cos(pitch)^2];
PD_Gin_P= PDGinP*t;

PDX3DotP= [-cos(yaw)*sin(pitch), cos(pitch)*cos(yaw)*sin(roll), cos(pitch)*cos(roll)*cos(yaw);
             -sin(pitch)*sin(yaw), cos(pitch)*sin(roll)*sin(yaw), cos(pitch)*cos(roll)*sin(yaw);
             -cos(pitch), -sin(pitch)*sin(roll), -cos(roll)*sin(pitch)];
PD_X3Dot_P=PDX3DotP*u;

PDX3DotR= [0, sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch), cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll);
             0, cos(roll)*sin(pitch)*sin(yaw) - cos(yaw)*sin(roll), - cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw); 
             0, cos(pitch)*cos(roll), -cos(pitch)*sin(roll)];
PD_X3Dot_R=PDX3DotR*u;

PDX3DotY= [-cos(pitch)*sin(yaw), - cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw), cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw);
              cos(pitch)*cos(yaw), cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw), sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch);
              0, 0, 0];
PD_X3Dot_Y= PDX3DotY*u;

PD_Gin_bg1= G_inverse*R*[-1;0 ; 0];
PD_Gin_bg2= G_inverse*R*[ 0;-1; 0];
PD_Gin_bg3= G_inverse*R*[ 0;0 ;-1];

PD_X3Dot_ba1= R*[-1;0 ; 0];

PD_X3Dot_ba2= R*[ 0;-1; 0];

PD_X3Dot_ba3= R*[ 0;0 ;-1];

% At matrix the elements were defined above and were calculated in background live script

At=[ 0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0;
     0 ,0 ,0 ,PD_Gin_R(1) ,PD_Gin_P(1) ,0 ,0 ,0 ,0 ,PD_Gin_bg1(1) ,PD_Gin_bg2(1) ,PD_Gin_bg3(1) ,0 ,0 ,0;
     0 ,0 ,0 ,PD_Gin_R(2) ,PD_Gin_P(2) ,0 ,0 ,0 ,0 ,PD_Gin_bg1(2) ,PD_Gin_bg2(2) ,PD_Gin_bg3(2) ,0 ,0 ,0;
     0 ,0 ,0 ,PD_Gin_R(3) ,PD_Gin_P(3) ,0 ,0 ,0 ,0 ,PD_Gin_bg1(3) ,PD_Gin_bg2(3) ,PD_Gin_bg3(3) ,0 ,0 ,0;
     0 ,0 ,0 ,PD_X3Dot_R(1) ,PD_X3Dot_P(1) ,PD_X3Dot_Y(1) ,0 ,0 ,0 ,0 ,0 ,0 ,PD_X3Dot_ba1(1) ,PD_X3Dot_ba2(2) ,PD_X3Dot_ba3(1);
     0 ,0 ,0 ,PD_X3Dot_R(2) ,PD_X3Dot_P(2) ,PD_X3Dot_Y(2) ,0 ,0 ,0 ,0 ,0 ,0 ,PD_X3Dot_ba1(2) ,PD_X3Dot_ba2(2) ,PD_X3Dot_ba3(2);
     0 ,0 ,0 ,PD_X3Dot_R(3) ,PD_X3Dot_P(3) ,PD_X3Dot_Y(3) ,0 ,0 ,0 ,0 ,0 ,0 ,PD_X3Dot_ba1(3) ,PD_X3Dot_ba2(3) ,PD_X3Dot_ba3(3);
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0; ];
 
%defining Vt matrix calculated in background live script

Vt= [ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,PD_Gin_bg1(1) ,PD_Gin_bg2(1) ,PD_Gin_bg3(1) ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,PD_Gin_bg1(2) ,PD_Gin_bg2(2) ,PD_Gin_bg3(2) ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,PD_Gin_bg1(3) ,PD_Gin_bg2(3) ,PD_Gin_bg3(3) ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,PD_X3Dot_ba1(1) ,PD_X3Dot_ba2(1) ,PD_X3Dot_ba3(1) ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,PD_X3Dot_ba1(2) ,PD_X3Dot_ba2(2) ,PD_X3Dot_ba3(2) ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,PD_X3Dot_ba1(3) ,PD_X3Dot_ba2(3) ,PD_X3Dot_ba3(3) ,0 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ; ];






%% calculating Ft Qd from At,Q
Ft= eye(15)+(dt*At);
Q=eye(15);
Qd= Q*dt;
%% implementing prediction step
uEst=uPrev+(dt*X_dot);
covarEst= (Ft*covarPrev*(Ft'))+ (Vt*Qd*(Vt'));

return;

end

