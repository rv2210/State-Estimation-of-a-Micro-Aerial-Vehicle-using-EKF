clear; % Clear variables
datasetNum = 4; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime] = init(datasetNum);
Z = sampledVicon(7:9,:);%all the measurements that you need for the update
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %J ust for saving state his.
prevTime = 0; %last time step in real time
%write your code here calling the pred_step.m and upd_step.m functions
for i = 1:length(sampledTime)
%%extracting position orientation linear velocity angular velocity data from sampledvicon
position= sampledVicon(1:3,i);              %Getting position data from studentdata sampledvicon
orientation= sampledVicon(4:6,i);           %Getting roll pitch yaw data from studentdata sampledvicon
linear_velocity= sampledVicon(7:9,i);       %Getting linear velocity data from studentdata sampledvicon
angular_velocity= sampledVicon(10:12,i);    %Getting angular velocity data from studentdata sampledvicon
rotation_matrix=eul2rotm(orientation'); % "  ' " is added after orientation to transpose the matrix to give as input




%% inputs parameters for prediction step and update step
acc= sampledData(i).acc;            %Acceleration data extracted from student data to be given as input to prediction step function
dt= sampledTime(i)-prevTime;        %Sample time value extracted from student data to give as input to prediction step function previous time should be subtracted to get discrete time steps for the model
angVel=sampledData(i).omg;          %Angular velocity data extracted from student data
z_t= Z(:,i);                        %z_t value taken from student data to be given as input to update step function

%% calling prediction and update steps
[covarEst,uEst]=pred_step(uPrev,covarPrev,angVel,acc,dt);       %Giving uPrev,covarPrev,angVel,acc,dt as input and getting covarEst,uEst as outputs
[uCurr,covar_curr]=upd_step(z_t,covarEst,uEst);                  %Giving z_t,covarEst,uEst as inputs and getting uCurr,covar_curr as outputs
%% Saving the values of ucurr and covarcurr in saved States
savedStates(:,i)= uCurr;        % current mean data is stored to savedstates to plot the data
uPrev= uCurr;                   % making current mean to prev mean to be used for next mean calculation
covarPrev=covar_curr;           % making current covariance to prev covariance to be used for next covariance calculation
prevTime=sampledData(i).t;      % making current sample time as previous time to be used for next iteration

end



plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);