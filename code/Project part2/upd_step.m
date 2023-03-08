function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively

%%Implementing Measurement model
Ct=[zeros(3) zeros(3) eye(3) zeros(3) zeros(3)]; %Order of matrix is 3x15 as velocity is only considered in this kalman filter part
R=eye(3)*0.001;
Kt= (covarEst*Ct')/((Ct*covarEst*Ct')+R);  %kalman gain calculation formula

uCurr= uEst+ (Kt*(z_t - (Ct*uEst)));
covar_curr= covarEst- (Kt*Ct*covarEst); %formula for kalman filter is used as the system is already linearised in process model


return;

end