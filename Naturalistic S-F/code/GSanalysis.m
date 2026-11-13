%GRAPH SIGNAL ANALYSIS 
data=zX_RS;
%% compute fMRI HF/LF portions
clear X_hat X_c X_d X_m N_c N_d X_all 
for s=1:size(data,3)
    X_hat(:,:,s)=M_all(:,:,s)'*data(:,:,s);
    X_c(:,:,s)=Vlow_all(:,:,s)*X_hat(:,:,s);
    X_d(:,:,s)=Vhigh_all(:,:,s)*X_hat(:,:,s);
    X_all(:,:,s)=M_all(:,:,s)*X_hat(:,:,s); 
end


