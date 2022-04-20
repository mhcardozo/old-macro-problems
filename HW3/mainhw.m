%%%%% HOMEWORK 3 - E622 Spring 2017 - Marcos Cardozo %%%%%

addpath C:\dynare\4.4.3\matlab\work
cd C:\dynare\4.4.3\matlab\work

%%%% PART 3 %%%%

dynare homework3 noclearall
load('homework3_results.mat','oo_')
irf3=oo_.irfs;
save irf3
load irf3

%Compute volatilities and correlations

volat3 = sqrt(diag(oo_.var));
rvolat3 = volat3./volat3(1);
save volat3
save rvolat3

X = [y c n k i r w];
correls3 = ones(length(var_list_),1)
for ii=1:length(var_list_)
correls3(ii) = (y-ones(1000,1)*oo_.mean(1))'*(X(:,ii)-ones(1000,1)*oo_.mean(ii))/(volat3(1)*volat3(ii));% Contemporaneous correlations vs Y
end


correls3t4 = oo_.autocorr{1,4}(:,1);
save correls3t4 % Correlations t-4 and t+4 with output

volat = table(var_list_,volat3) % Table with std devs
writetable(volat,'volat3.txt');

rvolat = table(var_list_,rvolat3) % Table with relative std devs vs Y
writetable(rvolat,'rvolat3.txt');

correls = table(var_list_,correls3) % Table with correlations vs Y t+-4
writetable(correls,'corrvsY.txt');



%%%% PART 4 %%%%

dynare homework4a noclearall
load('homework4a_results.mat','oo_')
irf4a=oo_.irfs;
save irf4a
load irf4a

dynare homework4b noclearall
load('homework4b_results.mat','oo_')
irf4b=oo_.irfs;
save irf4b
load irf4b

dynare homework4c noclearall
load('homework4c_results.mat','oo_')
irf4c=oo_.irfs;
save irf4c
load irf4c

dynare homework4d noclearall
load('homework4d_results.mat','oo_')
irf4d=oo_.irfs;
save irf4d
load irf4d

%Compute volatilities and correlations

volat4 = sqrt(diag(oo_.var));
rvolat4 = volat4./volat4(1);
save volat4
save rvolat4

%%%% PART 5 %%%%

dynare homework5 noclearall
load('homework5_results.mat','oo_')
irf5=oo_.irfs;
save irf5
load irf5

%Compute volatilities and correlations

volat5 = sqrt(diag(oo_.var));
rvolat5 = volat5./volat5(1);
save volat5
save rvolat5

% Create comparison of volatilities baseline vs Rogerson

volatRog = table(var_list_,volat3,volat5)
writetable(volatRog,'volatRog.txt');

%%%% PART 6 %%%%

dynare homework6 noclearall
load('homework6_results.mat','oo_')
irf6=oo_.irfs;
save irf6
load irf6

%Compute volatilities and correlations

volat6 = sqrt(diag(oo_.var));
rvolat6 = volat6./volat6(1);
save volat6
save rvolat6


%% Plot the Impulse-Response Functions (irf)

figure;         % Comparing different sigma
    x=[1:40];
    subplot(3,3,1);
    plot(x,irf3.c_e,'k',x,irf4c.c_e,'k:',x,irf4d.c_e,'k--');
    title('C');
 
    subplot(3,3,2);
    plot(x,irf3.k_e,'k',x,irf4c.k_e,'k:',x,irf4d.k_e,'k--');
    title('K');
    legend('show');
    legend('sigma=1','sigma=3','sigma=5','location','southeast');

    subplot(3,3,3);
    plot(x,irf3.n_e,'k',x,irf4c.n_e,'k:',x,irf4d.n_e,'k--');
    title('N');

    subplot(3,3,4);
    plot(x,irf3.y_e,'k',x,irf4c.y_e,'k:',x,irf4d.y_e,'k--');
    title('Y');

    subplot(3,3,5);
    plot(x,irf3.i_e,'k',x,irf4c.i_e,'k:',x,irf4d.i_e,'k--');
    title('I');

    subplot(3,3,6);
    plot(x,irf3.w_e,'k',x,irf4c.w_e,'k:',x,irf4d.w_e,'k--');
    title('W');

    subplot(3,3,7);
    plot(x,irf3.r_e,'k',x,irf4c.r_e,'k:',x,irf4d.r_e,'k--');
    title('r');

figure;                 % Comparing different psi
    x=[1:40];
    subplot(3,3,1);
    plot(x,irf3.c_e,'k',x,irf4a.c_e,'k:',x,irf4b.c_e,'k--');
    title('C');
 
    subplot(3,3,2);
    plot(x,irf3.k_e,'k',x,irf4a.k_e,'k:',x,irf4b.k_e,'k--');
    title('K');
    legend('show');
    legend('psi=1','psi=0.5','psi=2','location','southeast');

    subplot(3,3,3);
    plot(x,irf3.n_e,'k',x,irf4a.n_e,'k:',x,irf4b.n_e,'k--');
    title('N');

    subplot(3,3,4);
    plot(x,irf3.y_e,'k',x,irf4a.y_e,'k:',x,irf4b.y_e,'k--');
    title('Y');

    subplot(3,3,5);
    plot(x,irf3.i_e,'k',x,irf4a.i_e,'k:',x,irf4b.i_e,'k--');
    title('I');

    subplot(3,3,6);
    plot(x,irf3.w_e,'k',x,irf4a.w_e,'k:',x,irf4b.w_e,'k--');
    title('W');

    subplot(3,3,7);
    plot(x,irf3.r_e,'k',x,irf4a.r_e,'k:',x,irf4b.r_e,'k--');
    title('r');
    

figure;                 % Comparing baseline vs Rogerson model
    x = 1:40;
    subplot(3,3,1);
    plot(x,irf3.c_e,'k',x,irf5.c_e,'k:');
    title('C');

    subplot(3,3,2);
    plot(x,irf3.k_e,'k',x,irf5.k_e,'k:');
    title('K');
    legend('show');
    legend('Baseline','Rogerson','location','southeast');

    subplot(3,3,3);
    plot(x,irf3.n_e,'k',x,irf5.n_e,'k:');
    title('N');

    subplot(3,3,4);
    plot(x,irf3.y_e,'k',x,irf5.y_e,'k:');
    title('Y');

    subplot(3,3,5);
    plot(x,irf3.i_e,'k',x,irf5.i_e,'k:');
    title('I');

    subplot(3,3,6);
    plot(x,irf3.w_e,'k',x,irf5.w_e,'k:');
    title('W');

    subplot(3,3,7);
    plot(x,irf3.r_e,'k',x,irf5.r_e,'k:');
    title('r');
