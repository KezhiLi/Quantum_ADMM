%% initialize
result_rho_all= zeros(5,100);
result_resid_all = zeros(5,100);
t_095_all = zeros(5,2);
fide_095_all = zeros(5,2);
%% load results
load('q8_noise40_mine.mat')

result_rho_all(1,:)= result_rho;
result_resid_all(1,:) = result_resid;
t_095_all(1,:) = t_095;
fide_095_all(1,:) = fide_095;

load('q9_noise40_mine.mat')

result_rho_all(2,:)= result_rho;
result_resid_all(2,:) = result_resid;
t_095_all(2,:) = t_095;
fide_095_all(2,:) = fide_095;

load('q10_noise40_mine.mat')

result_rho_all(3,:)= result_rho;
result_resid_all(3,:) = result_resid;
t_095_all(3,:) = t_095;
fide_095_all(3,:) = fide_095;

load('q11_noise40_mine.mat')

result_rho_all(4,:)= result_rho;
result_resid_all(4,:) = result_resid;
t_095_all(4,:) = t_095;
fide_095_all(4,:) = fide_095;

load('q12_noise40_mine.mat')

result_rho_all(5,1:70)= result_rho;
result_resid_all(5,1:70) = result_resid;
t_095_all(5,:) = t_095;
fide_095_all(5,:) = fide_095;

%% draw figures

figure, hold on,
% for ii = 1:4;
%     plot(1-result_rho_all(ii,1:50),'-*');
%     hold on,
% end
plot(1-result_rho_all(1,1:50),'-o');
plot(1-result_rho_all(2,1:50),'-x');
plot(1-result_rho_all(3,1:50),'-s');
plot(1-result_rho_all(4,1:50),'-^');
plot(1-result_rho_all(5,1:50),'-d');
grid on,
grid minor
plot([1 50],[0.95 0.95],'--k')
hold off
xlabel('No. of iterations')
ylabel('1-(Hilbert Schmidt norm)')

legend('n=8','n=9','n=10','n=11','n=12');



