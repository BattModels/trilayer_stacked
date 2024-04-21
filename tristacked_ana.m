clear all;
cd '/Users/mbabar/Desktop/PhD/Analysis/TTLG/trilayer_stacked/'

load('aba_dos.mat')
plot(E_list, dos,'r-','linewidth', 2.0)
hold on
load('abc_dos.mat')
plot(E_list, dos,'b-','linewidth', 2.0)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',3,'GridLineStyle','--'); % ,'yscale','log'
xlabel('Energy (eV)','interpreter','latex')
ylabel('DOS','interpreter','latex')
lh = legend('ABA','ABC','northeast','FontSize',20);
legend boxoff
%title('3.9 deg nq=60')
hold off
