%%--------------------------------------------------------------------------------------------
% sensitivity of asymptotically optimal quantizer
%
% Author:
%   Ryo Hayakawa
% Article:
%   Ryo Hayakawa and Kazunori Hayashi,
%   "Asymptotic performance of discrete-valued vector reconstruction 
%    via box-constrained optimization with sum of l1 regularizers,"
%   IEEE Transactions on Signal Processing, vol. 68, pp. 4320-4335, 2020. 
%%--------------------------------------------------------------------------------------------


clear;
addpath('subfunctions');

% parameters
arrp0=0.1:0.02:0.9;
arrR=[-1 0 1];
Delta=0.8;
SNR=15;

L=length(arrR);
matOne=ones(L,L);

q2Iteration=100;
arrKappavsp0=zeros(size(arrp0));
for p0Index=1:length(arrp0)
  p0=arrp0(p0Index);
  disp(['p0=' num2str(p0)]);
  arrP=[(1-p0)/2 p0 (1-p0)/2];
  arrCDF=arrP*triu(matOne); % cumulative distribution

  % ternary search
  q2Min=0;
  q2Max=0.1;
  for q2IterationIndex=1:q2Iteration
    q2_1=(2*q2Min+q2Max)/3;
    q2_2=(q2Min+2*q2Max)/3;
    % q2_1
    arrCoef1=[1 q2_1 1];
    arrQ1=arrCoef1*[-1  1  1  1;...
                    -1 -1  1  1;...
                    -1 -1 -1  1;];
    arrQ1(1)=-Inf;
    arrQ1(L+1)=Inf;
    [~,~,alpha_opt1,beta_opt1]=get_theoretical(Delta,arrP,arrR,arrCoef1,arrQ1,SNR);
    % threshold of AOQ
    kappa2_1=softThr(-1/2+alpha_opt1^(2)/Delta*log((1-p0)/(2*p0)),alpha_opt1/(beta_opt1*sqrt(Delta)),arrQ1,arrR);
    arrThr1=[-Inf kappa2_1 -kappa2_1 Inf];
    SER_theo_opt1=SER_theo(alpha_opt1,beta_opt1,Delta,arrP,arrR,arrQ1,arrThr1);
    % q2_2
    arrCoef2=[1 q2_2 1];
    arrQ2=arrCoef2*[-1  1  1  1;...
                  -1 -1  1  1;...
                  -1 -1 -1  1;];
    arrQ2(1)=-Inf;
    arrQ2(L+1)=Inf;
    [~,~,alpha_opt2,beta_opt2]=get_theoretical(Delta,arrP,arrR,arrCoef2,arrQ2,SNR);
    % threshold of AOQ
    kappa2_2=softThr(-1/2+alpha_opt2^(2)/Delta*log((1-p0)/(2*p0)),alpha_opt2/(beta_opt2*sqrt(Delta)),arrQ2,arrR);
    arrThr2=[-Inf kappa2_2 -kappa2_2 Inf];
    SER_theo_opt2=SER_theo(alpha_opt2,beta_opt2,Delta,arrP,arrR,arrQ2,arrThr2);
    if SER_theo_opt1<SER_theo_opt2
      q2Max=q2_2;
    else
      q2Min=q2_1;
    end
    if abs(kappa2_2-kappa2_1)<1e-6
      if abs(q2_2-q2_1)<1e-6
        break;
      end
    end
  end
  arrKappavsp0(p0Index)=kappa2_1;
  
end

%% plot result
close all;
figure;
h=plot(arrp0,arrKappavsp0,'-k','LineWidth',1);
hold on;
grid on;
setLegend={['$\Delta=' num2str(Delta) '$']};
objLegend=legend(setLegend,'Location','northeast');
objLegend.Interpreter='latex';
objLegend.FontSize=18;
legend('boxoff');

fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
xlabel('$p_{0}$');
ylabel('$\kappa_{2}^{\ast}$');
fig.XLim=[0.1 0.9];
fig.YLim=[-1 0];
fig.XTick=0.1:0.1:0.9;
saveas(h, ['Kappavsp0_TS(Delta=' num2str(Delta) ',SNR=' num2str(SNR) ').eps'], 'epsc');
