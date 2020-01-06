%% DL AXMSE Replica Method 
% plot(Xsparsity,BiGAMPAMSE)
% hold on;
% plot(Xsparsity,NOBAMPAMSESmallSize);
% hold on;
% plot(Xsparsity,NOBAMPAMSEMiddleSize);
% hold on;
% plot(Xsparsity,NOBAMPAMSELargeSize);
% hold on;
% plot(TheoryrangeA,TheoryboundaryA);
% hold on;
% plot(Xsparsity,BiGAMPXMSE)
% hold on; 
% plot(Xsparsity,NOBAMPXMSESmallSize);
% hold on;
% plot(Xsparsity,NOBAMPXMSEMiddleSize);
% hold on;
% plot(Xsparsity,NOBAMPXMSELargeSize);
% hold on;
% plot(TheoryrangeX,TheoryboundaryX);
% legend('BiG-AMP','NO-BAMP Small Size','NO-BAMP Middle Size','NO-BAMP Large Size','Theory boundary');
%% MC Replica
% plot(VarName1,BiGAMPMSEA);
% hold on;
% plot(VarName1,NOBAMPMSEASmallSize);
% hold on;
% plot(VarName1,NOBAMPMSEAMiddleSize);
% hold on;
% plot(VarName1,NOBAMPMSEALargeSize);
% hold on;
% plot(VarName1,BiGAMPMSEX);
% hold on;
% plot(VarName1,NOBAMPMSEXSmallSize);
% hold on;
% plot(VarName1,NOBAMPMSEXMiddleSize);
% hold on;
% plot(VarName1,NOBAMPMSEXLargeSize);
% hold on;
% legend('BiG-AMP AMSE','NO-BAMP AMSE Small Size','NO-BAMP AMSE Middle Size','NO-BAMP AMSE Large Size','BiG-AMP XMSE','NO-BAMP XMSE Small Size','NO-BAMP XMSE Middle Size','NO-BAMP XMSE Large Size');
%% MC Compare Code
% plot(VarName1,BiGAMP);
% hold on;
% plot(VarName1,Grouse);
% hold on;
% plot(VarName1,IALM);
% hold on;
% plot(VarName1,LMaFit);
% hold on;
% plot(VarName1,NOBAMP);
% hold on;
% plot(VarName1,VBLR);
% legend('BiG-AMP','Grouse','IALM','LMaFit','NO-BAMP','VBLR');
%% RPCA Compare Code
plot(delta,BiGAMP1);
hold on;
plot(delta,IALM);
hold on;
plot(delta,SpaRCS);
hold on;
plot(delta,NOBAMP);
hold on;
plot(delta,VBLR);
hold on;
legend('BiG-AMP','IALM','SpaRCS','NO-BAMP','VBLR');
%% MC Uniform Variance Change Sparsity
% plot(FractionofY,NOBAMPNonUniformVariance);
% hold on;
% plot(FractionofY,NOBAMPUniformVariance);
% legend('NO-BAMP Non-Uniform Variance','NO-BAMP Uniform Variance');
%% MC Uniform Variance Iteration
% plot(Iteration,NOBAMPNonUniformVariance);
% hold on;
% plot(Iteration,NOBAMPUniformVariance);
% legend('NO-BAMP Non-Uniform Variance','NO-BAMP Uniform Variance');
%% MC Iteration 3 Fraction
% plot(VarName1,BigAMPfractionof02);
% hold on;
% plot(VarName1,BigAMPfractionof025);
% hold on;
% plot(VarName1,BigAMPfractionof03);
% hold on;
% plot(VarName1,BiVAMPfractionof02);
% hold on;
% plot(VarName1,BiVAMPfractionof025);
% hold on;
% plot(VarName1,BiVAMPfractionof03);
% hold on;
% legend('BiG-AMP Observed Fraction 0.2','BiG-AMP Observed Fraction 0.25','BiG-AMP Observed Fraction 0.3',...
%     'NO-BAMP Observed Fraction 0.2','NO-BAMP Observed Fraction 0.25','NO-BAMP Observed Fraction 0.3');
%% RPCA Replica Method
% plot(derta,BigAMPerrorA);
% hold on;
% plot(derta,BigAMPerrorX);
% hold on;
% plot(derta,KLerrorA);
% hold on;
% plot(derta,KLerrorX);
% hold on;
% plot(VarName6,errorAboundary);
% hold on;
% plot(VarName6,errorXboundary);
% hold on;
% legend('BigAMPerrorA','BigAMPerrorX','KLerrorA','KLerrorX','errorAboundary','errorXboundary');
%% DL SNR/Sparsity
% plot(SNR,BigAMPsparsity005);
% hold on;
% plot(SNR,BigAMPsparsity01);
% hold on;
% plot(SNR,BiGAMPsparsity03);
% hold on;
% plot(SNR,BiVAMPsparsity005);
% hold on;
% plot(SNR,BiVAMPsparsity01);
% hold on;
% plot(SNR,BiVAMPsparsity03SmallSize);
% hold on;
% plot(SNR,BiVAMPsparsity03MiddleSize);
% hold on;
% plot(SNR,BiVAMPsparsity03LargeSize);
% hold on;
% legend('BiG-AMP Sparsity 0.05','BiG-AMP Sparsity 0.1','BiG-AMP Sparsity 0.3','NO-BAMP Sparsity 0.05','NO-BAMP Sparsity 0.1','NO-BAMP Sparsity 0.3 Small Size','NO-BAMP Sparsity 0.3 Middle Size','NO-BAMP Sparsity 0.3 Large Size');
%% DL XMSE SNR Sparsity
% plot(SNR,BigAMPsparsity005);
% hold on;
% plot(SNR,BigAMPsparsity01);
% hold on;
% plot(SNR,BigAMPsparsity03);
% hold on;
% plot(SNR,BiVAMPsparsity005);
% hold on;
% plot(SNR,BiVAMPsparsity01);
% hold on;
% plot(SNR,BiVAMPsparsity03SmallSize);
% hold on;
% plot(SNR,BiVAMPsparsity03MiddleSize);
% hold on;
% plot(SNR,BiVAMPsparsity03LargeSize);
% hold on;
% legend('BiG-AMP Sparsity 0.05','BiG-AMP Sparsity 0.1','BiG-AMP Sparsity 0.3','NO-BAMP Sparsity 0.05','NO-BAMP Sparsity 0.1','NO-BAMP Sparsity 0.3 Small Size','NO-BAMP Sparsity 0.3 Middle Size','NO-BAMP Sparsity 0.3 Large Size');
%% SPCA MSE SNR
% plot(SNR,BigAMPMSEA);
% hold on;
% plot(SNR,BigAMPMSEX);
% hold on;
% plot(SNR,NBAMPMSEA);
% hold on;
% plot(SNR,NBAMPMSEX);
% legend('BiG-AMP MSE(A)','BiG-AMP MSE(X)','NO-BAMP MSE(A)','NO-BAMP MSE(X)');