% Script para probar la difusion numerica del esquema de diferencias
% finitas centradas de 2o orden - funcion periodica. Numero de onda aparente
close all;

%Longitud del dominio en x
Lx=2*pi;
%Max wavenumber to test ; lamda=2.pi/k
Kmax=8201;
Kmin=1;%pi/(2.0*Lx);
nKwave=83;
%Valores de numero de onda para la funcion de prueba f = exp( j Kx x)
Kwave=linspace(Kmin,Kmax,nKwave);%
KwaveNormalized=Kwave./Kmax;
nPointsK=size(Kwave,2);

%Numero de discretizaciones en potencias de 2 - nx_min=8 (2^3)
nPointsError=10;
%Matriz para almacenar las curvas de error por cada numero de onda probado
errorVec05CFDAbs=zeros(nPointsError,3,nPointsK);
errorVec05CFDIm=zeros(nPointsError,3,nPointsK);
errorVec05CFDRe=zeros(nPointsError,3,nPointsK);
%Matrices para almacenar curva de keff vs kreal por cada delta_x
vecKeffR=zeros(nPointsK,nPointsError);
vecKeffI=zeros(nPointsK,nPointsError);
vecKeffA=zeros(nPointsK,nPointsError);
vecKratioR=zeros(nPointsK,nPointsError);
vecKratioI=zeros(nPointsK,nPointsError);
vecKratioA=zeros(nPointsK,nPointsError);
%Vector de valores de nx y delta_x
vectorDeltaX=zeros(nPointsError,2);

for m = 1:nPointsError
    %N?mero de divisiones en x
    nx=2^(m+4);
    %Incrementos espaciales X
    delta_x=Lx/(nx-1);
    %Guardando valores
    vectorDeltaX(m,:)=[nx delta_x];
    %LLenado de la matriz de diferenciaci?n usando diferencias finitas
    %de segundo orden
    a=-0.5/delta_x; b=0.0; c=0.5/delta_x;
    % Para el esquema f'=D1x*f
    D1x=zeros(nx-1,nx-1);
    % LLenado de las matrices con los puntos interiores
    for l=2:nx-2
       D1x(l,l-1:l+1)= [a,b,c];
    end
    %OJO: Esta funcion es simetrica alrededor de i = 1, o i=end
    %Para el punto extremo "i=end"
    D1x(nx-1,nx-2:nx-1)=[a,b];
    D1x(nx-1,1)=c;
    %Para el punto extremo "i=1"
    D1x(1,1:2)=[b c];
    D1x(1,nx-1)=a;
    % Valores posicion X para computo de funcion prueba
    vecX=linspace(0,Lx-delta_x,nx-1)';%
    for k = 1:nPointsK
        %Seleccion del numero de onda
        Kx=Kwave(k);
        %%Crear vector para guardar los valores de la funcion
        %funcD0=zeros(nx,1);
        %%Crear vector para guardar los valores de la derivada teorica
        %funcD1The=zeros(nx,1);
        %%Crear vector para guardar los valores de la derivada numerica
        %funcD1Num=zeros(nx,1);
        %Llenando el vector de la funcion
        funcD0=exp(1j*Kx*vecX); % Para y=0
        %Llenando el vector de la derivada teorica
        funcD1The=1j*Kx*cos(Kx*vecX)-Kx*sin(Kx*vecX); % Para y=0
        %funcD1The=1j*Kx*exp(1j*Kx*vecX); % Para y=0

        %Registro de tiempo inicial
        t1 = cputime;
        %Llenando el vector de la derivada numerica. Solucion sistema directo
        funcD1Num=D1x*funcD0;
        texto=['Time der num:  ',num2str(cputime-t1),'; nx = ',num2str(nx),...
            '; Kx = ',num2str(Kx)];
        disp(texto)
        %Fin de calculo de derivada numerica 

        %Extraccion de la maxima amplitud de toda la funcion (absolute value)
        %Numero de onda effectivo absoluto
        vecKeffA(k,m)=max(abs(funcD1Num));

        %Extraccion de la parte real de la derivada numerica
        funcD1NumR=real(funcD1Num);
        %Extraccion de la maxima amplitud de la componente real
        %Numero de onda effectivo por parte real
        vecKeffR(k,m)=max(abs(funcD1NumR));
        %Extraccion de la parte imaginaria de la derivada numerica
        funcD1NumI=imag(funcD1Num);
        %Extraccion de la maxima amplitud de la componente imaginaria
        %Numero de onda effectivo por parte imaginaria
        vecKeffI(k,m)=max(abs(funcD1NumI));

        %Almacenamos valores de Numero de onda
        vecKratioR(k,m)=vecKeffR(k,m)/Kmax;
        vecKratioI(k,m)=vecKeffI(k,m)/Kmax;
        vecKratioA(k,m)=vecKeffA(k,m)/Kmax;

        %Almacenamos los valores de error (Max_norm) para esta discretizacion
        %y este valor de numero de onda
        %Abs
        errorVec05CFDA(m,1,k)=delta_x;
        errorVec05CFDA(m,2,k)=max(abs(funcD1Num-funcD1The));
        errorVec05CFDA(m,3,k)=Kx;
        %Real;
        errorVec05CFDR(m,1,k)=delta_x;
        errorVec05CFDR(m,2,k)=max(abs(real(funcD1Num)-real(funcD1The)));
        errorVec05CFDR(m,3,k)=Kx;
        %Imag;
        errorVec05CFDI(m,1,k)=delta_x;
        errorVec05CFDI(m,2,k)=max(abs(imag(funcD1Num)-imag(funcD1The)));
        errorVec05CFDI(m,3,k)=Kx;
    end
end

%--------------------
figure;ax1=gca;
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioA(:,round(nPointsError)),'ko-','DisplayName',texto1);
grid on; hold on;
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-1),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioA(:,round(nPointsError-1)),'kp:','DisplayName',texto1);
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-2),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioA(:,round(nPointsError-2)),'kd--','DisplayName',texto1);
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-3),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioA(:,round(nPointsError-3)),'ks-','DisplayName',texto1);
plot(KwaveNormalized,KwaveNormalized,'k-','DisplayName','K_{teorico}');
%ax1.YScale='log';
title('K_{eff} vs K_{real} - Usando Valores Absolutos');
legend('show')

figure;ax2=gca;
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioI(:,round(nPointsError)),'ko-','DisplayName',texto1);
grid on; hold on;
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-1),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioI(:,round(nPointsError-1)),'kp:','DisplayName',texto1);
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-2),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioI(:,round(nPointsError-2)),'kd--','DisplayName',texto1);
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-3),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioI(:,round(nPointsError-3)),'ks-','DisplayName',texto1);
plot(KwaveNormalized,KwaveNormalized,'k-','DisplayName','K_{teorico}');
%ax2.YScale='log';
title('K_{eff} vs K_{real} - Usando Valores Imaginarios');
legend('show')

figure;ax3=gca;
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioR(:,round(nPointsError)),'ko-','DisplayName',texto1);
grid on; hold on;
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-1),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioR(:,round(nPointsError-1)),'kp:','DisplayName',texto1);
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-2),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioR(:,round(nPointsError-2)),'kd--','DisplayName',texto1);
texto1=['\Delta x = ',...
    num2str(round(vectorDeltaX(round(nPointsError-3),2)*1e4)/1e4)];
plot(KwaveNormalized,vecKratioR(:,round(nPointsError-3)),'ks-','DisplayName',texto1);
plot(KwaveNormalized,KwaveNormalized,'k-','DisplayName','K_{teorico}');
%ax3.YScale='log';
title('K_{eff} vs K_{real} - Usando Valores Reales');
legend('show')

%plot(Kwave,vecKeffI,'kp:');

%figure; hold on;
%loglog(errorVec02(:,1),errorVec02(:,2),'ko-');grid on;

% figure;hold on;
% plot(vecX,funcD1The,'ko-')
% plot(vecX,funcD1Num,'rd:')
% figure;
% semilogy(vecX,error);

% subplot(1,2,1);
% spy(D1x);
% subplot(1,2,2);
% spy(D2x)
% figure;
% spy(inv(D1x));

