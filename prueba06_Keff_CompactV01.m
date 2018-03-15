% Script para probar la difusion numerica del esquema de diferencias
% finitas compactas UN. Numero de onda aparente

close all;

%Longitud del dominio en x
Lx=2*pi;
%Max wavenumber to test ; lamda=2.pi/k
Kmax=1001;
Kmin=1;%pi/(2.0*Lx);
nKwave=101;
%Valores de numero de onda para la funcion de prueba f = exp( j Kx x)
Kwave=linspace(Kmin,Kmax,nKwave);%
KwaveNormalized=Kwave./Kmax;
nPointsK=size(Kwave,2);

%Numero de discretizaciones en potencias de 2 - nx_min= 2^10
nPointsError=4;
%Matriz para almacenar las curvas de error por cada numero de onda probado
errorVec06Abs=zeros(nPointsError,3,nPointsK);
errorVec06Im=zeros(nPointsError,3,nPointsK);
errorVec06Re=zeros(nPointsError,3,nPointsK);
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
    nx=2^(m+9);
    %Incrementos espaciales X
    delta_x=Lx/(nx-1);
    %Guardando valores
    vectorDeltaX(m,:)=[nx delta_x];
    %LLenado de la matriz de diferenciacion usando diferencias finitas
    %compactas tipo Pade de sexto orden
    a=1.0/3.0; b=14.0/9.0; g=1.0/9.0;

    % Para el esquema D1*f'=D2*f
    D1x=zeros(nx,nx);
    D2x=zeros(nx,nx);

    % LLenado de las matrices con los puntos interiores
    for i=3:nx-2
       D1x(i,i-1:i+1)= [a,1,a];
       D2x(i,i-2:i+2)=[-g/(4*delta_x) -b/(2*delta_x) 0 b/(2*delta_x) g/(4*delta_x)];
    end

    % Para los puntos de los extremos se usa una diferencia finita totalmente 
    % descentrada hacia atr?s o hacia delante, respectivamente
    %Para el punto extremo "i=end"
    D1x(nx,nx-2:nx)=[10,10,1];
    D2x(nx,nx-5:nx)=[-2/(60*delta_x) 25/(60*delta_x) -200/(60*delta_x) -700/(60*delta_x) 650/(60*delta_x) 227/(60*delta_x)];
    %Para el punto extremo "i=end-1"
    D1x(nx-1,nx-3:nx-1)=[10,10,1];
    D2x(nx-1,nx-6:nx-1)=[-2/(60*delta_x) 25/(60*delta_x) -200/(60*delta_x) -700/(60*delta_x) 650/(60*delta_x) 227/(60*delta_x)];
    %Para el punto extremo "i=1"
    D1x(1,1:3)=[1 10 10];
    D2x(1,1:6)=[-227/(60*delta_x) -650/(60*delta_x) 700/(60*delta_x) 200/(60*delta_x) -25/(60*delta_x) 2/(60*delta_x)];
    %Para el punto extremo "i=2"
    D1x(2,2:4)=[1 10 10];
    D2x(2,2:7)=[-227/(60*delta_x) -650/(60*delta_x) 700/(60*delta_x) 200/(60*delta_x) -25/(60*delta_x) 2/(60*delta_x)];
    
    % Valores posicion X para computo de funcion prueba
    vecX=linspace(0,Lx,nx)';%
    
    for k = 1:nPointsK
        %Seleccion del numero de onda
        Kx=Kwave(k);
        
        % Reportando avance
        texto=['Calculando para deltaX = ',num2str(delta_x),...
            ', Kx = ',num2str(Kx)];
        disp(texto);
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
        
        %Inicio de calculo de derivada numerica
        vecRHS=D2x*funcD0;
        funcD1Num=(D1x)\vecRHS;
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

