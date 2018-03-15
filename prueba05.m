% Script para probar la difusi?n numerica del esquema de diferencias
% finitas


%LLenado de la matriz de diferenciaci?n usando diferencias finitas
%compactas tipo Pad? de sexto orden
a=1.0/3.0; b=14.0/9.0; g=1.0/9.0;


%Longitud del dominio en x
Lx=2*pi;
%N?mero de divisiones en x
nx=5001;
%Incrementos espaciales X
delta_x=Lx/(nx-1);

%Max wavenumber to test ; lamda=2.pi/k
Kmax=1250; 
Kmin=4*pi/Lx;
nKwave=51;
%Valores para la funcion de prueba f = exp( j Kx x)
Kwave=linspace(Kmin,Kmax,nKwave);%
   
nPointsK=size(Kwave,2);

vecKeffR=zeros(nPointsK,1);
vecKeffI=zeros(nPointsK,1);
vecKratioR=zeros(nPointsK,1);
vecKratioI=zeros(nPointsK,1);

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

% Computo de funcion prueba
vecX=linspace(0,Lx,nx)';%

for k = 1:nPointsK
    
    Kx=Kwave(k);
    
    funcD0=zeros(nx,1);
    funcD1The=zeros(nx,1);
    funcD1Num=zeros(nx,1);
    
    funcD0=exp(j*Kx*vecX); % Para y=0
    
    funcD1The=j*Kx*cos(Kx*vecX)-Kx*sin(Kx*vecX); % Para y=0
 
    t1 = cputime;
    % Solucion sistema directo
    vecRHS=D2x*funcD0;
    funcD1Num=(D1x)\vecRHS;
    texto=['Primera etapa con sol-directa:  ',num2str(cputime-t1),', en prueba 03. nx = ',num2str(nx)];
    disp(texto)
    
    %Se saca la parte real de la derivada aproximada con diferencias finitas
    funcD1NumR=real(funcD1Num);
    %Se obtiene la amplitud de la componente real
    vecKeffR(k)=max(abs(funcD1NumR));
    %Se saca la parte imaginaria de la derivada aproximada con diferencias finitas
    funcD1NumI=imag(funcD1Num);
    %Se obtiene la amplitud de la componente imaginaria
    vecKeffI(k)=max(abs(funcD1NumI));
    
    %Calculamos y Almacenamos valores de Error
    vecKratioR(k)=vecKeffR(k)/Kx;
    vecKratioI(k)=vecKeffI(k)/Kx;
    
end


figure;ax=gca;
plot(Kwave,vecKeffR,'ko-');
grid on; hold on;
plot(Kwave,vecKeffI,'kp:');
legend('KeffR','KeffI')

%figure; hold on;
%loglog(errorVec02(:,1),errorVec02(:,2),'ko-');grid on;


% figure;hold on;
% plot(vecX,funcD1The,'ko-')
% plot(vecX,funcD1Num,'rd:')
% figure;
% semilogy(vecX,error);

% for i = 1:nx
%     posX=vecX(i);
%     funcD0
%     %for j = 1:Ny
%         %posY=vecY(j);
%         %phi(j,i) = A*cos(Kx*posX)*cos(Ky*posY);
%     %end
% end



% subplot(1,2,1);
% spy(D1x);
% subplot(1,2,2);
% spy(D2x)
% figure;
% spy(inv(D1x));

