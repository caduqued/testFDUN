% Generando encabezado de documentacion


%Longitud del dominio en x
Lx=1;
%Longitud del dominio en y
Ly=0.5;
%Longitud del dominio en z
Lz=0.5;


%LLenado de la matriz de diferenciaci?n usando diferencias finitas
%compactas tipo Pad? de sexto orden
a=1.0/3.0; b=14.0/9.0; g=1.0/9.0;



%N?mero de divisiones en y
ny=50;
%N?mero de divisiones en z
nz=50;

%Incrementos espaciales
delta_x=Lx/(nx-1);
delta_y=Ly/(ny-1);
delta_z=Lz/(nz-1);
%Valores para la funcion de prueba f = A.cos(Kx X).cos(Ky Y)
A=0.5;%
Kx=pi/Lx;%
Ky=pi/Ly;%

nPointsError=11;
errorVec02=zeros(nPointsError,3);

for k = 3:nPointsError+2
    %N?mero de divisiones en x
    nx=2^k;
    %Incrementos espaciales X
    delta_x=Lx/(nx-1);
    
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
    D1x(nx,nx-1:nx)=[5,1];  
    D2x(nx,nx-5:nx)=[1/(20*delta_x) -25/(60*delta_x) 5/(3*delta_x) -5/delta_x 25/(60*delta_x) 197/(60*delta_x)];
    %Para el punto extremo "i=end-1"
    D1x(nx-1,nx-2:nx-1)=[5,1];
    D2x(nx-1,nx-6:nx-1)=[1/(20*delta_x) -25/(60*delta_x) 5/(3*delta_x) -5/delta_x 25/(60*delta_x) 197/(60*delta_x)];
    %Para el punto extremo "i=1"
    D1x(1,1:2)=[1 5];
    D2x(1,1:6)=[-197/(60*delta_x) -25/(60*delta_x) 5/delta_x -5/(3*delta_x) 25/(60*delta_x) -1/(20*delta_x)];
    %Para el punto extremo "i=2"
    D1x(2,2:3)=[1 5];
    D2x(2,2:7)=[-197/(60*delta_x) -25/(60*delta_x) 5/delta_x -5/(3*delta_x) 25/(60*delta_x) -1/(20*delta_x)];
  
    % Calculo de la matriz de diferenciacion
    matDX=inv(D1x)*D2x;

    % Computo de funcion prueba
    vecX=linspace(0,Lx,nx)';%
    vecY=linspace(0,Ly,ny)';%

    funcD0=zeros(nx,1);
    funcD1The=zeros(nx,1);
    funcD1Num=zeros(nx,1);

    funcD0=A*cos(Kx.*vecX); % Para y=0
    funcD1The=-A*Kx*sin(Kx.*vecX); % Para y=0
    funcD1Num=matDX*funcD0;
    error=abs(funcD1Num-funcD1The);
    
    %Almacenamos valores de Error
    errorVec02(k-2,1:3)=[nx delta_x max(error)];
end

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

