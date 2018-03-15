
%Longitud del dominio en x
Lx=1;
%Longitud del dominio en y
Ly=0.5;
%Longitud del dominio en z
Lz=0.5;


%LLenado de la matriz de diferenciacion usando diferencias finitas
%compactas tipo Pade de sexto orden
a=1.0/3.0; b=14.0/9.0; g=1.0/9.0;


%Numero de divisiones en y
ny=50;
%Numero de divisiones en z
nz=50;

%Incrementos espaciales
delta_y=Ly/(ny-1);
delta_z=Lz/(nz-1);
%Valores para la funcion de prueba f = A.cos(Kx X).cos(Ky Y)
A=0.5;%
Kx=pi/Lx;%
Ky=pi/Ly;%

nPointsError=11;
errorVec04=zeros(nPointsError,3);

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

    t1 = cputime;
    % Calculo de la matriz de diferenciacion
    matDX=inv(D1x)*D2x;
    
    texto=['Primera etapa:  ',num2str(cputime-t1),', en prueba 04. nx = ',num2str(nx)];
    disp(texto)    
    t1 = cputime;
    
    % Calculo de la matriz de diferenciacion de 2o Orden
    matD2X=matDX*matDX;
    
    texto=['Segunda etapa:  ',num2str(cputime-t1),', en prueba 04. nx = ',num2str(nx)];
    disp(texto) 

    % Computo de funcion prueba
    vecX=linspace(0,Lx,nx)';%
    vecY=linspace(0,Ly,ny)';%

    funcD0=zeros(nx,1);
    funcD2The=zeros(nx,1);
    funcD2Num=zeros(nx,1);

    funcD0=A*cos(Kx.*vecX); % Para y=0
    funcD2The=-A*(Kx^2)*cos(Kx.*vecX); % Para y=0
    funcD2Num=matD2X*funcD0;
    error=abs(funcD2Num-funcD2The);
    
    %Almacenamos valores de Error
    errorVec04(k-2,1:3)=[nx delta_x max(error)];
end

figure;ax=gca;
ax.XScale='log';ax.YScale='log';
loglog(errorVec04(:,2),errorVec04(:,3),'ko-');
grid on;