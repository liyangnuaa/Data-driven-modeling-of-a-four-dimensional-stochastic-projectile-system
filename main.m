clear;
clc;

zm=1;
zxmin=-zm;
zxmax=zm;
zymin=-zm;
zymax=zm;
zzmin=-zm;
zzmax=zm;
zwmin=-zm;
zwmax=zm;

Nz=1e4;
zxf0=(zxmax-zxmin)*rand(1,Nz)+zxmin;
zyf0=(zymax-zymin)*rand(1,Nz)+zymin;
zzf0=(zzmax-zzmin)*rand(1,Nz)+zzmin;
zwf0=(zwmax-zwmin)*rand(1,Nz)+zwmin;

%% parameters
h=0.0001;
rou=0.68;
kapa=50;
vbar=1012.3;
m00=44.2;
C00=0.12466;
A00=30.627;
S00=0.01269;
l00=2.7;
d00=0.122;
ommegaksi=79.0064;
cydot=6.5;
mzdot=-1.3546;
mzzdot=0.46;
my0=-5.8;
my2=30;

% %% Scale
% scale=100;
% zzf0=scale*zzf0;
% zwf0=scale*zwf0;

%% Generate data
Bh=sqrt(2*kapa*h)*randn(1,Nz);
zxf=(zwf0.*cos(zyf0)+zwf0.*sin(zyf0).*tan(zyf0)-1/(2*m00)*rou*vbar*S00*cydot*sin(zxf0).*cos(zyf0))*h-1/(2*m00)*rou*S00*cydot*sin(zxf0).*cos(zyf0).*Bh;
zyf=(-zzf0-1./(2*m00*cos(zxf0))*rou*vbar*S00*cydot.*sin(zyf0))*h-1./(2*m00*cos(zxf0))*rou*S00*cydot.*sin(zyf0).*Bh;
zzf=(-C00/A00*ommegaksi*zwf0+zwf0.^2.*tan(zyf0)+rou*S00*l00/(2*A00)*(-vbar^2*mzdot*sin(zyf0).*cos(zxf0)-...
    vbar*l00*mzzdot*zzf0+vbar*ommegaksi*d00.*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zxf0)))*h+...
    rou*S00*l00/(2*A00)*(-2*vbar*mzdot*sin(zyf0).*cos(zxf0)-l00*mzzdot*zzf0+ommegaksi*d00*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zxf0)).*Bh;
zwf=(C00/A00*ommegaksi*zzf0-zzf0.*zwf0.*tan(zyf0)+rou*S00*l00/(2*A00)*(vbar^2*mzdot*sin(zxf0)-vbar*l00*mzzdot*zwf0+...
    vbar*ommegaksi*d00.*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zyf0).*cos(zxf0)))*h+...
    rou*S00*l00/(2*A00)*(2*vbar*mzdot*sin(zxf0)-l00*mzzdot*zwf0+ommegaksi*d00*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zyf0).*cos(zxf0)).*Bh;

% %% Scale
% zzf=zzf/scale;
% zwf=zwf/scale;

%% Identify drift term
Ncoef=3;
zxinitial=zxf0;
zyinitial=zyf0;
zzinitial=zzf0;
zwinitial=zwf0;
x=zxf;
y=zyf;
z=zzf;
w=zwf;
n=Nz;
A=zeros(n,1);
A(:,1)=1;
EXP=[0,0,0,0];
for i=1:Ncoef
    for j=i:(-1):0
        for k=i-j:(-1):0
            for l=i-j-k:(-1):0
                A(:,i*(i+1)*(i+2)*(i+3)/24+(i-j)*(i-j+1)*(i-j+2)/6+(i-j-k)*(i-j-k+1)/2+i-j-k-l+1)=zxinitial'.^j.*zyinitial'.^k.*zzinitial'.^l.*zwinitial'.^(i-j-k-l);
                EXP=[EXP; j,k,l,i-j-k-l];
            end
        end
    end
end
A2=A;
Bx=x'/h;
By=y'/h;
Bz=z'/h;
Bw=w'/h;
lambda=0.01;

NA=length(A(1,:));
posx=1:1:NA;
for k=1:NA
    X=(A'*A)\(A'*Bx);
    I=abs(X)<lambda;
    A(:,I)=[];
    posx(I)=[];
    if isempty(X(I))
        break;
    end
end

A=A2;
posy=1:1:NA;
for k=1:NA
    Y=(A'*A)\(A'*By);
    I=abs(Y)<lambda;
    A(:,I)=[];
    posy(I)=[];
    if isempty(Y(I))
        break;
    end
end

A=A2;
posz=1:1:NA;
for k=1:NA
    Z=(A'*A)\(A'*Bz);
    I=abs(Z)<lambda;
    A(:,I)=[];
    posz(I)=[];
    if isempty(Z(I))
        break;
    end
end

A=A2;
posw=1:1:NA;
for k=1:NA
    W=(A'*A)\(A'*Bw);
    I=abs(W)<lambda;
    A(:,I)=[];
    posw(I)=[];
    if isempty(W(I))
        break;
    end
end

%% Identify diffusion term
Bxx=x.^2'/h;
Byy=y.^2'/h;
Bzz=z.^2'/h;
Bww=w.^2'/h;
lambdaBt=0.000001;

A=A2;
posBtxx=1:1:NA;
for k=1:NA
    X11=(A'*A)\(A'*Bxx);
    I=abs(X11)<lambdaBt;
    A(:,I)=[];
    posBtxx(I)=[];
    if isempty(X11(I))
        break;
    end
end

A=A2;
posBtyy=1:1:NA;
for k=1:NA
    X22=(A'*A)\(A'*Byy);
    I=abs(X22)<lambdaBt;
    A(:,I)=[];
    posBtyy(I)=[];
    if isempty(X22(I))
        break;
    end
end

A=A2;
posBtzz=1:1:NA;
for k=1:NA
    X33=(A'*A)\(A'*Bzz);
    I=abs(X33)<lambdaBt;
    A(:,I)=[];
    posBtzz(I)=[];
    if isempty(X33(I))
        break;
    end
end

A=A2;
posBtww=1:1:NA;
for k=1:NA
    X44=(A'*A)\(A'*Bww);
    I=abs(X44)<lambdaBt;
    A(:,I)=[];
    posBtww(I)=[];
    if isempty(X44(I))
        break;
    end
end

%% plot
%% plot drift 1
Nmesh=200;
x1=linspace(zxmin,zxmax,Nmesh);
x2=linspace(zymin,zymax,Nmesh);
x3=-0.5;
x4=-0.5;
[Xmesh1,Xmesh2]=meshgrid(x1,x2);
Xmesh3=x3;
Xmesh4=x4;
Flearn1=zeros(size(Xmesh1));
for i=1:length(posx)
    I=posx(i);
    Flearn1=Flearn1+X(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue1=zwf0.*cos(zyf0)+zwf0.*sin(zyf0).*tan(zyf0)-1/(2*m00)*rou*vbar*S00*cydot*sin(zxf0).*cos(zyf0);

figure;
mesh(Xmesh1,Xmesh2,Flearn1);

figure;
mesh(Xmesh1,Xmesh2,Ftrue1);


x1=-0.5;
x2=-0.5;
x3=linspace(zzmin,zzmax,Nmesh);
x4=linspace(zwmin,zwmax,Nmesh);
[Xmesh3,Xmesh4]=meshgrid(x3,x4);
Xmesh1=x1;
Xmesh2=x2;
Flearn1=zeros(size(Xmesh3));
for i=1:length(posx)
    I=posx(i);
    Flearn1=Flearn1+X(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue1=zwf0.*cos(zyf0)+zwf0.*sin(zyf0).*tan(zyf0)-1/(2*m00)*rou*vbar*S00*cydot*sin(zxf0).*cos(zyf0);

figure;
mesh(Xmesh3,Xmesh4,Flearn1);

figure;
mesh(Xmesh3,Xmesh4,Ftrue1);


%% plot drift 2
Nmesh=200;
x1=linspace(zxmin,zxmax,Nmesh);
x2=linspace(zymin,zymax,Nmesh);
x3=-0.5;
x4=-0.5;
[Xmesh1,Xmesh2]=meshgrid(x1,x2);
Xmesh3=x3;
Xmesh4=x4;
Flearn2=zeros(size(Xmesh1));
for i=1:length(posy)
    I=posy(i);
    Flearn2=Flearn2+Y(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue2=(-zzf0-1./(2*m00*cos(zxf0))*rou*vbar*S00*cydot.*sin(zyf0));

figure;
mesh(Xmesh1,Xmesh2,Flearn2);

figure;
mesh(Xmesh1,Xmesh2,Ftrue2);


x1=-0.5;
x2=-0.5;
x3=linspace(zzmin,zzmax,Nmesh);
x4=linspace(zwmin,zwmax,Nmesh);
[Xmesh3,Xmesh4]=meshgrid(x3,x4);
Xmesh1=x1;
Xmesh2=x2;
Flearn2=zeros(size(Xmesh3));
for i=1:length(posy)
    I=posy(i);
    Flearn2=Flearn2+Y(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue2=(-zzf0-1./(2*m00*cos(zxf0))*rou*vbar*S00*cydot.*sin(zyf0));

figure;
mesh(Xmesh3,Xmesh4,Flearn2);

figure;
mesh(Xmesh3,Xmesh4,Ftrue2);


%% plot drift 3
Nmesh=200;
x1=linspace(zxmin,zxmax,Nmesh);
x2=linspace(zymin,zymax,Nmesh);
x3=-0.5;
x4=-0.5;
[Xmesh1,Xmesh2]=meshgrid(x1,x2);
Xmesh3=x3;
Xmesh4=x4;
Flearn3=zeros(size(Xmesh1));
for i=1:length(posz)
    I=posz(i);
    Flearn3=Flearn3+Z(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue3=(-C00/A00*ommegaksi*zwf0+zwf0.^2.*tan(zyf0)+rou*S00*l00/(2*A00)*(-vbar^2*mzdot*sin(zyf0).*cos(zxf0)-...
    vbar*l00*mzzdot*zzf0+vbar*ommegaksi*d00.*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zxf0)));

figure;
mesh(Xmesh1,Xmesh2,Flearn3);

figure;
mesh(Xmesh1,Xmesh2,Ftrue3);


x1=-0.5;
x2=-0.5;
x3=linspace(zzmin,zzmax,Nmesh);
x4=linspace(zwmin,zwmax,Nmesh);
[Xmesh3,Xmesh4]=meshgrid(x3,x4);
Xmesh1=x1;
Xmesh2=x2;
Flearn3=zeros(size(Xmesh3));
for i=1:length(posz)
    I=posz(i);
    Flearn3=Flearn3+Z(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue3=(-C00/A00*ommegaksi*zwf0+zwf0.^2.*tan(zyf0)+rou*S00*l00/(2*A00)*(-vbar^2*mzdot*sin(zyf0).*cos(zxf0)-...
    vbar*l00*mzzdot*zzf0+vbar*ommegaksi*d00.*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zxf0)));

figure;
mesh(Xmesh3,Xmesh4,Flearn3);

figure;
mesh(Xmesh3,Xmesh4,Ftrue3);


%% plot drift 4
Nmesh=200;
x1=linspace(zxmin,zxmax,Nmesh);
x2=linspace(zymin,zymax,Nmesh);
x3=0.5;
x4=0.5;
[Xmesh1,Xmesh2]=meshgrid(x1,x2);
Xmesh3=x3;
Xmesh4=x4;
Flearn4=zeros(size(Xmesh1));
for i=1:length(posw)
    I=posw(i);
    Flearn4=Flearn4+W(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue4=(C00/A00*ommegaksi*zzf0-zzf0.*zwf0.*tan(zyf0)+rou*S00*l00/(2*A00)*(vbar^2*mzdot*sin(zxf0)-vbar*l00*mzzdot*zwf0+...
    vbar*ommegaksi*d00.*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zyf0).*cos(zxf0)));

figure;
mesh(Xmesh1,Xmesh2,Flearn4);

figure;
mesh(Xmesh1,Xmesh2,Ftrue4);


x1=0.5;
x2=0.5;
x3=linspace(zzmin,zzmax,Nmesh);
x4=linspace(zwmin,zwmax,Nmesh);
[Xmesh3,Xmesh4]=meshgrid(x3,x4);
Xmesh1=x1;
Xmesh2=x2;
Flearn4=zeros(size(Xmesh3));
for i=1:length(posw)
    I=posw(i);
    Flearn4=Flearn4+W(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue4=(C00/A00*ommegaksi*zzf0-zzf0.*zwf0.*tan(zyf0)+rou*S00*l00/(2*A00)*(vbar^2*mzdot*sin(zxf0)-vbar*l00*mzzdot*zwf0+...
    vbar*ommegaksi*d00.*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zyf0).*cos(zxf0)));

figure;
mesh(Xmesh3,Xmesh4,Flearn4);

figure;
mesh(Xmesh3,Xmesh4,Ftrue4);


%% plot diffusion 11
Nmesh=200;
x1=linspace(zxmin,zxmax,Nmesh);
x2=linspace(zymin,zymax,Nmesh);
x3=0.5;
x4=-0.5;
[Xmesh1,Xmesh2]=meshgrid(x1,x2);
Xmesh3=x3;
Xmesh4=x4;
Flearn11=zeros(size(Xmesh1));
for i=1:length(posBtxx)
    I=posBtxx(i);
    Flearn11=Flearn11+X11(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end
Flearn11=Flearn11/(2*kapa);

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue11=1/(2*m00)*rou*S00*cydot.*sin(zxf0).*cos(zyf0)*1/(2*m00)*rou*S00*cydot.*sin(zxf0).*cos(zyf0);

figure;
mesh(Xmesh1,Xmesh2,Flearn11);

figure;
mesh(Xmesh1,Xmesh2,Ftrue11);


%% plot diffusion 22
Nmesh=200;
x1=linspace(zxmin,zxmax,Nmesh);
x2=linspace(zymin,zymax,Nmesh);
x3=0.5;
x4=-0.5;
[Xmesh1,Xmesh2]=meshgrid(x1,x2);
Xmesh3=x3;
Xmesh4=x4;
Flearn22=zeros(size(Xmesh1));
for i=1:length(posBtyy)
    I=posBtyy(i);
    Flearn22=Flearn22+X22(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end
Flearn22=Flearn22/(2*kapa);

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue22=1./(2*m00*cos(zxf0))*rou*S00*cydot.*sin(zyf0)*1./(2*m00*cos(zxf0))*rou*S00*cydot.*sin(zyf0);

figure;
mesh(Xmesh1,Xmesh2,Flearn22);

figure;
mesh(Xmesh1,Xmesh2,Ftrue22);


%% plot diffusion 33
Nmesh=200;
x1=linspace(zxmin,zxmax,Nmesh);
x2=linspace(zymin,zymax,Nmesh);
x3=0.5;
x4=-0.5;
[Xmesh1,Xmesh2]=meshgrid(x1,x2);
Xmesh3=x3;
Xmesh4=x4;
Flearn33=zeros(size(Xmesh1));
for i=1:length(posBtzz)
    I=posBtzz(i);
    Flearn33=Flearn33+X33(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end
Flearn33=Flearn33/(2*kapa);

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue33=rou*S00*l00/(2*A00)*(-2*vbar*mzdot*sin(zyf0).*cos(zxf0)-l00*mzzdot*zzf0+ommegaksi*d00*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zxf0))*...
    rou*S00*l00/(2*A00).*(-2*vbar*mzdot*sin(zyf0).*cos(zxf0)-l00*mzzdot*zzf0+ommegaksi*d00*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zxf0));

figure;
mesh(Xmesh1,Xmesh2,Flearn33);

figure;
mesh(Xmesh1,Xmesh2,Ftrue33);


%% plot diffusion 44
Nmesh=200;
x1=linspace(zxmin,zxmax,Nmesh);
x2=linspace(zymin,zymax,Nmesh);
x3=0.5;
x4=-0.5;
[Xmesh1,Xmesh2]=meshgrid(x1,x2);
Xmesh3=x3;
Xmesh4=x4;
Flearn44=zeros(size(Xmesh1));
for i=1:length(posBtww)
    I=posBtww(i);
    Flearn44=Flearn44+X44(i)*Xmesh1.^EXP(I,1).*Xmesh2.^EXP(I,2).*Xmesh3.^EXP(I,3).*Xmesh4.^EXP(I,4);
end
Flearn44=Flearn44/(2*kapa);

zxf0=Xmesh1;
zyf0=Xmesh2;
zzf0=Xmesh3;
zwf0=Xmesh4;
Ftrue44=rou*S00*l00/(2*A00)*(2*vbar*mzdot*sin(zxf0)-l00*mzzdot*zwf0+ommegaksi*d00*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zyf0).*cos(zxf0))*...
    rou*S00*l00/(2*A00).*(2*vbar*mzdot*sin(zxf0)-l00*mzzdot*zwf0+ommegaksi*d00*(my0+my2*((sin(zxf0).^2.*(cos(zyf0).^2))+(sin(zyf0)).^2)).*sin(zyf0).*cos(zxf0));

figure;
mesh(Xmesh1,Xmesh2,Flearn44);

figure;
mesh(Xmesh1,Xmesh2,Ftrue44);

