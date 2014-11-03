function sumll = novoCIRMW(para,input, datescup,parametros,cupao,mw) %calculate log likelihood
format long
warning('off')
estado1 = input(:,3);
estado2 = input(:,4);

Z = input(:,2)/100;

datesobs = input(:,1);
datesobsmatlab = datesobs + datenum('30-Dec-1899');
datescupmatlab = datescup + datenum('30-Dec-1899');
[nrow,ncol] = size(datesobsmatlab); %numero de observações e de obrigacoes

%parametros a estimar
theta3 = para(1); k3 = para(2); sigma3 = para(3); eta3 = para(4);
theta5 = para(5); k5 = para(6); sigma5 = para(7); eta5 = para(8);
gama0 = para(9); chi0 = para(10); betad1 = para(11); betad2 = para(12); chi3 = para(13); 
sigmai = para(14);

%paramteros já estimados
theta1 = parametros(1,1) ; theta2 = parametros(1,2); alfa0 = parametros(1,3); 
k1 = parametros(2,1) ; k2 = parametros(2,2);
sigma1 = parametros(3,1) ; sigma2 = parametros(3,2); 
eta1 = parametros(4,1); eta2 = parametros(4,2);


%n de cupoes
ncupaux = size(datescup);
ncup = ncupaux(1);

h1=sqrt(((k1+eta1)^2)+2*(1+betad1)*sigma1*sigma1);
h2=sqrt(((k2+eta2)^2)+2*(1+betad2)*sigma2*sigma2);
h5=sqrt(((k5+eta5)^2)+2*sigma5*sigma5);

AffineAlpha1 = @(tau) ((2*k1*theta1)/(sigma1*sigma1))*log((2*h1*exp(0.5*(k1+eta1+h1)*tau))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau)));
AffineAlpha2 = @(tau) ((2*k2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(k2+eta2+h2)*tau))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau)));
AffineAlpha5 = @(tau) ((2*k5*theta5)/(sigma5*sigma5))*log((2*h5*exp(0.5*(k5+eta5+h5)*tau))/(h5-(k5+eta5)+(k5+eta5+h5)*exp(h5*tau)));

AffineBeta1 = @(tau) (2*(1+betad1)*(exp(h1*tau)-1))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau));
AffineBeta2 = @(tau) (2*(1+betad2)*(exp(h2*tau)-1))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau));
AffineBeta5 = @(tau) (2*(exp(h5*tau)-1))/(h5-(k5+eta5)+(k5+eta5+h5)*exp(h5*tau));

f = @(t,valor1,valor2) exp(-(alfa0+gama0)*t+AffineAlpha1(t)+AffineAlpha2(t)-AffineBeta1(t)*valor1-AffineBeta2(t)*valor2);
g = @(t,valor5) exp(-chi0*t+AffineAlpha5(t)-AffineBeta5(t)*valor5);

kimmela = ((2*(theta3^2)*(k3^2))/(sigma3^4))-((2*theta3*k3)/(sigma3^2))+((4*mw*chi3)/(sigma3^2))+(3/8);
kimmelb = 0.5*sqrt(((k3+eta3)^2)+2*(sigma3^2));
kimmeld = (-1)*((theta3*k3*(k3+eta3))/(sigma3^2));
kimmelalfa = ((2*theta3*k3)/(sigma3^2))-0.5;
kimmelgama = 0.5*(1-sqrt(1+8*kimmela));
kimmeldelta = @(t) 1-exp(-2*kimmelb*t);
kimmely = @(valor3) (2*sqrt(valor3))/(sigma3);
kimmelz = @(t,valor3) sqrt(2*kimmelb)*exp(-kimmelb*t)*kimmely(valor3);

auxiliarpi1= @(t,valor3) ((kimmelalfa-kimmelgama)*(kimmelalfa+kimmelgama-1)/(2*(kimmelz(t,valor3)^2)))+((2*kimmelalfa+1)/4)*(1-((k3+eta3)/(2*kimmelb)))+(kimmelz(t,valor3)^2)*(1/8)*((1-((k3+eta3)/(2*kimmelb)))^2);
auxiliarpi2 = @(t,valor3) ((kimmelalfa-kimmelgama)*(kimmelalfa-kimmelgama-2)*(kimmelalfa+kimmelgama-1)*(kimmelalfa+kimmelgama-3))/(4*(kimmelz(t,valor3)^4));
auxiliarpi3 = @(t,valor3) ((2*kimmelalfa-1)*(kimmelalfa-kimmelgama)*(kimmelalfa+kimmelgama-1)*(1-((k3+eta3)/(2*kimmelb))))/(4*kimmelz(t,valor3)^2);
auxiliarpi4 = @(t,valor3) ((1/16)*(2*kimmelalfa+3)*(2*kimmelalfa+1)+(1/8)*(kimmelalfa-kimmelgama)*(kimmelalfa+kimmelgama-1))*((1-((k3+eta3)/(2*kimmelb)))^2);
auxiliarpi5 = @(t,valor3) (1/16)*(2*kimmelalfa+3)*(kimmelz(t,valor3)^2)*((1-((k3+eta3)/(2*kimmelb)))^3)+(1/64)*(kimmelz(t,valor3)^4)*((1-((k3+eta3)/(2*kimmelb)))^4);

kimmelpi = @(t,valor3) 1+kimmeldelta(t)*auxiliarpi1(t,valor3)+0.5*(kimmeldelta(t)^2)*(auxiliarpi2(t,valor3)+auxiliarpi3(t,valor3)+auxiliarpi4(t,valor3)+auxiliarpi5(t,valor3));
kimmelwum = @(t,valor3) ((kimmelz(t,valor3)/(sqrt(2*kimmelb)))^(kimmelalfa-kimmelgama))*exp((1/4)*(1-((k3+eta3)/(2*kimmelb)))*(kimmelz(t,valor3)^2))*kimmelpi(t,valor3);
kimmelh = @(t,valor3) ((kimmelz(t,valor3)/(sqrt(2*kimmelb)))^kimmelgama)*exp(-0.5*kimmelb*(kimmely(valor3)^2)-(0.5*kimmelb+kimmeld)*t)*kimmelwum(t,valor3);

kimmelpifinal = @(t,valor3) ((4*valor3/(sigma3^2))^(0.25-((theta3*k3)/(sigma3^2))))*exp((valor3*(k3+eta3))/(sigma3^2))*kimmelh(t,valor3);

derivadaauxiliar1= @(t,valor3) ((-(kimmelalfa-kimmelgama)*(kimmelalfa+kimmelgama-1))/((kimmelz(t,valor3)^3)))+(kimmelz(t,valor3))*(1/4)*((1-((k3+eta3)/(2*kimmelb)))^2);
derivadaauxiliar2 = @(t,valor3) ((-(kimmelalfa-kimmelgama)*(kimmelalfa-kimmelgama-2)*(kimmelalfa+kimmelgama-1)*(kimmelalfa+kimmelgama-3)))/((kimmelz(t,valor3)^5));
derivadaauxiliar3 = @(t,valor3) ((2*kimmelalfa-1)*(kimmelalfa-kimmelgama)*(kimmelalfa+kimmelgama-1)*(1-((k3+eta3)/(2*kimmelb))))/(2*kimmelz(t,valor3)^3);
derivadaauxiliar4 = @(t,valor3) (1/8)*(2*kimmelalfa+3)*(kimmelz(t,valor3))*((1-((k3+eta3)/(2*kimmelb)))^3)+(1/16)*(kimmelz(t,valor3)^3)*((1-((k3+eta3)/(2*kimmelb)))^4);

derivadakimmelinicial = @(t,valor3) kimmeldelta(t)*derivadaauxiliar1(t,valor3)+0.5*(kimmeldelta(t)^2)*(derivadaauxiliar2(t,valor3)+derivadaauxiliar3(t,valor3)+derivadaauxiliar4(t,valor3));

derivadakimmelfinalauxiliar1= @(t,valor3) (4/(sigma3^2))*(0.25-((theta3*k3)/(sigma3^2)))*(((4*valor3)/(sigma3^2))^(-0.75-((theta3*k3)/(sigma3^2))))*exp((0.5*(k3+eta3)-((0.5*kimmelb)/(sigma3^2)))*valor3)*exp((-0.5*kimmelb-kimmeld)*t)*exp(0.25*(kimmelz(t,valor3)^2)*(1-((k3+eta3)/(2*kimmelb))))*kimmelpifinal(t,valor3);  
derivadakimmelfinalauxiliar2= @(t,valor3) (((4*valor3)/(sigma3^2))^(0.25-((theta3*k3)/(sigma3^2))))*(0.5*(k3*eta3)-2*(kimmelb/(sigma3^3)))*exp((0.5*(k3+eta3)-((0.5*kimmelb)/(sigma3^2)))*valor3)*exp((-0.5*kimmelb-kimmeld)*t)*exp(0.25*(kimmelz(t,valor3)^2)*(1-((k3+eta3)/(2*kimmelb))))*kimmelpifinal(t,valor3);
derivadakimmelfinalauxiliar3= @(t,valor3) (((4*valor3)/(sigma3^2))^(0.25-((theta3*k3)/(sigma3^2))))*exp((0.5*(k3+eta3)-((0.5*kimmelb)/(sigma3^2)))*valor3)*exp((-0.5*kimmelb-kimmeld)*t)*(sqrt(2*kimmelb)/(sigma3))*exp(-kimmelb*t)*(valor3^(-0.5))*exp(0.25*(kimmelz(t,valor3)^2)*(1-((k3+eta3)/(2*kimmelb))))*kimmelpifinal(t,valor3)*(kimmelz(t,valor3)/2)*(1-((k3+eta3)/(2*kimmelb)));
derivadakimmelfinalauxiliar4= @(t,valor3) (((4*valor3)/(sigma3^2))^(0.25-((theta3*k3)/(sigma3^2))))*exp((0.5*(k3+eta3)-((0.5*kimmelb)/(sigma3^2)))*valor3)*exp((-0.5*kimmelb-kimmeld)*t)*(sqrt(2*kimmelb)/(sigma3))*exp(-kimmelb*t)*(valor3^(-0.5))*exp(0.25*(kimmelz(t,valor3)^2)*(1-((k3+eta3)/(2*kimmelb))))*derivadakimmelinicial(t,valor3);

derivadakimmelpifinal = @(t,valor3) derivadakimmelfinalauxiliar1(t,valor3)+derivadakimmelfinalauxiliar2(t,valor3)+derivadakimmelfinalauxiliar3(t,valor3)+derivadakimmelfinalauxiliar4(t,valor3);

R = eye(ncol);
for i = 1:ncol
    R(i,i) = sigmai(i)^2;
end
dt = 1/12; %monthly data
initx = [theta3 ,theta5]';
initV = [sigma3^2*theta3/(2*k3),0;0,sigma5^2*theta5/(2*k5)];


% parameter setting for transition equation
C = [theta3*(1-exp(-k3*dt)),theta5*(1-exp(-k5*dt))]' ;
F = [exp(-k3*dt),0;0,exp(-k5*dt)];


% parameter setting for measurement equation
A = zeros(ncol, 1);
H = zeros(ncol,2);
 
%EKF
EY = initx;
VARY = initV;
ll = zeros(nrow,1); %log-likelihood

dates = [];
Aaux = 0;
Haux = 0;
Haux2 = 0;
for i=1:ncup
   if datescupmatlab(i) > datesobsmatlab(1)
       dates = [dates,datescupmatlab(i)];
   end

end

numdatesaux = size(dates);
numdates = numdatesaux(2);

Ezaux=0;
%Agora e ao contrário das treasury yields os A e H's vão ser diferentes
%para cada observação uma vez que estamos a actualizar diferente numero de
%cupões
for j=1:numdates
Aaux=Aaux+cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1))*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(2))+cupao*EY(2)*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1))*AffineBeta5((ceil((dates(j)-datesobsmatlab(1))/31))/12)*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(2))-cupao*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(2))*derivadakimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1))*EY(1);
Haux = Haux + cupao*derivadakimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1))*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(2)); 
Haux2 = Haux2 - AffineBeta5((ceil((dates(j)-datesobsmatlab(1))/31))/12)*cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1))*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(2));

%Ezaux = Ezaux+cupao*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(2))*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1));
end
A=Aaux+(1+AffineBeta5((ceil((dates(end)-datesobsmatlab(1))/31))/12)*EY(2))*kimmelpifinal((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(1))*f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(2))-cupao*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(2))*derivadakimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1))*EY(1);
H(1) = Haux + derivadakimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1))*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(2));
H(2) = Haux2 - AffineBeta5((ceil((dates(end)-datesobsmatlab(1))/31))/12)*kimmelpifinal((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(1))*f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(2));

% Aaux = 0;
Haux = 0;
Haux2 = 0;
%dates =[];




 EZ = A+H*EY;
% EZ = Ezaux + f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(2))*kimmelpifinal((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(1));
Aaux=0; %retirar
%Ezaux=0;
dates =[];%retirar
VARZ = H*VARY*H'+R;
Erro = Z(1,:)-EZ';
Erro=Erro';
K = VARY*H'*inv(VARZ);
EY = EY +K*Erro;

if EY(1) < 0
        EY(1) =0;
        
    else
        EY(1) =EY(1);
end
    
    if EY(2) < 0
        EY(2) =0;
        
    else
        EY(2) =EY(2);
    end

VARY = (eye(2)-K*H)*VARY;
ll(1) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;

estado3(1)=EY(1);
estado5(1)=EY(2);



for i = 2:nrow
    
    FUTUROY = C+F*EY;
    
    if FUTUROY(1) < 0
        FUTUROY(1) =0;
        
    else
        FUTUROY(1) =FUTUROY(1);
    end
    
    if FUTUROY(2) < 0
        FUTUROY(2) =0;
        
    else
        FUTUROY(2) =FUTUROY(2);
    end
      
for o=1:ncup
   if datescupmatlab(o) > datesobsmatlab(i)
       dates = [dates,datescupmatlab(o)];
   end

end

numdatesaux = size(dates);
numdates = numdatesaux(2);

for j=1:numdates
Aaux=Aaux+cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(2))+cupao*FUTUROY(2)*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))*AffineBeta5((ceil((dates(j)-datesobsmatlab(i))/31))/12)*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(2))-cupao*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(2))*derivadakimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))*FUTUROY(1);
 Haux = Haux + cupao*derivadakimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(2)); 
Haux2 = Haux2 - AffineBeta5((ceil((dates(j)-datesobsmatlab(i))/31))/12)*cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(2));
%Ezaux = Ezaux + cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(2));
end
A=Aaux+(1+AffineBeta5((ceil((dates(end)-datesobsmatlab(i))/31))/12)*FUTUROY(2))*kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY(1))*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY(2))-cupao*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(2))*derivadakimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))*FUTUROY(1);
H(1) = Haux + derivadakimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(2));
H(2) = Haux2 - AffineBeta5((ceil((dates(end)-datesobsmatlab(i))/31))/12)*kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY(1))*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY(2));

Haux = 0;
Haux2 = 0;
    
     
 Q1 = ((sigma3*sigma3)/k3)*(1-exp(-k3*dt))*((theta3/2)*(1-exp(-k3*dt))+exp(-k3*dt)*FUTUROY(1));
 Q2 = ((sigma5*sigma5)/k5)*(1-exp(-k5*dt))*((theta5/2)*(1-exp(-k5*dt))+exp(-k5*dt)*FUTUROY(2));
 Q = [Q1,0;0,Q2];   


VARY = VARY-F*VARY*F'+Q;
EZ = A+H*FUTUROY;
%EZ = Ezaux + kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY(1))*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY(2));
Aaux =0; %retirar
%Ezaux =0;
dates =[]; %retirar
    VARZ = H*VARY*H'+R;
    Erro = Z(i,:)-EZ';
    Erro=Erro';
    K = VARY*H'*inv(VARZ);
    FUTUROY = FUTUROY +K*Erro;
    
     if FUTUROY(1) < 0
        FUTUROY(1) =0;
        
    else
        FUTUROY(1) =FUTUROY(1);
    end
    
    if FUTUROY(2) < 0
        FUTUROY(2) =0;
        
    else
        FUTUROY(2) =FUTUROY(2);
    end
    
    
    VARY = (eye(2)-K*H)*VARY;
    EY = FUTUROY;
    
    
    estado3(i) = EY(1);
    estado5(i) = EY(2);
    

    
    ll(i) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;
    
   
end

estado = [estado3' estado5'];
csvwrite('makewhole.csv',estado)

sumll = -sum(ll)
 
end



