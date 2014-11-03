function sumll = CIRCall(para,input, datescup,parametros,cupao) %calculate log likelihood
format long
warning('off')
estado1 = input(:,3);
estado2 = input(:,4);
estado3 = input(:,5);

Z = input(:,2)/100;

datesobs = input(:,1);
datesobsmatlab = datesobs + datenum('30-Dec-1899');
datescupmatlab = datescup + datenum('30-Dec-1899');
[nrow,ncol] = size(datesobsmatlab); %numero de observações e de obrigacoes

%parametros a estimar
theta4 = para(1); k4 = para(2); sigma4 = para(3); eta4 = para(4);
csi0 = para(5); csi3 = para(6); 
sigmai = para(7);

%paramteros já estimados
theta1 = parametros(1,1) ; theta2 = parametros(1,2); theta3 = parametros(1,3); 
k1 = parametros(2,1) ; k2 = parametros(2,2); k3 = parametros(2,3);
sigma1 = parametros(3,1) ; sigma2 = parametros(3,2); sigma3 = parametros(3,3);
eta1 = parametros(4,1); eta2 = parametros(4,2); eta3 = parametros(4,3);
alfa0 = parametros(1,4);gama0 = parametros(2,4); betad1 = parametros(3,4); betad2 = parametros(4,4); 


%n de cupoes
ncupaux = size(datescup);
ncup = ncupaux(1);

h1=sqrt(((k1+eta1)^2)+2*(1+betad1)*sigma1*sigma1);
h2=sqrt(((k2+eta2)^2)+2*(1+betad2)*sigma2*sigma2);
h4=sqrt(((k4+eta4)^2)+2*sigma4*sigma4);

AffineAlpha1 = @(tau) ((2*k1*theta1)/(sigma1*sigma1))*log((2*h1*exp(0.5*(k1+eta1+h1)*tau))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau)));
AffineAlpha2 = @(tau) ((2*k2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(k2+eta2+h2)*tau))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau)));
AffineAlpha4 = @(tau) ((2*k4*theta4)/(sigma4*sigma4))*log((2*h4*exp(0.5*(k4+eta4+h4)*tau))/(h4-(k4+eta4)+(k4+eta4+h4)*exp(h4*tau)));

AffineBeta1 = @(tau) (2*(1+betad1)*(exp(h1*tau)-1))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau));
AffineBeta2 = @(tau) (2*(1+betad2)*(exp(h2*tau)-1))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau));
AffineBeta4 = @(tau) (2*(exp(h4*tau)-1))/(h4-(k4+eta4)+(k4+eta4+h4)*exp(h4*tau));

f = @(t,valor1,valor2) exp(-(alfa0+gama0)*t+AffineAlpha1(t)+AffineAlpha2(t)-AffineBeta2(t)*valor2-AffineBeta1(t)*valor1);
g = @(t,valor4) exp(-csi0*t+AffineAlpha4(t)-AffineBeta4(t)*valor4);

kimmela = ((2*(theta3^2)*(k3^2))/(sigma3^4))-((2*theta3*k3)/(sigma3^2))+((4*cupao*csi3)/(sigma3^2))+(3/8);
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

R = eye(ncol);
for i = 1:ncol
    R(i,i) = sigmai(i)^2;
end
dt = 1/12; %monthly data
initx = theta4;
initV = sigma4^2*theta4/(2*k4);


% parameter setting for transition equation
C = theta4*(1-exp(-k4*dt)) ;
F = exp(-k4*dt);


% parameter setting for measurement equation
A = zeros(ncol, 1);
H = zeros(ncol,1);
 
%EKF
EY = initx;
VARY = initV;
ll = zeros(nrow,1); %log-likelihood

dates = [];
Aaux = 0;
Haux = 0;

for i=1:ncup
   if datescupmatlab(i) > datesobsmatlab(1)
       dates = [dates,datescupmatlab(i)];
   end

end

numdatesaux = size(dates);
numdates = numdatesaux(2);


%Agora e ao contrário das treasury yields os A e H's vão ser diferentes
%para cada observação uma vez que estamos a actualizar diferente numero de
%cupões
for j=1:numdates
Aaux=Aaux+cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado3(1))*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1))+cupao*EY(1)*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado3(1))*AffineBeta4((ceil((dates(j)-datesobsmatlab(1))/31))/12)*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1));
 Haux = Haux - AffineBeta4((ceil((dates(j)-datesobsmatlab(1))/31))/12)*cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado3(1))*f((ceil((dates(j)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(j)-datesobsmatlab(1))/31))/12,EY(1));
end
A=Aaux+(1+AffineBeta4((ceil((dates(end)-datesobsmatlab(1))/31))/12)*EY(1))*kimmelpifinal((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado3(1))*f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(1));
H(1) = Haux - AffineBeta4((ceil((dates(end)-datesobsmatlab(1))/31))/12)*kimmelpifinal((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado3(1))*f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY(1));
% Aaux = 0;
Haux = 0;
%dates =[];




 EZ = A+H*EY;
% EZ = Aaux + f((ceil((dates(end)-datesobsmatlab(1))/31))/12,estado1(1),estado2(1),estado3(1))*g((ceil((dates(end)-datesobsmatlab(1))/31))/12,EY);
Aaux=0; %retirar
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
    
    

VARY = (eye(1)-K*H)*VARY;
ll(1) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;

estado4(1)=EY(1);



for i = 2:nrow
    
    FUTUROY = C+F*EY;
    
    if FUTUROY(1) < 0
        FUTUROY(1) =0;
        
    else
        FUTUROY(1) =FUTUROY(1);
    end
    
    
for o=1:ncup
   if datescupmatlab(o) > datesobsmatlab(i)
       dates = [dates,datescupmatlab(o)];
   end

end

numdatesaux = size(dates);
numdates = numdatesaux(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:numdates
Aaux=Aaux+cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado3(i))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1))+cupao*FUTUROY(1)*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado3(i))*AffineBeta4((ceil((dates(j)-datesobsmatlab(i))/31))/12)*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1));
Haux = Haux - AffineBeta4((ceil((dates(j)-datesobsmatlab(i))/31))/12)*cupao*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado3(i))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,FUTUROY(1));
 
end
A=Aaux+(1+AffineBeta4((ceil((dates(end)-datesobsmatlab(i))/31))/12)*FUTUROY(1))*kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado3(i))*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY(1));
H(1) = Haux - AffineBeta4((ceil((dates(end)-datesobsmatlab(i))/31))/12)*kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado3(i))*f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,FUTUROY(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Haux = 0;
    
     

 Q = ((sigma4*sigma4)/k4)*(1-exp(-k4*dt))*((theta4/2)*(1-exp(-k4*dt))+exp(-k4*dt)*FUTUROY(1));   


VARY = VARY-F*VARY*F'+Q;
EZ = A+H*FUTUROY;
Aaux =0; %retirar
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
    
    
    
    VARY = (eye(1)-K*H)*VARY;
    EY = FUTUROY;
    
    
  
    estado4(i) = EY(1);
    

    
    ll(i) = -(ncol/2)*log(2*pi)-0.5*log(det(VARZ))-0.5*Erro'*inv(VARZ)*Erro;
    
   
end

estado = [estado4'];
csvwrite('call.csv',estado)

sumll = -sum(ll)
 
end



