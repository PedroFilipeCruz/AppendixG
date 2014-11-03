function [para, sumll] = CIRDefaultumtent
format long
%1st columns have the dates of the observations. In the
%2nd the numerical value given by excel to that data, in the 3rd the bond
%prices, and in  the 4th, 5th and 6th the value of Y1 and Y2 for those
%dates
input = xlsread('boadefaultprecos.xlsx'); 

%1st columns the dates of coupon payments and in the 2nd the numerical
%value given by excel
datescup = xlsread('boacupoesdefault.xlsx');

%1st column have theta, k, sigma and eta for Y1 and in the second the same for Y2. 
%In the third column, the first value is alfa0 and the rest are zeros
parametros = xlsread('parametros2cir.xlsx');

%coupon value of the bond
cupao = 0.0605/2;

%trial parameter vector if we want to use fmin con
para0 = [0.023, 0.3115, 0.119, -0.201, 0.627,-0.627,-0.611, 0.001]; %valores iniciais para começar a minimização
optionss = struct('MaxFunEvals',300); %for fmincon
options = gaoptimset('TimeLimit',400); %for the genetic algorithm

%optimization
[x, fval] = fmincon(@(x) CIRDefaultent(x,input, datescup,parametros,cupao), para0,[],[],[],[],[0.02,0.31,0.11, -0.30,0.62,-0.68,-0.62,0.00000001],[0.03,0.32,0.12,-0.29,0.63,-0.67,-0.61,0.01],@inequa,optionss);
%[x, fval] = ga(@(x) CIRDefaultent(x,input, datescup,parametros,cupao),8,[],[],[],[],[0.001,0.01,0.001, -0.999,0,-0.999,-0.999,0.0000000001],[0.5,1,0.5,0,1,0.7,0.7,0.01],@inequa,options);
para = x;
sumll = fval;

%state vectors
M = csvread('default.csv');
estado3 = M(:);
estado1 = input(:,3);
estado2 = input(:,4);



datesobs = input(:,1);
datesobsmatlab = datesobs + datenum('30-Dec-1899');
datescupmatlab = datescup + datenum('30-Dec-1899');
[nrow,ncol] = size(datesobsmatlab);


ncupaux = size(datescup);
ncup = ncupaux(1);

%observation prices
Z = input(:,2)/100;

%parameters
theta3 = para(1); k3 = para(2); sigma3 = para(3); eta3 = para(4); gama0 = para(5);
beta1 = para(6); beta2 = para(7);

theta1 = parametros(1,1) ; theta2 = parametros(1,2);
k1 = parametros(2,1) ; k2 = parametros(2,2); 
sigma1 = parametros(3,1) ; sigma2 = parametros(3,2); 
eta1 = parametros(4,1); eta2 = parametros(4,2); 
alfa0 = parametros(1,3);

%functions nedded to obtain bond prices

h1=sqrt(((k1+eta1)^2)+2*(1+beta1)*sigma1*sigma1);
h2=sqrt(((k2+eta2)^2)+2*(1+beta2)*sigma2*sigma2);

h3=sqrt(((k3+eta3)^2)+2*sigma3*sigma3);

AffineAlpha1 = @(tau) ((2*k1*theta1)/(sigma1*sigma1))*log((2*h1*exp(0.5*(k1+eta1+h1)*tau))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau)));
AffineAlpha2 = @(tau) ((2*k2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(k2+eta2+h2)*tau))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau)));
AffineAlpha3 = @(tau) ((2*k3*theta3)/(sigma3*sigma3))*log((2*h3*exp(0.5*(k3+eta3+h3)*tau))/(h3-(k3+eta3)+(k3+eta3+h3)*exp(h3*tau)));


AffineBeta1 = @(tau) (2*(1+beta1)*(exp(h1*tau)-1))/(h1-(k1+eta1)+(k1+eta1+h1)*exp(h1*tau));
AffineBeta2 = @(tau) (2*(1+beta2)*(exp(h2*tau)-1))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau));

AffineBeta3 = @(tau) (2*(exp(h3*tau)-1))/(h3-(k3+eta3)+(k3+eta3+h3)*exp(h3*tau));

f = @(t,valor1,valor2) exp(-alfa0*t+AffineAlpha1(t)+AffineAlpha2(t)-AffineBeta1(t)*valor1-AffineBeta2(t)*valor2);
g = @(t,valor3) exp(-gama0*t+AffineAlpha3(t)-AffineBeta3(t)*valor3);



Precoaux = 0;
dates =[];  
for i=1:nrow 
  for o=1:ncup 
   if datescupmatlab(o) > datesobsmatlab(i)
       dates = [dates,datescupmatlab(o)];
   end

end

numdatesaux = size(dates);
numdates = numdatesaux(2);

%price of the bond by the method for each observation
for j=1:numdates 
 Precoaux = Precoaux + cupao*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado3(i));
end

Preco(i)=Precoaux+f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado1(i),estado2(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado3(i));

Precoaux = 0;
dates =[];  

end 


dt = 1/12;
deltat = 0:dt:dt*(nrow-1);

%implicit yield-to.maturity from observational prices and from the prices given by the
%method
cashflow = [];
for i=1:nrow %para cada observação
  for o=1:ncup %para a observação que temos 'seleccionada' vemos quais são as datas de cupão a partir daí
   if datescupmatlab(o) > datesobsmatlab(i)
       dates = [dates,datescupmatlab(o)];
       cashflow = [cashflow, cupao];
   end
  end
  cashflow(end) = cashflow(end)+1;
Yieldpedro(i) = ytmpedro(2*cupao,(ceil((dates(end)-datesobsmatlab(i))/31))/12,Z(i),1);
Yieldfitpedro(i) = ytmpedro(2*cupao,(ceil((dates(end)-datesobsmatlab(i))/31))/12,Preco(i),1);
dates =[]; 
cashflow = [];
end




figure(1)
plot(deltat,100*Z)
hold on
plot(deltat,100*Preco,'Color','g')

figure(2)
plot(deltat,Yieldpedro)
hold on
plot(deltat,Yieldfitpedro,'Color','g')


%RMSE
diferencaspreco=[];
diferencaspedro=[];
for i=1:nrow
      AuxiliarDiferencaspedro = (Yieldpedro(i)-Yieldfitpedro(i)).^2;
    diferencaspedro=[diferencaspedro,AuxiliarDiferencaspedro];
      AuxiliarDiferencaspreco = (Z(i)-Preco(i)).^2;
    diferencaspreco=[diferencaspreco,AuxiliarDiferencaspreco];

end

RMSEpreco = 10000*sqrt(sum(diferencaspreco)/(nrow))
RMSEpedro = 10000*sqrt(sum(diferencaspedro)/(nrow))
likelihood=fval

end

%Feller Condition
function [c,ceq]=inequa(x)


 c =  x(3)^2 - 2*x(1)*x(2);
 ceq = [];   


end
    


    