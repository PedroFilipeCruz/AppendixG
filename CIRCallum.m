function [para, sumll] = CIRCallum
%ficheiro que na 1ª coluna tem datas mensais das observações. Depois tem na
%2ª um valor numerico dado pelo excel a essa data, na 3ª tem os preços, na
%4ª e 5ª os valores dos vectores de estado que formam a risk-free short
%rate para as datas das observações (mas no total tenho para entre
%Jan82 e Jun2004) 
input = xlsread('EI406362.xlsx'); 

%ficheiro que na 1ª coluna tem os meses de cupão e na seguinte valor
%númerico dado pelo excel a essa data
datescup = xlsread('cupoesEI406362.xlsx');

%primeira coluna tem os parametros do primeiro factor e na segunda o mesmo. A
%3º na primeira linha tem o alfa0 e o resto da linha tem 0's. 
parametros = xlsread('parametros3cir.xlsx');

%valor do cupão que a bond paga. Se é 10% semi anual paga 5% a cada 6
%meses. Para o código funcionar tenho que dar os 5% e não os 10%. 
cupao = 0.05/2;

%valor do make-whole premium

para0 = [0.02, 0.01, 0.0016, -1,-0.9,0.000000000008,0.1, 0.01]; %valores iniciais para começar a minimização
options = struct('FinDiffRelStep',1,'MaxFunEvals',100);
optionss = gaoptimset('TimeLimit',500);


%optimização
[x, fval] = fmincon(@(x) CIRCall(x,input, datescup,parametros,cupao), para0,[],[],[],[],[0.001,0.001,0.001, -1.5,-1.5,0.0000000000001,0.000000000001, 0.00001],[0.1,0.1,0.1,0,1,1,1,0.01],@inequa,options);
%[x, fval] = ga(@(x) CIRCall(x,input, datescup,parametros,cupao),8,[],[],[],[],[0.0000001,0.0000001,0.00000001,-1,-1,0,0,0.000001],[1,1,1,0,1,1,1,0.01],@inequa,optionss)

para = x;
sumll = fval;

%ler todos os vectores de estado incluindo o que agora já temos para
%default. 
M = csvread('call.csv');
estado4 = M(:,1);
estado1 = input(:,3);
estado2 = input(:,4);
estado3 = input(:,5);


%passar as datas que estão com valores do excel para a forma como o matlab
%lida com datas. 
datesobs = input(:,1);
datesobsmatlab = datesobs + datenum('30-Dec-1899');
datescupmatlab = datescup + datenum('30-Dec-1899');
[nrow,ncol] = size(datesobsmatlab);

%obter o numero de cupoes a ser pagos considerando a data da 1ª observação
ncupaux = size(datescup);
ncup = ncupaux(1);

%valores de preços(observações) considerando FV=1
Z = input(:,2)/100;

%parâmetros obtidos já depois da minimização
theta4 = para(1); k4 = para(2); sigma4 = para(3); eta4 = para(4);
csi0 = para(5); csi3 = para(6); 
csi1 = para(7); 
sigmai = para(8);


%paramteros já estimados
theta1 = parametros(1,1) ; theta2 = parametros(1,2); theta3 = parametros(1,3); 
k1 = parametros(2,1) ; k2 = parametros(2,2); k3 = parametros(2,3);
sigma1 = parametros(3,1) ; sigma2 = parametros(3,2); sigma3 = parametros(3,3);
eta1 = parametros(4,1); eta2 = parametros(4,2); eta3 = parametros(4,3);
alfa0 = parametros(1,4);gama0 = parametros(2,4); betad1 = parametros(3,4); betad2 = parametros(4,4); 

%funções necessárias para escrever o preço das bonds

h2=sqrt(((k2+eta2)^2)+2*(1+betad2)*sigma2*sigma2);
h4=sqrt(((k4+eta4)^2)+2*sigma4*sigma4);

AffineAlpha2 = @(tau) ((2*k2*theta2)/(sigma2*sigma2))*log((2*h2*exp(0.5*(k2+eta2+h2)*tau))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau)));
AffineAlpha4 = @(tau) ((2*k4*theta4)/(sigma4*sigma4))*log((2*h4*exp(0.5*(k4+eta4+h4)*tau))/(h4-(k4+eta4)+(k4+eta4+h4)*exp(h4*tau)));

AffineBeta2 = @(tau) (2*(1+betad2)*(exp(h2*tau)-1))/(h2-(k2+eta2)+(k2+eta2+h2)*exp(h2*tau));
AffineBeta4 = @(tau) (2*(exp(h4*tau)-1))/(h4-(k4+eta4)+(k4+eta4+h4)*exp(h4*tau));

f = @(t,valor2) exp(-(alfa0+gama0)*t+AffineAlpha2(t)-AffineBeta2(t)*valor2);
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


kimmela1 = ((2*(theta1^2)*(k1^2))/(sigma1^4))-((2*theta1*k1)/(sigma1^2))+((4*cupao*csi1)/(sigma1^2))+(3/8);
kimmelb1 = 0.5*sqrt(((k1+eta1)^2)+2*(1+betad1)*(sigma1^2));
kimmeld1 = (-1)*((theta1*k1*(k1+eta1))/(sigma1^2));
kimmelalfa1 = ((2*theta1*k1)/(sigma1^2))-0.5;
kimmelgama1 = 0.5*(1-sqrt(1+8*kimmela1));
kimmeldelta1 = @(t) 1-exp(-2*kimmelb1*t);
kimmely1 = @(valor1) (2*sqrt(valor1))/(sigma1);
kimmelz1 = @(t,valor1) sqrt(2*kimmelb1)*exp(-kimmelb1*t)*kimmely1(valor1);

auxiliarpi11= @(t,valor1) ((kimmelalfa1-kimmelgama1)*(kimmelalfa1+kimmelgama1-1)/(2*(kimmelz1(t,valor1)^2)))+((2*kimmelalfa1+1)/4)*(1-((k1+eta1)/(2*kimmelb1)))+(kimmelz1(t,valor1)^2)*(1/8)*((1-((k1+eta1)/(2*kimmelb1)))^2);
auxiliarpi21 = @(t,valor1) ((kimmelalfa1-kimmelgama1)*(kimmelalfa1-kimmelgama1-2)*(kimmelalfa1+kimmelgama1-1)*(kimmelalfa1+kimmelgama1-3))/(4*(kimmelz1(t,valor1)^4));
auxiliarpi31 = @(t,valor1) ((2*kimmelalfa1-1)*(kimmelalfa1-kimmelgama1)*(kimmelalfa1+kimmelgama1-1)*(1-((k1+eta1)/(2*kimmelb1))))/(4*kimmelz1(t,valor1)^2);
auxiliarpi41 = @(t,valor1) ((1/16)*(2*kimmelalfa1+3)*(2*kimmelalfa1+1)+(1/8)*(kimmelalfa1-kimmelgama1)*(kimmelalfa1+kimmelgama1-1))*((1-((k1+eta1)/(2*kimmelb1)))^2);
auxiliarpi51 = @(t,valor1) (1/16)*(2*kimmelalfa1+3)*(kimmelz1(t,valor1)^2)*((1-((k1+eta1)/(2*kimmelb1)))^3)+(1/64)*(kimmelz1(t,valor1)^4)*((1-((k1+eta1)/(2*kimmelb1)))^4);

kimmelpi1 = @(t,valor1) 1+kimmeldelta1(t)*auxiliarpi11(t,valor1)+0.5*(kimmeldelta1(t)^2)*(auxiliarpi21(t,valor1)+auxiliarpi31(t,valor1)+auxiliarpi41(t,valor1)+auxiliarpi51(t,valor1));
kimmelwum1 = @(t,valor1) ((kimmelz1(t,valor1)/(sqrt(2*kimmelb1)))^(kimmelalfa1-kimmelgama1))*exp((1/4)*(1-((k1+eta1)/(2*kimmelb1)))*(kimmelz1(t,valor1)^2))*kimmelpi1(t,valor1);
kimmelh1 = @(t,valor1) ((kimmelz1(t,valor1)/(sqrt(2*kimmelb1)))^kimmelgama1)*exp(-0.5*kimmelb1*(kimmely1(valor1)^2)-(0.5*kimmelb1+kimmeld1)*t)*kimmelwum1(t,valor1);

kimmelpifinal1 = @(t,valor1) ((4*valor1/(sigma1^2))^(0.25-((theta1*k1)/(sigma1^2))))*exp((valor1*(k1+eta1))/(sigma1^2))*kimmelh1(t,valor1);


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





Precoaux = 0;
dates =[];  
for i=1:nrow %para cada observação
  for o=1:ncup %para a observação que temos 'seleccionada' vemos quais são as datas de cupão a partir daí
   if datescupmatlab(o) > datesobsmatlab(i)
       dates = [dates,datescupmatlab(o)];
   end

end

numdatesaux = size(dates);
numdates = numdatesaux(2); %numero de cupões que faltam para essa observação


for j=1:numdates %obtem-se os valores descontados dos cupões
 Precoaux = Precoaux + cupao*kimmelpifinal1((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i))*kimmelpifinal((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado3(i))*f((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado2(i))*g((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado4(i));
end
%juntam-se o valor descontado do FV=1
Preco(i)=Precoaux+f((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado2(i))*kimmelpifinal1((ceil((dates(j)-datesobsmatlab(i))/31))/12,estado1(i))*kimmelpifinal((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado3(i))*g((ceil((dates(end)-datesobsmatlab(i))/31))/12,estado4(i));

Precoaux = 0;
dates =[];  

end 

%temos observações mensais
dt = 1/12;
deltat = 0:dt:dt*(nrow-1);

%yield to maturity das observações e do fit. Está duas vezes cupão porque é
%com o valor do cupão anual que a função funciona
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

%plots dos preços e fit de preços e das yields e fit das yields
figure(1)
plot(deltat,Z)
hold on
plot(deltat,Preco,'Color','g')

figure(2)
plot(deltat,Yieldpedro)
hold on
plot(deltat,Yieldfitpedro,'Color','g')

%para calcular o RMSE nas yield to maturity
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

%impor condição de feller
function [c,ceq]=inequa(x)


 c(1) =  x(3)^2 - 2*x(1)*x(2);
 ceq = [];   


end

    


    