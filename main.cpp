// System includes
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <libalglib/ap.h>
#include <libalglib/specialfunctions.h>

int A = 1;  //коэффициент при правой части
int C = 0;  //левый край отрезка по x
int D = 1;  //правый край отрезка по x
int T1 = 0;  //начальное время t
int T2 = 1;  //конечное время t
int n = 0;  //номер слоя по времени

bool flag = false;  //Выбор задачи(false - основная задача, true - задача на разрывную функцию)

double alpha = 0.5;  //задайте альфа на полуинтервале (0;1]
double dx = 0.1;  //шаг по х
double dt = 0.001 * pow(dx, 3.0 / alpha) / A;
double s = 0.1;  //самая первая ступенька (тесная зависимость с alpha, по какой формуле?)

size_t stepAmountX=(size_t)((D-C)/dx+1);		//кол-во узлов сеточной функции
int stepAmountT=(int)((T2-T1)/dt);		//полные временные слои
int myStep=100;					//число временных слоев необходимых для вычисления
int dodo=2;						//заглушка для явной схемы, показывающая с КАКОГО временного слоя начинать неявную схему

int threadNum = 1;
//------------------------------------------------------------------------------------------------
//массивы для прогонки
double* fiveUp=(double*)calloc(stepAmountX,sizeof(double));
double* fiveDown=(double*)calloc(stepAmountX,sizeof(double));
double* A1=(double*)calloc(stepAmountX,sizeof(double));//поддиагональ
double* B1=(double*)calloc(stepAmountX,sizeof(double));//диагональ
double* C1=(double*)calloc(stepAmountX,sizeof(double));//наддиагональ
double* D1=(double*)calloc(stepAmountX,sizeof(double));//здесь правая часть системы
double* X =(double*)calloc(stepAmountX,sizeof(double));//метод прогонки записывает решение
double* ksi=(double*)calloc((stepAmountX+1),sizeof(double));//прогоночные коэфициенты
double* eta=(double*)calloc((stepAmountX+1),sizeof(double));//прогоночные коэфициенты
double* dzeta=(double*)calloc((stepAmountX+1),sizeof(double));//прогоночные коэфициенты
//------------------------------------------------------------------------------------------------
double* func = (double*)calloc(stepAmountX, sizeof(double));		//массив значений косинуса
double* funcDer = (double*)calloc(stepAmountX, sizeof(double));		//1 производная косинуса
double* funcDer2 = (double*)calloc(stepAmountX, sizeof(double));	//2 производная косинуса
double* funcDer3 = (double*)calloc(stepAmountX, sizeof(double));	//3 производная 

double** U=(double**)calloc(myStep, sizeof(double*));	//временный массив для значений V
double** Fin=(double**)calloc(myStep, sizeof(double*));	//конечный массив результатов (нормированный)
double* Us=(double*)calloc(stepAmountX, sizeof(double));

double* Coeff = (double*)calloc(myStep, sizeof(double));

//-------------------------------Необходимые функции--------------------------------------------
//метод прогонки
//sloy - номер временного слоя на котором выполняется прогонка
void Shuttle(int sloy){
	int N = stepAmountX - 1;
	A1[0] = 0;
	C1[N-1] = 0;
	ksi[0] = 0;
	eta[0] = 0;
	for (int i = 0; i < N; i++)
	{
		ksi[i+1] = C1[i] / (-B1[i] - A1[i] * ksi[i]);
		eta[i+1] = (A1[i] * eta[i] - D1[i]) / (-B1[i] - A1[i] * ksi[i]);
	}
	X[N-1] = -eta[N];
	for(int i = N - 2; i > -1; i--)
	{
		X[i] = ksi[i+1] * X[i+1] - eta[i+1];
	}
	#pragma omp parallel for
		for (int j = 1; j < N; j++)
		{
			U[sloy][j] = -X[j];
		}
}
//-----------------------------------------------------------------------------
//гамма функция
//input - аргумент гамма функции
double gamma(double input)
{
	return alglib::gammafunction(input);
}
//-----------------------------------------------------------------------------------------
//записывает результаты в файл, с нормировкой или без неё
void toFile(int t,double d)
{
	double gam = gamma(alpha);
	double shag=0.;
	std::fstream str;
	str.open("out.txt", std::ios::out);
	//for(int n=0;n<t;n+=49)	//для вывода малого количества графиков (чтобы не захламлять графиками картинку)
	for(int n=0;n<t;n++)
	{
		str<<"#"<<n<<std::endl;	//для gnuplot
		shag=0;
		for(int i=0;i<stepAmountX;i++)
		{
			if(flag)
			{
				str<< n * dt << "\t" <<shag<<"\t"<<Fin[n][i]/gam<<endl;	//нормировка по Гамма(альфа)
			}else{
				str<< n * dt << "\t" <<shag<<"\t"<<Fin[n][i]<<endl;
			}
			shag=shag+d;
		}
		str<<endl;
	}
	str.close();
}
//------------------------------------------------------------------------------------------
//гипергеометрическая функция
double hyperGeometric(double a, double b, double c, double z)
{
	if(z==0) return 1.;//подстановка 0

	double res=1., res2=1.;
	int hyperIter=20;		//кол-во слагаемых для гипергеометрической функции
	for(int i=1;i<hyperIter;i++)
	{
		res2=1;
		for(int l=0;l<i-1;l++)
		{
			res2*=(a+l)*(b+l)/((1+l)*(c+l));//par dont need
		}
		res+=pow(z,i)*res2;
	}
	return res;
}
//----------------------------------------------------------------------------------------
//функция - начальное условие
double fCos(double x)
{
	if(flag==false)
	{
		return x*x*(1-x);
	}
	else
	{
		double res=0;
		if ( x>=0 && x<0.5 )  res=1;
		if ( x==0.5 )		  res=0.5;
		if ( x>0.5 && x<=1 )  res=0;
		return res;
	}
}
//-----------------------------------------------------------------------------------------
//функция - граничное условие
double g(double t)
{
	return 0.;
	//return gamma(alpha);
	//return 1.;
}
//-------------------------------------------------------------------------------------
//центральная производная
//sloy - временной слой на котором берется производная
//nextTochka - X(i+1)
//currTochka - X(i-1)
double centrDeriv(int sloy, int nextTochka, int currTochka)
{
	return (U[sloy][nextTochka]-U[sloy][currTochka])/(2*dx);
}
//-----------------------------------------------------------------------------------------
// dt/()
//вычисляет интеграл
//podstanovka - Tn
//up - T(j+1)
//down - Tj
double newIntegral1(double podstanovka, int check)
{
	double ret=0;
	if (podstanovka>check)
	{
		ret=(1/(alpha-1))*(pow(s+podstanovka*dt-check*dt,1-alpha)-pow(s+podstanovka*dt-(check-1)*dt,1-alpha));
	}
	return ret;
}
//----------------------------------------------------------------------------------------
// t*dt/()
//podstanovka - Tn
//up - T(j+1)
//down - Tj
double newIntegral2(double podstanovka, int check)
{
	double ret=0;
	if (podstanovka>check)
	{
		ret=((1.)/((alpha-2)*(alpha-1)))*(pow(s+podstanovka*dt-check*dt,1-alpha)*(alpha*check*dt-s-podstanovka*dt-check*dt)-pow(s+podstanovka*dt-(check-1)*dt,1-alpha)*(alpha*(check-1)*dt-s-podstanovka*dt-(check-1)*dt));
	}
	return ret;
}
//---------------------------------------------------------------------------------------
//подсчет дробного интеграла на n слое
//down - номер временного слоя
//order -  порядок дробного интеграла
//step - номер шага по Ox
double fracInt(double down, double order, int step)
{
	//случай выхода за границу сетки
	if(step>=stepAmountX) return 0.;

	double result=0;
	for (int j=1;j<down;j++)
	{
		result+=newIntegral1(n-1,j+1)*(U[j][step]+((s+(j-1)*dt)/dt)*(-U[j+1][step]+U[j][step])) + 1/(dt)*newIntegral2(n-1,j+1)*(U[j+1][step]-U[j][step]);
	}
	return result/gamma(order);
}
//------------------MAIN---------------------------------------------------------
int main()
{
	omp_set_dynamic(0);		// запретить библиотеке openmp менять число потоков во время исполнения программы
	omp_set_num_threads(threadNum); 
	double start,end;	//переменные для подсчета потраченного времени
	start=omp_get_wtime();
	fstream temp;
	temp.open("temp.txt", ios::out);
	cout.precision(15);
	cout<<"Number of threads: "<<threadNum<<endl;
	cout<<"dt: "<<dt<<endl;
	cout<<"dx: "<<dx<<endl;
	cout<<"Step amount over x: "<<stepAmountX<<endl;
	cout<<"Step amount over t: "<<stepAmountT<<endl;
	cout<<"How much to count: "<<myStep<<endl;
	cout<<"Task type: "<<flag<<endl;
	double shag=0.,sum1=0.,sum2=0.,sum3=0.,sum4=0.,sum5=0.,sum6=0.,sum7=0.,sum8=0.,allSum=0.,sumGroup=0,sumGroup2=0;
	
	double time=0.;
	//-------------------проинициализируем начальные условия------------------------------------------
	for(int i = 0; i < myStep; i++)
	{
		U[i] = (double*)calloc(stepAmountX, sizeof(double));//need par
		Fin[i] = (double*)calloc(stepAmountX, sizeof(double));
	}
	double constUS = pow(s,alpha-1)/gamma(alpha);
	double constFive = (-A*pow(dt, alpha))/(2*pow(dx, 3));
	#pragma omp parallel for
		for(int i=0;i<stepAmountX;i++)
		{
			func[i]=fCos(i*dx);//начальное условие
			Us[i]=func[i]*pow(s,alpha-1)/gamma(alpha);
			fiveUp[i] = constFive;
			fiveDown[i] = -constFive;
			U[0][i]=0;	//начальные слои
			U[1][i]=0;
		}
	Coeff[0]=pow(-1.,0)*1.;
	Coeff[1]=pow(-1.,1)*alpha;
	for(int i=2;i<myStep;i++)
	{
		Coeff[i]=pow(-1.,i)*Coeff[i-1]*(alpha-i+1)/(i);
		U[i][0]=g(s+n*dt)-fCos(0)*pow(s,alpha-1)/gamma(alpha);//левый край
	}
	double constFunc = 2*dx;
	double constFunc2 = dx*dx;
	#pragma omp parallel for
		for(int i=1;i<stepAmountX-1;i++)
		{
			funcDer[i]=(func[i+1]-func[i-1])/constFunc;
			funcDer2[i]=(func[i+1]-2*func[i]+func[i-1])/constFunc2;
			funcDer3[i]=(funcDer2[i]-funcDer2[i-1])/(dx);
		}
	funcDer[0]=funcDer2[0]=funcDer3[0]=func[0];
	funcDer[stepAmountX-1]=funcDer2[stepAmountX-1]=funcDer3[stepAmountX-1]=func[stepAmountX-1];
	//-----------------------------------------------------------------------------------------
	//Явная схема
	
	int i = 0;
	for (n = 1; n < dodo; n++)
	{
		#pragma omp parallel
		{
			#pragma omp for
			for (i = 2; i < stepAmountX - 1; i++)//need par
			{
				sum2=sin(M_PI*alpha)*pow(s,alpha)*func[i]/(M_PI*(s+n*dt)*pow(n*dt,alpha))-pow(s,alpha-1)*sin(M_PI*alpha)*func[i]/(pow(n*dt,alpha)*M_PI);
				sum3=U[n][i]+func[i]*pow(s,alpha-1)/gamma(alpha);
				sum4=sin(M_PI*alpha)/(M_PI*alpha)*funcDer[i]*pow(s,alpha)*hyperGeometric(alpha,alpha,alpha+1,s/(s+(n-1)*dt))/pow(s+(n-1)*dt,alpha)+funcDer[i]*pow(s,alpha-1)/gamma(alpha);
				sum5=1/(2*dx)*(fracInt(n,1-alpha,i+1)-fracInt(n,1-alpha,i-1));
				sum6=A*(U[n][i-2]-3*U[n][i-1]+3*U[n][i]-U[n][i+1])/(dx*dx*dx)+A*pow(s,alpha-1)*funcDer3[i]/gamma(alpha);
				U[n+1][i]=-sum1+pow(dt,alpha)*(sum2-sum3*(sum4+sum5)+sum6);
			}
		}
		if(flag==true){
			U[n][0]=U[n][1]+Us[0]-Us[1];
		}else{
			U[n][0]=U[n][1]+Us[0]+Us[1];
		}
		U[n][stepAmountX-2]=U[n][stepAmountX-3]+Us[stepAmountX-3]-Us[stepAmountX-2];
		U[n][stepAmountX-1]=U[n][stepAmountX-2]+Us[stepAmountX-2]-Us[stepAmountX-1];
	}

	//----------------------------------------------------------------------------------------
	//Неявная схема
	double const1 = (-A*pow(dt,alpha)) / pow(dx,2.);
	double const2 = pow(dt,alpha)/2.;
	double const3 = (1./(2.*dx))*(pow(dt,1-alpha)/gamma(3-alpha));
	double const6 = ((A*pow(s,alpha-1)*pow(dt,alpha))/(gamma(alpha)));
	for(int n = dodo + 1; n < myStep-1; n++)
	{
		C1[0]=0;
		B1[0]=1;
		D1[0]=s+n*dt;
		double const4 = (M_PI*alpha*pow(s+(n-1)*dt,alpha));
		double const5 = (M_PI*pow(n,alpha)*(s+n*dt));
		double hyperConst1 = hyperGeometric(alpha,alpha,alpha+1,s/(s+(n-1)*dt));
		double hyperConst2 = hyperGeometric(alpha,alpha,alpha+1,s/(s+n*dt));
		double const7 = pow(s,alpha-1)/gamma(alpha);
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 2; i < stepAmountX - 2; i++)
		{
			sum1 = 0;
			for(int j = 1; j <= n; j++)
			{
				sum1 += Coeff[j] * U[n-j+1][i];
			}
			//i-1
			A1[i] = const1  -  const2*(U[n][i]+func[i]*const7)*const3;
			//i+1
			C1[i-1]= -A1[i] + 2*const1;
			//i
			B1[i-1]=1 + (const2)*((sin(M_PI*alpha)*pow(s,alpha)*funcDer[i]*hyperConst1)/const4+(funcDer[i]*const7) + (1./(2*dx))*(fracInt(n,1-alpha,i+1)-fracInt(n,1-alpha,i-1))) - 2 * const1;
			sumGroup=0;
			sumGroup2=0;
			//слагаемые - остатки от n+1 fracInt'а
			for (int j=1;j<=n-1;j++)
			{
				sumGroup+=U[j][i+1]*(pow(n+1-j,1-alpha)-pow(n-j,1-alpha))/(1-alpha)  +  (U[j+1][i+1]-U[j][i+1])*(pow(n+1-j,2-alpha)-pow(n-j,2-alpha))/((1-alpha)*(2-alpha))-pow(n-j,1-alpha)/(1-alpha);
			}
			for (int j=1;j<=n-1;j++)
			{
				sumGroup2+=U[j][i-1]*(pow(n+1-j,1-alpha)-pow(n-j,1-alpha))/(1-alpha)  +  (U[j+1][i-1]-U[j][i-1])*(pow(n+1-j,2-alpha)-pow(n-j,2-alpha))/((1-alpha)*(2-alpha))-pow(n-j,1-alpha)/(1-alpha);
			}
			D1[i-1]=-(sum1  -  ((sin(M_PI*alpha)*pow(s,alpha)*func[i])/const5) + ((pow(s,alpha-1)*func[i]*sin(M_PI*alpha))/(M_PI*pow(n,alpha)))  +  ((pow(dt,alpha)*func[i])*const7)*((sin(M_PI*alpha)*pow(s,alpha)*funcDer[i]*hyperConst1)/const4+(funcDer[i]*const7)+((1.)/(2*dx))*(fracInt(n,1-alpha,i+1)-fracInt(n,1-alpha,i-1)))  +  ((pow(dt,alpha))/2.)*(U[n][i]+func[i]*const7)*(sin(M_PI*alpha)*pow(s,alpha)*funcDer[i]*hyperConst2/(const5*alpha)+(funcDer[i]*const7)+(1./(2*dx))*(sumGroup-sumGroup2))  -  funcDer2[i]*const6);
		}
		//---(ищем A10, C9, B9, B10, D9, D10, при 11 узлах) формула универсальна для всех не найденных A B C D--------
		for(int i=stepAmountX-2;i<stepAmountX-1;i++)
		{
			for(int j=1;j<=n;j++)
			{
				sum1+=Coeff[j]*U[n-j][i];
			}
			//i-1
			A1[i]=-1;
			//i+1
			C1[i-1]=(-A*pow(dt,alpha))/(pow(dx,2)) + (pow(dt,alpha))/2.*(U[n][i]+func[i]*pow(s,alpha-1)/gamma(alpha))*((1)/(2.*dx))*(pow(dt,1-alpha)/gamma(3-alpha));
			//i
			B1[i-1]=1 + ((pow(dt,alpha))/2.)*((sin(M_PI*alpha)*pow(s,alpha)*funcDer[i]*hyperGeometric(alpha,alpha,alpha+1,s/(s+(n-1)*dt)))/(M_PI*alpha*pow(s+(n-1)*dt,alpha))+(funcDer[i]*pow(s,alpha-1)/(gamma(alpha))) + (1./(2*dx))*(fracInt(n,1-alpha,i+1)-fracInt(n,1-alpha,i-1))) + ((2*A*pow(dt,alpha))/pow(dx,2.));
			//последняя B (B10)
			B1[i]=1;
			sumGroup=0;
			sumGroup2=0;
			for (int j=1;j<=n-1;j++)
			{
				sumGroup+=U[j][i+1]*(pow(n+1-j,1-alpha)-pow(n-j,1-alpha))/(1-alpha)  +  (U[j+1][i+1]-U[j][i+1])*(pow(n+1-j,2-alpha)-pow(n-j,2-alpha))/((1-alpha)*(2-alpha))-pow(n-j,1-alpha)/(1-alpha);
			}
			for (int j=1;j<=n-1;j++)
			{
				sumGroup2+=U[j][i-1]*(pow(n+1-j,1-alpha)-pow(n-j,1-alpha))/(1-alpha)  +  (U[j+1][i-1]-U[j][i-1])*(pow(n+1-j,2-alpha)-pow(n-j,2-alpha))/((1-alpha)*(2-alpha))-pow(n-j,1-alpha)/(1-alpha);
			}
			
			D1[i-1]=-(sum1  -  ((sin(M_PI*alpha)*pow(s,alpha)*func[i])/(M_PI*pow(n,alpha)*(s+n*dt))) + ((pow(s,alpha-1)*func[i]*sin(M_PI*alpha))/(M_PI*pow(n,alpha)))  +  ((pow(dt,alpha)*func[i]*pow(s,alpha-1))/(2*gamma(alpha)))*((sin(M_PI*alpha)*pow(s,alpha)*funcDer[i]*hyperGeometric(alpha,alpha,alpha+1,s/(s+(n-1)*dt)))/(M_PI*alpha*pow(s+(n-1)*dt,alpha))+(funcDer[i]*pow(s,alpha-1)/gamma(alpha))+((1.)/(2*dx))*(fracInt(n,1-alpha,i+1)-fracInt(n,1-alpha,i-1)))  +  ((pow(dt,alpha))/2.)*(U[n][i]+func[i]*pow(s,alpha-1)/gamma(alpha))*(sin(M_PI*alpha)*pow(s,alpha)*funcDer[i]*hyperGeometric(alpha,alpha,alpha+1,s/(s+n*dt))/(M_PI*alpha*pow(s+n*dt,alpha))+(funcDer[i]*pow(s,alpha-1)/gamma(alpha))+(1./(2*dx))*(sumGroup-sumGroup2))  -  ((A*pow(s,alpha-1)*funcDer2[i]*pow(dt,alpha))/(gamma(alpha))));
			sumGroup=0;
			D1[i]=0;
		}
	}
		//----------------------------------------------------------------------------

		Shuttle(n);

		if(flag==true){
			U[n][0]=U[n][1]+Us[0]-Us[1];
		}else{
			U[n][0]=U[n][1]+Us[0]+Us[1];
		}
		U[n][stepAmountX-2]=U[n][stepAmountX-3]+Us[stepAmountX-3]-Us[stepAmountX-2];
		U[n][stepAmountX-1]=U[n][stepAmountX-2]+Us[stepAmountX-2]-Us[stepAmountX-1];
	}
	
	int m=myStep-1;
	if(flag==true){
		U[m][0]=U[m][1]+Us[0]-Us[1];
	}else{
		U[m][0]=U[m][1]+Us[0]+Us[1];
	}
	U[m][stepAmountX-2]=U[m][stepAmountX-3]+Us[stepAmountX-3]-Us[stepAmountX-2];
	U[m][stepAmountX-1]=U[m][stepAmountX-2]+Us[stepAmountX-2]-Us[stepAmountX-1];

	for(int l=0;l<myStep;l++)
	{
		#pragma omp parallel for
		for(int p=0;p<stepAmountX;p++)
		{
			Fin[l][p]=U[l][p]+func[p]*constUS;
		}
	}
	end=omp_get_wtime()-start;
	cout<<"time= "<<end<<" seconds\n";
	toFile(myStep,dx);
	system("testOMP.plt");
	system("plot.plt");
	return 0;
}