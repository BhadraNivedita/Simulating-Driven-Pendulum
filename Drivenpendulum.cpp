		//Simulating Pendulum equation


		#include <stdio.h>

		#include <iostream>

		#include <math.h>

		#include <fstream>

		using namespace std;


		class pendulum

		{

		private:

		double Q, A_roof, omega_0, omega_roof,g; //

		double y[2];//for the initial-values of phi and v

		int n;// how many steps

		double delta_t,delta_t_roof;


		public:

		void derivatives(double,double*,double*);

		void initialise();

		void rk4_step(double,double*,double*,double); // we need it in function rk4() and asc()

		void rk4(); //runge-kutta-fourth-order

		};


		void pendulum::derivatives(double t, double* in, double* out)

		{

		out[0]=in[1];//out[0] = (phi)' = v

		if(Q)

		out[1]=-in[1]/((double)Q)-sin(in[0])+A_roof*cos(omega_roof *t); //out[1] = (phi)''

		else
		
		out[1]=-sin(in[0])+A_roof*cos(omega_roof*t); //out[1] = (phi)''

		}

		//Initialise

		void pendulum::initialise()

		{

		double m,l,omega,A,viscosity,phi_0,v_0,t_end;

		m=1;	l=1;	omega=1;

		A=1;	viscosity=.01;	y[0]=.01;

		y[1]=.001;	n=10000;
	
		t_end=10;	g=9.81;	

		t_end *= acos(-1.);

		omega_0=sqrt(g/((double)l));


		if (viscosity) Q= m*g/((double)omega_0*viscosity);

		else Q=0;

		A_roof=A/((double)m*g);
	
		omega_roof=omega/((double)omega_0);

		delta_t_roof=omega_0*t_end/((double)n); //delta_t without dimension

		delta_t=t_end/((double)n);
		
		}


		void pendulum::rk4_step(double t,double *yin,double *yout,double delta_t)

		{

		double k1[2],k2[2],k3[2],k4[2],y_k[2];

		derivatives(t,yin,yout);

		k1[1]=yout[1]*delta_t;

		k1[0]=yout[0]*delta_t;

		y_k[0]=yin[0]+k1[0]*0.5;

		y_k[1]=yin[1]+k1[1]*0.5;


		derivatives(t+delta_t*0.5,y_k,yout);

		k2[1]=yout[1]*delta_t;

		k2[0]=yout[0]*delta_t;

		y_k[0]=yin[0]+k2[0]*0.5;

		y_k[1]=yin[1]+k2[1]*0.5;


		derivatives(t+delta_t*0.5,y_k,yout);

		k3[1]=yout[1]*delta_t;

		k3[0]=yout[0]*delta_t;

		y_k[0]=yin[0]+k3[0];

		y_k[1]=yin[1]+k3[1];


		derivatives(t+delta_t,y_k,yout);

		k4[1]=yout[1]*delta_t;

		k4[0]=yout[0]*delta_t;

		yout[0]=yin[0]+1.0/6.0*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);

		yout[1]=yin[1]+1.0/6.0*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
	
		}


		void pendulum::rk4()

		{

		int i;

		double t_h;

		double yout[2],y_h[2]; //k1[2],k2[2],k3[2],k4[2],y_k[2];

		t_h=0;

		y_h[0]=y[0]; //phi

		y_h[1]=y[1]; //v


		ofstream fout("rk4.out");

		fout.setf(ios::scientific);

		fout.precision(20);


		for(i=1; i<=n; i++){
		
		rk4_step(t_h,y_h,yout,delta_t_roof);

		fout<<i*delta_t<<"\t\t"<<yout[0]<<"\t\t"<<yout[1]<<"\n";

		t_h+=delta_t_roof;

		y_h[0]=yout[0];
		
		y_h[1]=yout[1];

		}

		fout.close();

		}
//////////////////////////////////////////////////////////////////////////////////////////

		int main()

		{

		pendulum testcase;

		testcase.initialise();

		testcase.rk4();

		return 0;

		}
