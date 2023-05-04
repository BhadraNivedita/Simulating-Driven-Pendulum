//Simulating Pendulum equation

//Simulation thetadotdot+omega0sin(theta)=0

#include <stdio.h>

#include <iostream>

#include <math.h>

#include <fstream>

using namespace std;


	class pendelum

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


				void pendelum::derivatives(double t, double* in, double* out)
					{

						out[0]=in[1];//out[0] = (phi)' = v

						out[1]=-omega_0*sin(in[0]); //out[1] = (phi)''=-sin(x)

					}

					//Initialise

				void pendelum::initialise()
					
					{
	
					double m,l,omega,A,viscosity,phi_0,v_0,t_end;

					m=1;	l=1;	omega=1;	A=1;	viscosity=.01;		n=10000;	t_end=10;		
					
					g=9.81;

					y[0]=.01;

					y[1]=.001;

					t_end *= acos(-1.);
	
					omega_0=sqrt(g/((double)l));

					delta_t_roof=omega_0*t_end/((double)n); //delta_t without dimension
		
					delta_t=t_end/((double)n);
					
					}


				void pendelum::rk4_step(double t,double *yin,double *yout,double delta_t)
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

				void pendelum::rk4()

					{

					int i;

					double t_h;

					double yout[2],y_h[2]; //k1[2],k2[2],k3[2],k4[2],y_k[2];

					t_h=0;

					y_h[0]=y[0]; //phi

					y_h[1]=y[1]; //v

					ofstream fout("simplependulum.dat");

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

					int main()

						{

						pendelum testcase;

						testcase.initialise();
  
						testcase.rk4();

						return 0;
			
						}
