#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


double functionx(double rex, double vex, double rey, double vey, double t);
double (*funcx)(double rex, double vex, double rey, double vey, double t);
double functiony(double rex, double vex, double rey, double vey, double t);
double (*funcy)(double rex, double vex, double rey, double vey, double t);

//double functionS(double x, double v, double t);
//double (*funcS)(double x, double v, double t);

double rkn4(double (*funcx)(double rex, double vex, double rey, double vey, double t),double (*funcy)(double rex, double vex, double rey, double vey, double t) ,double rex[], double rex0, double vex[],double vex0,double rey[], double rey0, double vey[],double vey0,  double t, double h, long steps);

double functionx(double rex, double vex, double rey, double vey, double t){
	double me, ms,d, G, r0, rex0, rsx0, a, e;
	me     = 5.9736e24;          //[kg]
	ms     = 1.9891e30 ;         //[kg]
//	a      = 149598261e3 ;       //[m]
//	e      = 0.0167112303531389 ;
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
//	r0 = a*(1+e);
	//rex0 = r0;
//	rsx0 = -(me/ms)*r0;
	d = sqrt(pow(rex,2)+pow(rey,2));
	return -G/pow(d,3) * ms *rex;
        
}

double functiony(double rex, double vex, double rey, double vey, double t){
	double d, ms, G;
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
	ms     = 1.9891e30 ;         //[kg]	
	d = sqrt(pow(rex,2)+pow(rey,2));
	return -G/pow(d,3) * ms *rey;

}

int main() {

	//define constants to be used
        double me, ms, a, e, G, Tday, I1, I3, thet0, pi;
	pi = 3.14159;
	me     = 5.9736e24;          //[kg]
	ms     = 1.9891e30 ;         //[kg]
	a      = 149598261e3 ;       //[m]
	e      = 0.0167112303531389 ;
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
	Tday   = 86164.1    ;        //[s] -> SIDEREAL DAY
	I1     = 8.008e37   ;        //[kg m^2]
	I3     = 8.034e37 ;	    //[kg m^2]
	thet0 = 23.45*pi/180.;      //[radians]

	//initial conditions
	//vec(Re(0)) = {r0, 0, 0} ; vec(ve(0)) = {0,v0,0}
	double r0, KE, v0, Tyear,  ps0, ph0, dotTh0, dotPh0, dotPs0;
	double rex0, rey0, rez0, vex0, vey0, vez0;
	//earth radius (initially at aphelion)
	r0 = a*(1+e);
	rex0 = r0;
	rey0 = 0;
	rez0 = 0;
	//kinetic energy
	KE = (G*ms*me)/a * (1/(1+e) - 0.5);
	//earth velocity
	v0 = sqrt((2*KE)/(me*(1+me/ms)));
	vex0 = 0;
	vey0 = v0;
	vez0 = 0;
	//euler angles
	ph0 = 0	;	//phi_0
	ps0 = 0	;	//psi_0
	dotPh0 = 0;	//dot(phi_0)
	dotTh0 = 0;	//dot(theta_0)
	dotPs0 = 2*pi/Tday; //dot(psi_0)
	//initial conditions for sun
/*	double rsx0, rsy0, rsz0, vsx0, vsy0, vsz0;
	rsx0 = -(me/ms)*r0;
	rsy0 = 0;
	rsz0 = 0;
	vsx0 = 0;
	vsy0 = -(me/ms)*v0;
	vsz0 = 0;
	//constant based on init. conditions
	double c0 = dotPs0 + dotPh0*cos(thet0);*/
	//dist btwn earth and sun 
	double d;
	//d = sqrt(pow(rex,2)+pow(rey,2));





	long steps = 10000000; 				//Number of steps.
	double t=0,h=1000; 			//time, step size.
	double *rex = malloc(steps * sizeof(double)); 	//Position array, x, earth
        double *rsx = malloc(steps * sizeof(double)); 	//position array, x, sun
	double *vex = malloc(steps * sizeof(double)); 	//velocity array, x, earth
        double *vsx = malloc(steps * sizeof(double)); 	//Velocity array, x sun

	double *rey = malloc(steps * sizeof(double)); 	//Position array, y, earth
        double *rsy = malloc(steps * sizeof(double)); 	//position array, y, sun
	double *vey = malloc(steps * sizeof(double)); 	//velocity array, y, earth
        double *vsy = malloc(steps * sizeof(double)); 	//Velocity array, y, sun
	
	/*double *rez = malloc(steps * sizeof(double)); 	//Position array, z, earth
        double *rsz = malloc(steps * sizeof(double)); 	//position array, z, sun
	double *vez = malloc(steps * sizeof(double)); 	//velocity array, z, earth
        double *vsz = malloc(steps * sizeof(double)); 	//Velocity array, z sun*/



// Set pointer to point at the function to integrate. 
	funcx = functionx;
	funcy = functiony;

// Do integration. 
	rkn4(funcx,funcy, rex,rex0,vex,vex0,rey, rey0, vey, vey0, t,h,steps);
	//rkn4(func,rey,rey0,vey,vey0,rsy, rsy0, vsy, vsy0, t,h,steps);
	//rkn4(func,rez,rez0,vez,vez0,rsz, rsz0, vsz, vsz0, t,h,steps);

	
// Print results to STDOUT 
	long int i;
        for ( i=0; i<steps; ++i){
		t += h;
		if (i%100 ==0){printf(" %f %e %e\n",t,rex[i],rey[i]);}
        }



	free(rex);
	free(rsx);
	free(vex);
	free(vsx);
	free(rey);
	free(rsy);
	free(vey);
	free(vsy);
	/*free(rez);
	free(rsz);
	free(vez);
	free(vsz);*/

	
	return 0;
}


/*	//define derivative functions
	//dot(psi)
	double dotPsi(double c0, double dotPhi, double theta)
{	
	return c0 - dotPhi*cos(theta);
}


	//ddot(phi)
	double ddotPhi( double c0, double dotPsi, double dotThet, double dotPhi, double phi, double theta, double dx, double dy, double d3)
{
	return 1/sin(theta) * (-2*dotPhi*dotThet + I1/I3*c0*dotThet + ((3*G*ms)/d^5)*((I1-I3)/I1)*d3*(dx*cos(phi)+dy*sin(phi)));
}
	
	//ddot(theta)
	double ddotTheta()
{
}
*/

	//runge kutta nystrom 4

	double rkn4(double (*funcx)(double rex, double vex, double rey, double vey, double t),double (*funcy)(double rex, double vex, double rey, double vey, double t), double rex[], double rex0, double vex[],double vex0,double rey[], double rey0, double vey[],double vey0,  double t, double h, long steps){
	double k1x,k2x,k3x,k4x, k1y, k2y, k3y, k4y, me, ms;
	me     = 5.9736e24;          //[kg]
	ms     = 1.9891e30 ;         //[kg]

	
	rex[0] = rex0; //Initial position
	vex[0] = vex0; //Initial velocity
	rey[0] = rey0; //Initial position
	vey[0] = vey0; //Initial velocity

	
	long i;
	for ( i=1; i<steps; ++i){
		k1x = funcx(rex[i-1],vex[i-1],rey[i-1], vey[i-1],t);
		k1y = funcy(rex[i-1],vex[i-1],rey[i-1], vey[i-1],t);		

		k2x = funcx(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k1x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k1y/2.,t+h/2.);
		k2y = funcy(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k1x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k1y/2.,t+h/2.);

                k3x = funcx(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k2x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k2y/2.,t+h/2.);
                k3y = funcy(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k2x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k2y/2.,t+h/2.);

		k4x = funcx(rex[i-1]+h*vex[i-1]+h*h*k3x/2.,vex[i-1]+h*k3x, rey[i-1]+h*vey[i-1]+h*h*k3y/2.,vey[i-1]+h*k3y,t+h);
		k4y = funcy(rex[i-1]+h*vex[i-1]+h*h*k3x/2.,vex[i-1]+h*k3x, rey[i-1]+h*vey[i-1]+h*h*k3y/2.,vey[i-1]+h*k3y,t+h);

		


		vex[i] = vex[i-1] + h*(k1x + 2.*k2x + 2.*k3x + k4x)/6.;
		vey[i] = vey[i-1] + h*(k1y + 2.*k2y + 2.*k3y + k4y)/6.;
		rex[i] = rex[i-1] + h*vex[i-1] + h*h*(k1x+k2x+k3x)/6.;	
		rey[i] = rey[i-1] + h*vey[i-1] + h*h*(k1y+k2y+k3y)/6.;
		t+=h;
		//vex[i] = vex[i-1] + h*(k1 + 2.*k2 + 2.*k3 + k4)/6.;
		//vsx[i] = -me/ms * vex[i];
	//	rex[i] = rex[i-1] + h*vex[i-1] + h*h*(k1+k2+k3)/6.;
		//rsx[i] = -me/ms * rex[i];
		//t+=h;
	}
	
	return 0;
}

/*	//Runge Kutta 4
	double rk4(double(*f)(double,double), double tau, double t, double yn)
{
	double tn, tau, k1, k2, k3, k4, yn;
	//len = timeTot / tau
	//for(int n = 0; n < len; n++){

	//tau = 
	k1 = tau * f(tn, yn);
	k2 = tau * f(tn + tau/2, yn + 0.5*k1);
	k3 = tau * f(tn + tau/2, yn + 0.5*k2);
	k4 = tau * f(tn + tau, yn + k3);

	return yn + (1./6.)*(k1 + 2*k2 + 2*k3 + k4);
	//tnplus = tn + tau
}
*/

