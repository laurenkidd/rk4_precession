#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


double functionx(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t);
double (*funcx)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t);
double functiony(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t);
double (*funcy)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t);
double functiony(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t);
double (*funcy)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t);
double functionz(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t);
double (*funcz)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t);

double functiontheta(double rex, double rey, double rez, double rsx,double rsy,  double rsz, double theta, double vtheta, double phi, double vphi);
double (*funct)(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi);
double functionphi(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi);
double (*funcph)(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi);


//double functionS(double x, double v, double t);
//double (*funcS)(double x, double v, double t);

double rkn4(double (*funcx)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t),
double (*funcy)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t),
 double (*funcz)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t),
 double (*funct)(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi),
double (*funcph)(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi),
double rex[], double rex0, double vex[],double vex0,double rey[], double rey0, double vey[],double vey0,  double rez[], double rez0, double vez[], double vez0, double rsx[], double rsx0,
double vsx[],  double rsy[], double rsy0, double vsy[], double rsz[], double rsz0, double vsz[], double phi[], double phi0, double vphi[], double vphi0, double theta[], double theta0, double vtheta[], double vtheta0, double psi[], double t, double h, long steps);

double functionx(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t){
	//x position in CM-centered frame
	double me, ms,d, G, I3, I1;
	//double pi = 3.14159;
	me     = 5.9736e24;          //[kg]
	ms     = 1.9891e30 ;         //[kg]
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
	I1     = 8.008e37   ;        //[kg m^2]
	I3     = 8.034e37 ;	    //[kg m^2]
	d = fabs(sqrt(pow(rsx,2)+pow(rsy,2)+pow(rsz,2)) - sqrt(pow(rex,2)+pow(rey,2)+pow(rez,2)));
	double dx, dy, dz, d3;
	dx = fabs(rsx - rex);
	dy = fabs(rsy - rey);
	dz = fabs(rsz - rez);
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	printf("%f  %f %f %f %f\n", phi, theta, dy, dz, d3);
	double fx;
	fx = (-(G*me*ms)/(d*d*d))*dx + ((3*G*ms)/pow(d,5))*(I1-I3)*(dx/2. - (2.5*d3*d3*dx)/(d*d) + d3*(sin(phi)*sin(theta)));
	return -fx/me;
        
}

double functiony(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t){
	//y position in CM-centered frame
	double d, ms, G, I3, I1, me;
	//double pi = 3.14159;
	me     = 5.9736e24;          //[kg]
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
	ms     = 1.9891e30 ;         //[kg]
	I1     = 8.008e37   ;        //[kg m^2]
	I3     = 8.034e37 ;	    //[kg m^2]	
	d = fabs(sqrt(pow(rsx,2)+pow(rsy,2)+pow(rsz,2)) - sqrt(pow(rex,2)+pow(rey,2)+pow(rez,2)));
	double dx, dy, dz, d3;
	dx = fabs(rsx - rex);
	dy = fabs(rsy - rey);
	dz = fabs(rsz - rez);
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	double fy;
	fy = (-(G*me*ms)/(d*d*d))*dy + ((3*G*ms)/pow(d,5))*(I1-I3)*(dy/2. - (2.5*d3*d3*dy)/(d*d) + d3*(-cos(phi)*sin(theta)));
	return -fy/me;
}

double functionz(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t){
	//y position in CM-centered frame
	double d, ms, G, I3, I1, me;
	//double pi = 3.14159;
	me     = 5.9736e24;          //[kg]
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
	ms     = 1.9891e30 ;         //[kg]
	I1     = 8.008e37   ;        //[kg m^2]
	I3     = 8.034e37 ;	    //[kg m^2]	
	d = fabs(sqrt(pow(rsx,2)+pow(rsy,2)+pow(rsz,2)) - sqrt(pow(rex,2)+pow(rey,2)+pow(rez,2)));
	double dx, dy, dz, d3;
	dx = fabs(rsx - rex);
	dy = fabs(rsy - rey);
	dz = fabs(rsz - rez);
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	double fz;
	fz = (-(G*me*ms)/(d*d*d))*dz + ((3.0*G*ms)/pow(d,5))*(I1-I3)*(dz/2. - (2.5*d3*d3*dz)/(d*d) + d3*(cos(theta)));
	return -fz/me;
}

double functionphi(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi){
	double I1, I3, d, G, ms, Tday;
	double pi = 3.14159;
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
	ms     = 1.9891e30 ;         //[kg]
	I1     = 8.008e37   ;        //[kg m^2]
	I3     = 8.034e37 ;	    //[kg m^2]	
	Tday   = 86164.1    ;        //[s] -> SIDEREAL DAY
	d = fabs(sqrt(pow(rsx,2)+pow(rsy,2)+pow(rsz,2)) - sqrt(pow(rex,2)+pow(rey,2)+pow(rez,2)));
	double dx, dy, dz, d3, c0;
	dx = fabs(rsx - rex);
	dy = fabs(rsy - rey);
	dz = fabs(rsz - rez);
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	c0 =  2*pi/Tday;
	//(1.0/sin(theta)*(-2.0*vphi*vtheta*cos(theta)+(I3/I1)*c0*vtheta+((3*G*ms)/pow(d,5))*((I1-I3)/I1)*d3*(dx*cos(phi)+dy*sin(phi)))
	return (1.0/sin(theta))*(((3*G*ms)/pow(d,5))*((I1-I3)/I1)*d3*(dx*cos(phi)+dy*sin(phi)));

}

double functiontheta(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi){
	double I1, I3, d, G, ms, Tday;
	double pi = 3.14159;
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
	ms     = 1.9891e30 ;         //[kg]
	I1     = 8.008e37   ;        //[kg m^2]
	I3     = 8.034e37 ;	    //[kg m^2]	
	Tday   = 86164.1    ;        //[s] -> SIDEREAL DAY
	d = fabs(sqrt(pow(rsx,2)+pow(rsy,2)+pow(rsz,2)) - sqrt(pow(rex,2)+pow(rey,2)+pow(rez,2)));
	double dx, dy, dz, d3, c0;
	dx = fabs(rsx - rex);
	dy = fabs(rsy - rey);
	dz = fabs(rsz - rez);
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	c0 =  2*pi/Tday;
	return vphi*vphi*sin(theta)*cos(theta)-(I3/I1)*c0*sin(theta)+(3*G*ms)/pow(d,5)*((I1-I3)/I1)*d3*((dx*sin(phi)-dy*cos(phi))*cos(theta)-dz*sin(theta));
}


int main() {

	//define constants to be used
        double me, ms, a, e, G, Tday, I1, I3, theta0, pi;
	pi = 3.14159;
	me     = 5.9736e24;          //[kg]
	ms     = 1.9891e30 ;         //[kg]
	a      = 149598261e3 ;       //[m]
	e      = 0.0167112303531389 ;
	G      = 6.67428e-11  ;      //[N m^2 kg^-2]
	Tday   = 86164.1    ;        //[s] -> SIDEREAL DAY
	I1     = 8.008e37   ;        //[kg m^2]
	I3     = 8.034e37 ;	    //[kg m^2]
	theta0 = 23.45*pi/180.;      //[radians]

	//initial conditions
	//vec(Re(0)) = {r0, 0, 0} ; vec(ve(0)) = {0,v0,0}
	double r0, KE, v0,  psi0, phi0, vtheta0, vphi0, vpsi0;
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
	phi0 = 0	;	//phi_0
	psi0 = 0	;	//psi_0
	vphi0 = 0;	//dot(phi_0)
	vtheta0 = 0;	//dot(theta_0)
	vpsi0 = 2*pi/Tday; //dot(psi_0)
	//initial conditions for sun
	double rsx0, rsy0, rsz0, vsx0, vsy0, vsz0;
	rsx0 = -(me/ms)*rex0;
	rsy0 = 0;
	rsz0 = 0;
	vsx0 = 0;
	vsy0 = -(me/ms)*v0;
	vsz0 = 0;
	//constant based on init. conditions
	//double c0 = 2*pi/Tday;
	//dist btwn earth and sun 
	




	long steps = 100000; 		//00		//Number of steps.
	double t=0,h=10000; 			//time, step size.
	//define arrays to be used
	double *rex = malloc(steps * sizeof(double)); 	//Position array, x, earth
        double *rsx = malloc(steps * sizeof(double)); 	//position array, x, sun
	double *vex = malloc(steps * sizeof(double)); 	//velocity array, x, earth
        double *vsx = malloc(steps * sizeof(double)); 	//Velocity array, x sun

	double *rey = malloc(steps * sizeof(double)); 	//Position array, y, earth
        double *rsy = malloc(steps * sizeof(double)); 	//position array, y, sun
	double *vey = malloc(steps * sizeof(double)); 	//velocity array, y, earth
        double *vsy = malloc(steps * sizeof(double)); 	//Velocity array, y, sun
	
	double *rez = malloc(steps * sizeof(double)); 	//Position array, z, earth
        double *rsz = malloc(steps * sizeof(double)); 	//position array, z, sun
	double *vez = malloc(steps * sizeof(double)); 	//velocity array, z, earth
        double *vsz = malloc(steps * sizeof(double)); 	//Velocity array, z sun
	
	double *theta = malloc(steps * sizeof(double)); 	
        double *phi = malloc(steps * sizeof(double)); 	
	double *psi = malloc(steps * sizeof(double)); 
        double *vtheta = malloc(steps * sizeof(double)); 
	double *vphi = malloc(steps * sizeof(double)); 	



// Set pointer to point at the function to integrate. 
	funcx = functionx;
	funcy = functiony;
	funcz = functionz;
	funct = functiontheta;
	funcph = functionphi;

// Do integration. 

/*double rex[], double rex0, double vex[],double vex0,double rey[], double rey0, double vey[],double vey0,  double rez[], double rez0, double vez[], double vez0, double rsx[], double rsx0,
double vsx[],  double rsy[], double rsy0, double vsy[], double rsz[], double rsz0, double vsz[], double phi[], double phi0, double vphi[], double vphi0, double theta[], double theta0, double vtheta[], double vtheta0, double psi[], double t, double h, long steps);*/

	rkn4(funcx,funcy,funcz,funct,funcph, rex,rex0,vex,vex0,rey, rey0, vey, vey0, rez, rez0, vez, vez0, rsx, rsx0, vsx, rsy, rsy0, vsy, rsz, rsz0, vsz, phi, phi0, vphi, vphi0, theta,
		 theta0, vtheta, vtheta0, psi, t,h,steps);

	
// Print results to STDOUT 
	/*
	long int i;
        for ( i=0; i<steps; ++i){
		t += h;
		if (i%100 ==0){printf(" %f %f %f\n",t,rex[i],rey[i]);}
        }*/



	free(rex);
	free(rsx);
	free(vex);
	free(vsx);
	free(rey);
	free(rsy);
	free(vey);
	free(vsy);
	free(rez);
	free(rsz);
	free(vez);
	free(vsz);
	free(theta);
	free(phi);
	free(psi);
	free(vtheta);
	free(vphi);

	
	return 0;
}



	//runge kutta nystrom 4

	double rkn4(double (*funcx)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t),
double (*funcy)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t),
 double (*funcz)(double rex, double vex, double rey, double vey, double rez, double vez, double rsx, double rsy, double rsz, double phi, double theta, double t),
 double (*funct)(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi),
double (*funcph)(double rex, double rey, double rez, double rsx, double rsy, double rsz, double theta, double vtheta, double phi, double vphi),
double rex[], double rex0, double vex[],double vex0,double rey[], double rey0, double vey[],double vey0,  double rez[], double rez0, double vez[], double vez0, double rsx[], double rsx0,
double vsx[],  double rsy[], double rsy0, double vsy[], double rsz[], double rsz0, double vsz[], double phi[], double phi0, double vphi[], double vphi0, double theta[], double theta0, double vtheta[], double vtheta0, double psi[], double t, double h, long steps){
	double k1x,k2x,k3x,k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z, me, ms;
	double k1t, k2t, k3t, k4t, k1ph, k2ph, k3ph, k4ph, k1sx, k1sy, k1sz, k2sx, k2sy, k2sz, k3sx, k3sy, k3sz, k4sx, k4sy, k4sz;
	me     = 5.9736e24;          //[kg]
	ms     = 1.9891e30 ;         //[kg]
	double pi = 3.14159;
	double Tday   = 86164.1    ;        //[s] -> SIDEREAL Day
	double c0 = 2*pi/Tday;

	
	rex[0] = rex0; //Initial position
	vex[0] = vex0; //Initial velocity
	rey[0] = rey0; //Initial position
	vey[0] = vey0; //Initial velocity
	rez[0] = rez0;
	vez[0] = vez0;
	rsx[0] = rsx0;
	rsy[0] = rsy0;
	rsz[0] = rsz0;
	theta[0] = theta0;
	phi[0] = phi0;
	vtheta[0] = vtheta0;
	vphi[0] = vphi0;
	vsx[0] = -vex0*me/ms;
	vsy[0] = -vey0*me/ms;
	vsz[0] = -vez0*me/ms;


	
	long i;
	for ( i=1; i<steps; ++i){
		//double vsx, vsy, vsz;
		
		//1//

		k1x = funcx(rex[i-1],vex[i-1],rey[i-1], vey[i-1], rez[i-1], vez[i-1], rsx[i-1], rsy[i-1], rsz[i-1], phi[i-1], theta[i-1], t);
		//printf("%f",k1x);
		k1sx = -k1x*me/ms;
		k1y = funcy(rex[i-1],vex[i-1],rey[i-1], vey[i-1], rez[i-1], vez[i-1], rsx[i-1], rsy[i-1], rsz[i-1], phi[i-1], theta[i-1], t);
		k1sy = -k1y*me/ms;
		k1z = funcx(rex[i-1],vex[i-1],rey[i-1], vey[i-1], rez[i-1], vez[i-1], rsx[i-1], rsy[i-1], rsz[i-1], phi[i-1], theta[i-1], t);
		k1sz = -k1z*me/ms;
		k1t = funct(rex[i-1], rey[i-1], rez[i-1], rsx[i-1], rsy[i-1], rsz[i-1], theta[i-1], vtheta[i-1], phi[i-1], vphi[i-1] );
		k1ph = funcph(rex[i-1], rey[i-1], rez[i-1], rsx[i-1], rsy[i-1], rsz[i-1], theta[i-1], vtheta[i-1], phi[i-1], vphi[i-1]);		

		//2//

		k2x = funcx(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k1x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k1y/2.,rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k1z/2., 				rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,rsy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,theta[i-1]+
			h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k2sx = -k2x*me/ms;
		k2y = funcy(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k1x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k1y/2.,rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k1z/2., 				rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,rsy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,theta[i-1]+
			h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k2sy = -k2y*me/ms;
		k2z = funcz(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k1x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k1y/2.,rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k1z/2., 				rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,rsy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,theta[i-1]+
			h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k2sz = -k2z*me/ms;		
		k2t = funct(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8., rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8., rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8., rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8., 						rsy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8.,theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., vtheta[i-1]+h*k1t/2.,  phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8., 				vphi[i-1]+h*k1ph/2. );
		k2ph = funcph(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8., rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8., rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8., rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8., 						rsy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8.,theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., vtheta[i-1]+h*k1t/2.,  phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8., 				vphi[i-1]+h*k1ph/2. );

		//3//
	
                k3x = funcx(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k2x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k2y/2.,rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k2z/2., 				rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,rsy[i-1]+h*vsy[i-1]/2.+h*h*k2sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,
			theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k3sx = -k3x*me/ms;
                k3y = funcy(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k2x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k2y/2.,rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k2z/2., 				rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,rsy[i-1]+h*vsy[i-1]/2.+h*h*k2sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,
			theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k3sy = -k3y*me/ms;
		k3z = funcz(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k2x/2.,rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k2y/2.,rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k2z/2., 				rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,rsy[i-1]+h*vsy[i-1]/2.+h*h*k2sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,
			theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k3sz = -k3z*me/ms;
		k3t = funct(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8., rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8., rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8., rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8., 						rsy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8.,theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., vtheta[i-1]+h*k2t/2.,  
			phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8., vphi[i-1]+h*k2ph/2. );
		k3ph = funcph(rex[i-1]+h*vex[i-1]/2.+h*h*k1x/8., rey[i-1]+h*vey[i-1]/2.+h*h*k1y/8., rez[i-1]+h*vez[i-1]/2.+h*h*k1z/8., rsx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8., 						rsy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., rsz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8.,theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., vtheta[i-1]+h*k2t/2.,  
			phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8., vphi[i-1]+h*k2ph/2. );
	

		//4//
		k4x = funcx(rex[i-1]+h*vex[i-1]+h*h*k3x/2.,vex[i-1]+h*k3x, rey[i-1]+h*vey[i-1]+h*h*k3y/2.,vey[i-1]+h*k3y, rez[i-1]+h*vez[i-1]+h*h*k3z/2.,vez[i-1]+h*k3z,
			rsx[i-1]+h*vsx[i-1]+h*h*k3sx/2., rsy[i-1]+h*vsy[i-1]+h*h*k3y/2., rsz[i-1]+h*vsz[i-1]+h*h*k3z/2., phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., 
			theta[i-1]+h*vtheta[i-1]+h*h*k3t/2., t+h);
		k4sx = -k4x*me/ms;
		k4y = funcy(rex[i-1]+h*vex[i-1]+h*h*k3x/2.,vex[i-1]+h*k3x, rey[i-1]+h*vey[i-1]+h*h*k3y/2.,vey[i-1]+h*k3y, rez[i-1]+h*vez[i-1]+h*h*k3z/2.,vez[i-1]+h*k3z,
			rsx[i-1]+h*vsx[i-1]+h*h*k3sx/2., rsy[i-1]+h*vsy[i-1]+h*h*k3y/2., rsz[i-1]+h*vsz[i-1]+h*h*k3z/2., phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., 
			theta[i-1]+h*vtheta[i-1]+h*h*k3t/2., t+h);
		k4sy = -k4y*me/ms;
		k4z = funcz(rex[i-1]+h*vex[i-1]+h*h*k3x/2.,vex[i-1]+h*k3x, rey[i-1]+h*vey[i-1]+h*h*k3y/2.,vey[i-1]+h*k3y, rez[i-1]+h*vez[i-1]+h*h*k3z/2.,vez[i-1]+h*k3z,
			rsx[i-1]+h*vsx[i-1]+h*h*k3sx/2., rsy[i-1]+h*vsy[i-1]+h*h*k3y/2., rsz[i-1]+h*vsz[i-1]+h*h*k3z/2., phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., 
			theta[i-1]+h*vtheta[i-1]+h*h*k3t/2., t+h);
		k4sz = -k4z*me/ms;
		k4t = funct(rex[i-1]+h*vex[i-1]+h*h*k3x/2., rey[i-1]+h*vey[i-1]+h*h*k3y/2., rez[i-1]+h*vez[i-1]+h*h*k3z/2.,rsx[i-1]+h*vsx[i-1]+h*h*k3sx/2., rsy[i-1]+h*vsy[i-1]+h*h*k3y/2., 			rsz[i-1]+h*vsz[i-1]+h*h*k3z/2., theta[i-1]+h*vtheta[i-1]+h*h*k3t/2.,vtheta[i-1]+h*k3t, phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., vphi[i-1]+h*k3ph);
		k4ph = funcph(rex[i-1]+h*vex[i-1]+h*h*k3x/2., rey[i-1]+h*vey[i-1]+h*h*k3y/2., rez[i-1]+h*vez[i-1]+h*h*k3z/2.,rsx[i-1]+h*vsx[i-1]+h*h*k3sx/2., 
			rsy[i-1]+h*vsy[i-1]+h*h*k3y/2., rsz[i-1]+h*vsz[i-1]+h*h*k3z/2., theta[i-1]+h*vtheta[i-1]+h*h*k3t/2.,vtheta[i-1]+h*k3t, phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., 					vphi[i-1]+h*k3ph);


		//next velocity term
		vex[i] = vex[i-1] + h*(k1x + 2.*k2x + 2.*k3x + k4x)/6.;
		vsx[i] = vsx[i-1] + h*(k1sx + 2.*k2sx + 2.*k3sx + k4sx)/6.;

		vey[i] = vey[i-1] + h*(k1y + 2.*k2y + 2.*k3y + k4y)/6.;
		vsy[i] = vsy[i-1] + h*(k1sy + 2.*k2sy + 2.*k3sy + k4sy)/6.;

		vez[i] = vez[i-1] + h*(k1z + 2.*k2z + 2.*k3z + k4z)/6.;
		vsz[i] = vsz[i-1] + h*(k1sz + 2.*k2sz + 2.*k3sz + k4sz)/6.;
		
		//again testing const theta
		vtheta[i] = 0.;
		//vtheta[i] = vtheta[i-1] + h*(k1t + 2.*k2t + 2.*k3t + k4t)/6.;
		vphi[i] = vphi[i-1] + h*(k1ph + 2.*k2ph + 2.*k3ph + k4ph)/6.;

		//next position term
		rex[i] = rex[i-1] + h*vex[i-1] + h*h*(k1x+k2x+k3x)/6.;	
		rsx[i] = -rex[i]*me/ms;

		rey[i] = rey[i-1] + h*vey[i-1] + h*h*(k1y+k2y+k3y)/6.;
		rsy[i] = -rey[i]*me/ms;

		rez[i] = rez[i-1] + h*vez[i-1] + h*h*(k1z+k2z+k3z)/6.;
		rsz[i] = -rez[i]*me/ms;

		//theta[i] = theta[i-1] + h*vtheta[i-1] + h*h*(k1t+k2t+k3t)/6.;
		//test const theta, because theta was blowing up 
		theta[i] = 23.45*pi/180.;
		phi[i] = phi[i-1] + h*vphi[i-1] + h*h*(k1ph+k2ph+k3ph)/6.;
	
		//psi integration 
		double k1, k2, k3, k4;
		k1 = c0 + vphi[i-1]*cos(theta[i-1]); 
       		k2 = c0 + (vphi[i-1]+h*k1ph/2.)*cos(theta[i-1]+h*k1t/2.); 
        	k3 = c0 + (vphi[i-1]+h*k2ph/2.)*cos(theta[i-1]+h*k2t/2.); 
        	k4 = c0 + (vphi[i-1]+h*k3ph)*cos(theta[i-1]+h*k3t); 
        	psi[i]= (h/6.)*(k1+2*k2+2*k3+k4);


		t+=h;
		}	
		return 0;
}


