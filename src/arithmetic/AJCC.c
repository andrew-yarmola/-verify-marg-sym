#include "AJCC.h"

const AJCC operator*(const AJCC&x,const AJCC&y) {

	double xdist = size(x);
	double ydist = size(y);
	double ax = absUB(x.f), ay = absUB(y.f);
	AComplex r_f = x.f*y.f;
	AComplex r_z0 = x.f*y.z0+x.z0*y.f;
	AComplex r_z1 = x.f*y.z1+x.z1*y.f;
	AComplex r_w0 = x.f*y.w0+x.w0*y.f;
	AComplex r_w1 = x.f*y.w1+x.w1*y.f;
	double A = (xdist + x.e) * (ydist + y.e);
	double B = ax * y.e + ay * x.e;
	double C = (r_f.e + (r_z0.e + r_w0.e)) + (r_z1.e + r_w1.e);
	double r_error = (1 + 3 * EPS) * ((A + B) + C);
	return AJCC(r_f.z, r_z0.z, r_z1.z, r_w0.z, r_w1.z, r_error);

}

const AJCC operator/(const AJCC&x,const AJCC&y) {

	double xdist = size(x);
	double ydist = size(y);
	double ax = absUB(x.f), ay = absLB(y.f);
	double D = ay-(1+EPS)*(y.e+ydist);
	if(!(D > 0)) return AJCC(0,0,0,0,0,infinity());
	AComplex den = (y.f*y.f);
	AComplex r_f = x.f/y.f;
	AComplex r_z0 = (x.z0*y.f-x.f*y.z0)/den;
	AComplex r_z1 = (x.z1*y.f-x.f*y.z1)/den;
	AComplex r_w0 = (x.w0*y.f-x.f*y.w0)/den;
	AComplex r_w1 = (x.w1*y.f-x.f*y.w1)/den;
	double A = (ax + (xdist + x.e))/D;
	double B = (ax/ay + xdist/ay) + (ydist * ax)/(ay * ay);
	double C = (r_f.e + (r_z0.e + r_w0.e)) + (r_z1.e + r_w1.e);
	double r_error = (1 + 3 * EPS) * (((1 + 3 * EPS) * A-(1-3 * EPS) * B) + C);
	return AJCC(r_f.z, r_z0.z, r_z1.z, r_w0.z, r_w1.z, r_error);

}

const AJCC operator/(double x,const AJCC&y) {

	double ydist = size(y);
	double ax = fabs(x), ay = absLB(y.f);
	double D = ay-(1+EPS)*(y.e+ydist);
	if(!(D > 0))return AJCC(0,0,0,0,0,infinity());
	AComplex den = (y.f*y.f);
	AComplex r_f = x/y.f;
	AComplex r_z0 = (-x*y.z0)/den;
	AComplex r_z1 = (-x*y.z1)/den;
	AComplex r_w0 = (-x*y.w0)/den;
	AComplex r_w1 = (-x*y.w1)/den;
	double B = ax/ay + (ydist * ax)/(ay * ay);
	double C = (r_f.e + (r_z0.e + r_w0.e)) + (r_z1.e + r_w1.e);
	double r_error = (1 + 3 * EPS) * (((1 + 2 * EPS) * (ax/D)-(1-3 * EPS) * B) + C);
	return AJCC(r_f.z, r_z0.z, r_z1.z, r_w0.z, r_w1.z, r_error);

}

const AJCC sqrt(const AJCC&x) {

	double xdist = size(x);
	double ax = absUB(x.f);
	double D = ax-(1 + EPS) * (xdist + x.e);
	if(!(D > 0)) {
		return AJCC(0,0,0,0,0,
        (1 + 2 * EPS) * sqrt(ax + (xdist + x.e)));
	} else {
		AComplex r_f = sqrt(x.f);
		AComplex t = r_f+r_f;
		AComplex r_z0 = AComplex(x.z0.re,x.z0.im,0)/t;
		AComplex r_z1 = AComplex(x.z1.re,x.z1.im,0)/t;
		AComplex r_w0 = AComplex(x.w0.re,x.w0.im,0)/t;
		AComplex r_w1 = AComplex(x.w1.re,x.w1.im,0)/t;
		double C = (r_f.e + (r_z0.e + r_w0.e)) + (r_z1.e + r_w1.e);
		double r_error = (1 + 3 * EPS) * (
							 ((1 + EPS) * sqrt(ax) - (1 - 3 * EPS) * (xdist/(2 * sqrt(ax)) + sqrt(D))) + C);
		return AJCC(r_f.z,r_z0.z,r_z1.z,r_w0.z,r_w1.z,r_error);
	}

}
