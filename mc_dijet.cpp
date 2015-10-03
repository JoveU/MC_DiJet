#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multimin.h>
#include "interp2d/interp2d.h"
#include "interp2d/interp2d_spline.h"
#include <iostream>
#include <fstream>
#include <string>
#include <random>

/*
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/factorials.hpp>
*/
using namespace std;
//using namespace blitz;

typedef mt19937 RNGType;
const size_t Y_size=21, Pt_size=178;
const size_t number_of_events=1000000;
const double sqrts=100.0;

unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;
    a=a-c;
    a=a^(c >> 13);
    b=b-c;
    b=b-a;
    b=b^(a << 8);
    c=c-a;
    c=c-b;
    c=c^(b >> 13);
    a=a-b;
    a=a-c;
    a=a^(c >> 12);
    b=b-c;
    b=b-a;
    b=b^(a << 16);
    c=c-a;
    c=c-b;
    c=c^(b >> 5);
    a=a-b;
    a=a-c;
    a=a^(c >> 3);
    b=b-c;
    b=b-a;
    b=b^(a << 10);
    c=c-a;
    c=c-b;
    c=c^(b >> 15);
    return c;
}




class TMD
{
    interp2d_spline* interp_xG;
    interp2d_spline* interp_xH;
    gsl_interp_accel *xa, *ya;
    const interp2d_type* T;

    double* Y_arr;
    double* Pt_arr;
    double* xG_arr;
    double* xH_arr;

    double sqrtS;
    double s;
    double prefactor;
    static const double ratio_Q_over_Pt = 2.0;
    static const double sum_charge2 = .666666;
    static const double alpha_em = .00729927;
    static const double alpha_s = .15;
    static const double S_perp = 1300.0; //mb
    static const double S_perp_JIMWLK =  2704.0; //mb

    void load_data(void);
    double get_Q(double Pt);
    double get_epsf2(double Q, double z);
    double get_Xsection(double Pt, double qt, double z, double phi);
    bool Out_of_Kinematic_Constrains(double Pt, double qt, double x);


public:
    static const double x0 = 1e-2;
    TMD(double sqrtSin);
    double get_x0(void)
    {
        return x0;
    }
    double get_xG_at_Y_qt(double Y, double qt);
    double get_xH_at_Y_qt(double Y, double qt);
    double operator() (double Pt, double qt, double z, double phi)
    {
        return get_Xsection(Pt,qt,z,phi);
    }

    double get_x(double Pt, double qt, double z);
};

TMD::TMD(double sqrtSin): sqrtS(sqrtSin)
{
    s = sqrtS*sqrtS;
    prefactor = alpha_em * alpha_s * sum_charge2;

    T=interp2d_bilinear;
    xa = gsl_interp_accel_alloc();
    ya = gsl_interp_accel_alloc();
    interp_xG = interp2d_spline_alloc(T, Pt_size, Y_size);
    interp_xH = interp2d_spline_alloc(T, Pt_size, Y_size);

    Y_arr = new double[Y_size];
    Pt_arr = new double[Pt_size];
    xG_arr = new double[Pt_size*Y_size];
    xH_arr = new double[Pt_size*Y_size];

    load_data();

    interp2d_spline_init(interp_xG, Pt_arr, Y_arr, xG_arr, Pt_size, Y_size);
    interp2d_spline_init(interp_xH, Pt_arr, Y_arr, xH_arr, Pt_size, Y_size);
}

double TMD::get_x(double Pt, double qt, double z)
{
    double x = ( qt*qt + Pt*Pt/(z*(1.0-z)) ) / s;
    return x;
}

double TMD::get_Q(double Pt)
{
    double Q = ratio_Q_over_Pt * Pt;
    return Q;
}

double TMD::get_epsf2(double Q, double z)
{
    double epsf2 = z*(1.0-z)*Q*Q;
    return epsf2;
}

bool TMD::Out_of_Kinematic_Constrains(double Pt, double qt, double x)
{
    if(x>x0) return true;
    if(qt>Pt) return true;
    if(qt>30.0) return true;
    if(x<.000129349) return true;

    return false;
}


double TMD::get_Xsection(double Pt, double qt, double z, double phi)
{
    double x = get_x(Pt, qt, z);
    if(Out_of_Kinematic_Constrains(Pt, qt, x)) return 0.0;

    double Q = get_Q(Pt);
    double epsf2= get_epsf2(Q, z);
    double Y = log(x0/x);

    double xG = get_xG_at_Y_qt(Y,qt);
    double xH = get_xH_at_Y_qt(Y,qt);

    double Conversion_prefactor= S_perp/S_perp_JIMWLK * 2 * M_PI; //one angle is integrated out

    double transverse_Xsection =
        prefactor * (z*z  + (1.0-z)*(1.0-z) ) * ( epsf2*epsf2 + pow(Pt,4) ) / pow(Pt*Pt+epsf2, 4) *
        ( xG - 2.0 * epsf2 * Pt * Pt/ ( epsf2*epsf2+pow(Pt,4) ) * cos (2.0*phi) * xH );

    double long_Xsection =
        prefactor * ( z*(1.0-z) ) * 8.0 * ( epsf2*pow(Pt,2) ) / pow(Pt*Pt+epsf2, 4) *
        ( xG + cos (2.0*phi) * xH );

    return Conversion_prefactor * (transverse_Xsection + long_Xsection) ;
}



void TMD::load_data(void)
{
    ifstream data ("misc.dat");


    for(size_t i=0; i<Pt_size; i++)
    {
        data>>Pt_arr[i];
    }

    for(size_t i=0; i<Y_size; i++)
    {
        data>>Y_arr[i];
    }

    for(size_t i=0; i<Pt_size*Y_size; i++)
    {
        data>>xG_arr[i] ;
    }

    for(size_t i=0; i<Pt_size*Y_size; i++)
    {
        data>>xH_arr[i];
    }
    data.close();
}

double TMD::get_xG_at_Y_qt(double Y, double qt)
{
    double result = interp2d_spline_eval(interp_xG,  qt, Y,  xa, ya);
    return result;
}

double TMD::get_xH_at_Y_qt(double Y, double qt)
{
    double result = interp2d_spline_eval(interp_xH, qt, Y,  xa, ya);
    return result;
}


class DiJetEvent
{
    RNGType *rng;
    TMD *Xsection;
    std::uniform_real_distribution<> *z_sample;
    std::uniform_real_distribution<> *Pt_sample;
    std::uniform_real_distribution<> *qt_sample;
    std::uniform_real_distribution<> *phi_sample;
    std::uniform_real_distribution<> *r_sample;

		static double f_minimize(const gsl_vector * x, void * params);
    double z_min, z_max, qt_min, qt_max, x0, Xsmax;
public:
    DiJetEvent(TMD* Xs);
    vector<double>  operator() (void);

};

vector<double> DiJetEvent::operator() (void)
{
    double r;
    double xs;

    double phi;
    double qt;
    double Pt;
    double z;


    do
    {
        r = (*r_sample)(*rng);
        phi=(*phi_sample)(*rng);
        z=(*z_sample)(*rng);
        qt=(*qt_sample)(*rng);
        Pt=(*Pt_sample)(*rng);
        xs = (*Xsection)(Pt, qt, z, phi) / Xsmax;
        //cout << xs << "\n";
        if (xs>1.0) cerr << "something wrong " << Xsmax << " " <<  xs<< "\n";
    }
    while(r>xs);
    vector<double> output;
    output.push_back(Pt);
    output.push_back(qt);
    output.push_back(z);
    output.push_back(phi);
    return output;
}


double DiJetEvent::f_minimize(const gsl_vector * x, void * params)
{
	double Pt = gsl_vector_get (x, 0);
	double qt = gsl_vector_get (x, 1);
	double z = gsl_vector_get (x, 2);


	DiJetEvent *instance = (DiJetEvent *) params;
	if(  (z-instance->z_min)*(z-instance->z_max) >0 ) return 0.0;
	if(  (qt-instance->qt_min)*(qt-instance->qt_max) >0 ) return 0.0;

	return ( - (*(instance->Xsection))(Pt, qt, z, 0.0) );
}

DiJetEvent::DiJetEvent(TMD* Xs): Xsection(Xs)
{
    z_min = 0.001;
    z_max = 1.0-0.001;
    qt_min = 1.0;
    qt_max = 10.0;

    x0 = Xsection->get_x0();


    const gsl_multimin_fminimizer_type *T =
        gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point */
    x = gsl_vector_alloc (3);
    gsl_vector_set (x, 0, qt_min);
    gsl_vector_set (x, 1, qt_min);
    gsl_vector_set (x, 2, 0.5);

    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (3);
    gsl_vector_set_all (ss, 0.1);

    /* Initialize method and iterate */
    minex_func.n = 3;
    minex_func.f = this->f_minimize;
    minex_func.params = this;

    s = gsl_multimin_fminimizer_alloc (T, 3);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);

        if (status == GSL_SUCCESS)
        {
            cerr<< "converged to maximim at \n";
        }

    }
    while (status == GSL_CONTINUE && iter < 100);

     Xsmax = - 2.0*s->fval ;


    cerr << Xsmax << "\n";

    unsigned long seed = mix(clock(), time(NULL), time(NULL));

    rng = new RNGType(seed);

    z_sample = new uniform_real_distribution<> ( z_min, z_max );
    phi_sample = new uniform_real_distribution<> ( 0.0, 2.0*M_PI );
    Pt_sample = new uniform_real_distribution<> ( qt_min, qt_max );
    qt_sample = new uniform_real_distribution<> ( qt_min, qt_max );
    r_sample = new uniform_real_distribution<> ( 0.0, 1.0 );

}


int main(int argc, char** argv)
{
    gsl_ieee_env_setup();
    TMD  generator = TMD(sqrts);
    DiJetEvent DJ (&generator);
    for(int i=0; i<number_of_events; i++)
    {
        vector<double>  event = DJ();
        cout << event[0] << " "  << event[1] << " "
             << event[2] << " "
             << event[3] << "\n"<<flush;
    }
    //6.33333 4.5 7.05232
    //cout << generator.get_xG_at_Y_qt(6.33333,4.5) << "\n";
    //cout << generator.get_xH_at_Y_qt(6.33333,4.5) << "\n";
}
