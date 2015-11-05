#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_integration.h>
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
typedef mt19937 RNGType;

const size_t number_of_events=1000000;
const double sqrts=100.0;
const int A_target=197;


//input tables
const size_t Y_size=21, Pt_size=178;


const double qt_min = 1.0;
const double qt_max = 30.0;
const double Pt_max = 30.0;
const double z_reg = 1e-3;
const double epsrel = 1e-1;
constexpr double X0 = 1e-2;


//Target: Parameters we use for Au
const double S_perp0 = 1300;
const double Qs0 = 1; // in GeV
const int A0 = 197;
//

const double M = 1; //nuclon mass

unsigned long seed;
RNGType *rng;

double Flux_T(double Q, double W)
{
    return 1.0;
}

double Flux_L(double Q, double W)
{
    return 1.0;
}




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
    double Q;
    double S_perp; //mb
    double Qs; //mb
    //constexpr static double ratio_Q_over_Pt = 2.0;

    constexpr static double sum_charge2 = .666666;
    constexpr static double alpha_em = .00729927;
    constexpr static double alpha_s = .15;
    constexpr static double S_perp_JIMWLK =  2704.0; //mb

    void load_data(void);
    double get_Q(double Pt);
    double get_epsf2(double Q, double z);
    double get_Xsection(double Pt, double qt, double z, double phi, int pol);
    bool Out_of_Kinematic_Constrains(double Pt, double qt, double x);


public:
    constexpr static  double x0 = X0;
    TMD(double sqrtSin, double Qin, int A);
    vector<double> get_Xsection_components(double Pt, double qt, double z, double phi);
    double get_x0(void)
    {
        return x0;
    }

    double get_xG_at_Y_qt(double Y, double qt);
    double get_xH_at_Y_qt(double Y, double qt);
    double operator() (double Pt, double qt, double z, double phi, int pol)
    {
        return get_Xsection(Pt,qt,z,phi,pol);
    }

    double get_x(double Pt, double qt, double z);
};

TMD::TMD(double sqrtSin, double Qin, int A): sqrtS(sqrtSin), Q(Qin)
{
    S_perp = S_perp0 * pow( ((double) A)/((double) A0), 2.0/3.0);
    Qs = Qs0 * pow( ((double) A)/((double) A0), 1.0/6.0);
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
    if(qt>qt_max) return true;
    if(x<.000129349) return true;

    return false;
}


double TMD::get_Xsection(double Pt, double qt, double z, double phi, int pol)
{
    double x = get_x(Pt, qt, z);
    if(Out_of_Kinematic_Constrains(Pt, qt, x)) return 0.0;

    //double Q = get_Q(Pt);
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
    if (pol==0) return  Pt*qt*Conversion_prefactor * (transverse_Xsection) ;
    if (pol==1) return  Pt*qt*Conversion_prefactor * (long_Xsection);
    if (pol==2) return  Pt*qt*Conversion_prefactor * (long_Xsection+transverse_Xsection);
    return 0.0;
}

vector<double> TMD::get_Xsection_components(double Pt, double qt, double z, double phi)
{
    double x = get_x(Pt, qt, z);
    vector<double> out;
    if(Out_of_Kinematic_Constrains(Pt, qt, x))
    {
        out.push_back(0.0);
        out.push_back(0.0);
        return out;
    }

    //double Q = get_Q(Pt);
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

    out.push_back(Pt*qt*Conversion_prefactor * transverse_Xsection);
    out.push_back(Pt*qt*Conversion_prefactor * long_Xsection);
    return out;
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
    double result = interp2d_spline_eval(interp_xG,  qt/Qs, Y,  xa, ya);
    return result;
}

double TMD::get_xH_at_Y_qt(double Y, double qt)
{
    double result = interp2d_spline_eval(interp_xH, qt/Qs, Y,  xa, ya);
    return result;
}


class DiJetEvent
{
    TMD *Xsection;
    std::uniform_real_distribution<> *z_sample;
    std::uniform_real_distribution<> *Pt_sample;
    std::uniform_real_distribution<> *qt_sample;
    std::uniform_real_distribution<> *phi_sample;
    std::uniform_real_distribution<> *r_sample;

    static double f_minimize(const gsl_vector * x, void * params);
    double z_min, z_max, x0, Xsmax;
public:
    DiJetEvent(TMD* Xs);
    vector<double>  operator() (int pol);
    vector<double> k1k2f(vector<double>  params);

};

vector<double> DiJetEvent::k1k2f(vector<double>  params)
{
    double Pt = params.at(0);
    double qt = params.at(1);
    double z = params.at(2);
    double phi = params.at(3);
    double phi_pt = params.at(4);

    double Ptx = Pt * cos(phi_pt);
    double Pty = Pt * sin(phi_pt);

    double qtx = qt * cos(phi_pt+phi);
    double qty = qt * sin(phi_pt+phi);

    double k1x =  Ptx + z*qtx;
    double k1y =  Pty + z*qty;

    double k2x =  -Ptx + (1.0-z)*qtx;
    double k2y =  -Pty + (1.0-z)*qty;

    vector<double> k1k2;
    k1k2.push_back(k1x);
    k1k2.push_back(k1y);
    k1k2.push_back(k2x);
    k1k2.push_back(k2y);

    return k1k2;
}

vector<double> DiJetEvent::operator() (int pol)
{
    double r;
    double xs;

    double phi;
    double phi_Pt;
    double qt;
    double Pt;
    double z;

    size_t count=0;

    do
    {
        r = (*r_sample)(*rng);
        phi=(*phi_sample)(*rng);
        phi_Pt=(*phi_sample)(*rng);
        z=(*z_sample)(*rng);
        qt=(*qt_sample)(*rng);
        Pt=(*Pt_sample)(*rng);
        xs = (*Xsection)(Pt, qt, z, phi, pol) / Xsmax;
        //cout << xs << "\n";
        count++;
        if (xs>1.0) cerr << "something wrong " << Xsmax << " " <<  xs<< "\n";
        if((count>500000)||(std::isnan(xs)))
        {
            vector<double> out;
            //cerr << "it is taking too long for this configuration \n";
            return out; //empty vector
        }
    }
    while(r>xs);
    vector<double> output;
    output.push_back(Pt);
    output.push_back(qt);
    output.push_back(z);
    output.push_back(phi);
    output.push_back(phi_Pt);
    output.push_back(xs*Xsmax);
    return output;
}


double DiJetEvent::f_minimize(const gsl_vector * x, void * params)
{
    double Pt = gsl_vector_get (x, 0);
    double qt = gsl_vector_get (x, 1);
    double z = gsl_vector_get (x, 2);


    DiJetEvent *instance = (DiJetEvent *) params;
    if(  (z-instance->z_min)*(z-instance->z_max) >0 ) return 0.0;
    if(  (qt-qt_min)*(qt-qt_max) >0 ) return 0.0;

    return ( - (*(instance->Xsection))(Pt, qt, z, 0.0,2) );
}

DiJetEvent::DiJetEvent(TMD* Xs): Xsection(Xs)
{
    z_min =  z_reg;
    z_max = 1.0-z_reg;

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
        /*
        if (status == GSL_SUCCESS)
        {
        cerr<< "converged to maximim at \n";
        }
        */

    }
    while (status == GSL_CONTINUE && iter < 10000);

    Xsmax = - 2.0*s->fval ;


    //cerr << Xsmax << "\n";



    z_sample = new uniform_real_distribution<> ( z_min, z_max );
    phi_sample = new uniform_real_distribution<> ( 0.0, 2.0*M_PI );
    Pt_sample = new uniform_real_distribution<> ( qt_min, qt_max );
    qt_sample = new uniform_real_distribution<> ( qt_min, qt_max );
    r_sample = new uniform_real_distribution<> ( 0.0, 1.0 );

}


struct tmd_parameters
{
    TMD* gen;
    int* pol;
    double* Pt;
    double* qt;
    double* z;
};


class DIS
{
    std::uniform_real_distribution<> *Q_sample;
    std::uniform_real_distribution<> *x_sample;
    std::uniform_real_distribution<> *r_sample;


    double sqrtS;
    double Q_min, Q_max;
    double x0;
    int  A;

    double* Xs_L;
    double* Xs_T;
    double* vec_Q;
    double* vec_W;

    size_t ind_Q, ind_W;

    interp2d_spline* interp_Xs_L;
    interp2d_spline* interp_Xs_T;

    gsl_interp_accel *xa, *ya;
    const interp2d_type* T;


    TMD* generator;
    DiJetEvent* DJ;

    double integrated_Xs(double Q, double W, int pol);
public:
    DIS(double sqrtSin, int Ain);
    vector<double> operator() (void);
};



double u_Int_qt(double qt, void* params)
{
    double Pt = *(((struct tmd_parameters *) params )->Pt);
    double z = *(((struct tmd_parameters *) params )->z);
    int pol = *(((struct tmd_parameters *) params )->pol);
    TMD* gen  = ((struct tmd_parameters *) params )->gen;

    vector<double> X_section = gen->get_Xsection_components(Pt, qt, z, M_PI/4.0); //so that cos(2 phi) = 0  --  the angular integration results in this


    //cout << z << " " << qt << " " << X_section.at(pol) << "\n";
    return X_section.at(pol);
}


double u_Int_Pt(double Pt, void* params)
{
    double result, error;
    gsl_function F;

    F.function = &u_Int_qt;

    tmd_parameters  pass = *(struct tmd_parameters *) params;
    pass.Pt = &Pt;
    F.params = &pass;

    gsl_integration_workspace * gsl_int_ws = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag(&F, qt_min, Pt, 0, epsrel, 1000, 1, gsl_int_ws, &result, &error); //the upper integration limit is Pt
    gsl_integration_workspace_free (gsl_int_ws);

    return result;
}


double u_Int_z(double z, void* params)
{
    double result, error;
    gsl_function F;

    F.function = &u_Int_Pt;

    //cout << z << "zz \n";

    tmd_parameters  pass = *(struct tmd_parameters *) params;
    pass.z = &z;
    //cout << *pass.z << " " << z <<  " z \n";
    F.params = &pass;
    gsl_integration_workspace * gsl_int_ws = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag(&F, qt_min, Pt_max, 0, epsrel, 1000, 1, gsl_int_ws, &result, &error); //the upper integration limit is Pt
    gsl_integration_workspace_free (gsl_int_ws);

    return result;
}

double DIS::integrated_Xs(double Q, double W, int pol)
{
    TMD generator = TMD(W, Q, A); //create the TMD class

    double result, error;
    gsl_function F;

    F.function = &u_Int_z;

    tmd_parameters  pass;
    pass.gen = &generator;
    pass.pol = &pol;
    F.params = &pass;

    //cout << z_reg << " reg \n";

    gsl_integration_workspace * gsl_int_ws = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag(&F, z_reg, 1.0-z_reg, 0, epsrel, 1000, 1, gsl_int_ws, &result, &error); //the upper integration limit is Pt
    gsl_integration_workspace_free (gsl_int_ws);

    return result;
}

DIS::DIS(double sqrtSin, int Ain):sqrtS(sqrtSin),A(Ain)
{
    ind_Q = 50;
    ind_W = 50;

    Xs_L = new double[ind_Q*ind_W]; //matrices to be populated with the integrated x-section
    Xs_T = new double[ind_Q*ind_W];

    vec_Q = new double[ind_Q];
    vec_W = new double[ind_W];

    double Q_max = sqrt(sqrtS*sqrtS-M*M)*sqrt(X0);
    double Q_min = 2;
    double Q_step = (Q_max-Q_min)/(ind_Q-1);

    double W_max = sqrtS;
    double W_min = sqrt(M*M+Q_min*Q_min*(1.0-X0)/X0);
    double W_step = (W_max-W_min)/(ind_W-1);

    cout << "# generating interpolating cache\n";
    for(int i=0; i<ind_Q; i++)
    {
        double Q = Q_min + Q_step*i;
        vec_Q[i] =  Q;
        for(int j=0; j<ind_W; j++)
        {
            double W = W_min + W_step*j;
            vec_W[j] =  W;
            Xs_T[i+j*ind_Q] = Flux_T(Q, W) * integrated_Xs(Q, W, 0);
            Xs_L[i+j*ind_Q] = Flux_L(Q, W) * integrated_Xs(Q, W, 1);
            //cout << Q << " " <<  W  <<  " "  << Xs_T[i+j*ind_Q] << " " << Xs_L[i+j*ind_Q] <<   " \n";
        }
    }
    cout << "# done generating interpolating cache\n";


    T=interp2d_bilinear;
    xa = gsl_interp_accel_alloc();
    ya = gsl_interp_accel_alloc();
    interp_Xs_L = interp2d_spline_alloc(T, ind_Q, ind_W);
    interp_Xs_T = interp2d_spline_alloc(T, ind_Q, ind_W);

    interp2d_spline_init(interp_Xs_L, vec_Q, vec_W, Xs_L, ind_Q, ind_W);
    interp2d_spline_init(interp_Xs_T, vec_Q, vec_W, Xs_T, ind_Q, ind_W);

    //unsigned long seed = mix(clock(), time(NULL), time(NULL));
    //rng = new RNGType(seed);
    Q_sample = new uniform_real_distribution<> ( Q_min, Q_max );
    r_sample = new uniform_real_distribution<> ( 0, 1 );
}


struct event_output
{

};

vector<double> DIS::operator() (void)
{

    double Q, x, W, r, Xs_L, Xs_T, ratio;
    bool longit;
    int pol;
    TMD*  generator;
    DiJetEvent* DJ;
    vector<double> event;

    do
    {
        Q = (*Q_sample)(*rng);

        x_sample = new uniform_real_distribution<> ( Q*Q/(sqrtS*sqrtS-M*M), X0);
        x = (*x_sample)(*rng);

        W = sqrt(M*M + (1.0-x)*Q*Q/x);

        r = (*r_sample)(*rng);

        Xs_L =  interp2d_spline_eval(interp_Xs_L,  Q, W,  xa, ya);
        Xs_T =  interp2d_spline_eval(interp_Xs_T,  Q, W,  xa, ya);

        ratio  = Xs_L/( Xs_L + Xs_T );

        longit = false;
        if(r<ratio) longit = true;

        //Step 5 from the notes


        TMD*  generator = new TMD(W,Q,A);
        DJ = new DiJetEvent (generator);

        pol = 0; //transverse
        if (longit) pol=1; //longit

        event = (*DJ)(pol);
    }
    while(event.size()<1);

    vector<double> k1k2 = DJ->k1k2f(event);

    vector<double> output = event;
    output.push_back(W);
    output.push_back(Q);
    output.push_back(pol);
    for(int i=0; i<k1k2.size(); i++)
    {
        output.push_back(k1k2.at(i));
    }

    delete(x_sample);

    return output;
}

int main(int argc, char** argv)
{
    gsl_ieee_env_setup();
    //TMD  generator = TMD(sqrts);
    //DiJetEvent DJ (&generator);
    seed = mix(clock(), time(NULL), time(NULL));
    rng = new RNGType(seed);

    DIS generator = DIS(sqrts, A_target);
    cout << "#  Pt\t qt\t z\t phi\t phiT\t xS\t W\t Q\t Pol\t k1x\t k1y\t k2x\t k2y\n";
    for(int i=0; i<number_of_events; i++)
    {
        vector<double> event = generator();
        for(int j=0; j<event.size(); j++)
        {
            cout << event[j] << " ";
        }
        cout << "\n" << flush;
    }
    //6.33333 4.5 7.05232
    //cout << generator.get_xG_at_Y_qt(6.33333,4.5) << "\n";
    //cout << generator.get_xH_at_Y_qt(6.33333,4.5) << "\n";


}
