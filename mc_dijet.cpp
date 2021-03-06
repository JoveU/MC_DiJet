//==============================================================================
//  mc_dijet.cpp
//
//  Copyright (C) 2015 Adrian Dumitru and Vladimir Skokov
//
//  Please cited hep_number if used
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation.
//  This program is distributed in the hope that it will be useful,
//  but without any warranty; without even the implied warranty of
//  merchantability or fitness for a particular purpose. See the
//  GNU General Public License for more details.
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  Author: Adrian Dumitru and Vladimir Skokov
//  Last update:
//  $Date: November 15, 2015
//  $Authors: vladi@skokov.net
//==============================================================================



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
// Older versions of gcc may require explicit inclusion of the boost librariries
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/factorials.hpp>
*/

using namespace std;
typedef mt19937 RNGType;

const size_t number_of_events=1000000; // The total number of events to be generated
const double sqrts=100.0; // CM energy √s
const double alpha_em = 1.0/137.0; // α = e^2/4π
const int A_target=197; // The total number of nucleons in the nucleus target

const size_t Y_size=21, Pt_size=178; // The sizes of the input tables, which were generated by the external
// program solving JIMWLK evolution in Y. Do not change!


const double min_qt_to_Pt_ratio = 0.666;

const double qt_min = 1.0; // The minimum value of q⊥
const double qt_max = 30.0; // The maximum value of q⊥. Do not change due to finite range of the JIMWLK tables
const double Pt_max = 30.0; // The maximum value of \tilde{P⊥}
const double z_reg = 1e-3; // z lies in the range [z_reg,1-z_reg]. 0 and 1 are excluded to remove spurious
// divergent contributions proportional to 1/z and 1/(1-z)

const double epsrel = 1e-1; // The aimed relative error for computing the numericl integrals
constexpr double X0 = 1e-2; // The maximal value of Bjorken x and x_g

// Below I the parameters for the Gold Nucleus are defined
const double S_perp0 = 1300; // The transverse area S⊥ in mb
const double Qs0 = 1; // The saturation scale in GeV
const int A0 = 197; // The number of nucleons

const double M = 1; //The nuclon mass

unsigned long seed; // The seed and the pointer to the random number generator object.
// The initialization is in main
RNGType *rng;

double yf(double Q2, double W2, double S)
/*
 y = (W^2 - M^2 + Q^2) / (s - M^2)
 */
{
    return (W2-M*M+Q2)/(S-M*M);
}


double Flux_T(double Q2, double W2, double S)
/*
 The photonic flux F(Q^2, W^2, S) for the transverse polariazation, see Eq. (???)
 */
{
    double y = yf(Q2,W2,S);
    double redy = 1-y;
    return  alpha_em/(2*M_PI*Q2*S*y) * (1+redy*redy) ;
}

double Flux_L(double Q2, double W2, double S)
/*
 The photonic flux F(Q^2, W^2, S) for the longitudonal polariazation, see Eq. (???)
 */
{
    double y = yf(Q2,W2,S);
    double redy = 1-y;
    return  alpha_em/(M_PI*Q2*S*y) * redy  ;
}


unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
/*
 Somewhat complicated procedure to generate a random seed.
 It makes possible the massive submission of the several instances of the
 program to a cluster. The seed is guranteed to be different even if the time
 of exceution is (almost) the same.
 */
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
/*
 The class computing the x-section for given P⊥, q⊥, z, φ and the polarization. The inputs are W, Q, A.
 sqrtSin is W.
 */
{
    interp2d_spline* interp_xG; //The interpolating function for xG(Y, q⊥).
    interp2d_spline* interp_xH; //The interpolating function for xH(Y, q⊥).
    gsl_interp_accel *xa, *ya;
    const interp2d_type* T;

    double* Y_arr; // The array of Y's; read from the external file
    double* Pt_arr; // The array of q⊥'s; read from the external file
    double* xG_arr; // The array of xG; red from the external file
    double* xH_arr; // The array of xH; red from the external file

    double sqrtS; // W
    double s; // W^2
    double prefactor; //Prefactor that includes α_em α_s ∑_f q_f^2, see Eqs (1-2) in 1508.04438
    double Q; // Q
    double S_perp; // The transverse area in mb
    double Qs; // The saturation scele in GeV

    constexpr static double sum_charge2 = .666666; // The sum of quark electric charges squared; ∑_f q_f^2
                                          constexpr static double alpha_s = .15; // α_s Can be changed.
    constexpr static double S_perp_JIMWLK =  2704.0; // Do not change! The conversion factor.
    // JIMWLK is computed on a lattice of a certain transverse size.

    void load_data(void); //The function to load the external data file from the JIMWLK evolution
    double get_Q(double Pt); // Return Q (it is Pt independent, the Pt is kept due to historical reasons)
    double get_epsf2(double Q, double z); //ϵ_f^2 = z(1-z) Q^2, see the top of the page 2 in 1508.04438
    double get_Xsection(double Pt, double qt, double z, double phi, int pol); // Returns the crossection
    bool Out_of_Kinematic_Constrains(double Pt, double qt, double x); // Retruns true, if the variables are
    //out of the kinematical constraints and violate the "correlations limit"


public:
    constexpr static  double x0 = X0;
    TMD(double sqrtSin, double Qin, int A); // The constructor
		virtual ~TMD(); //The destructor 
    vector<double> get_Xsection_components(double Pt, double qt, double z, double phi);
    double get_x0(void)
    {
        return x0;
    }

    double get_xG_at_Y_qt(double Y, double qt); // Returns xG at Y and q⊥
    double get_xH_at_Y_qt(double Y, double qt);  // Returns xH at Y and q⊥
    double operator() (double Pt, double qt, double z, double phi, int pol) //The operator () introduced
    // to simplifie the code, returns the x-section.
    {
        return get_Xsection(Pt,qt,z,phi,pol);
    }

    double get_x(double Pt, double qt, double z); // Returns x, see Eq. 5 in 1508.04438, note that s->W^2.
};

TMD::TMD(double sqrtSin, double Qin, int A): sqrtS(sqrtSin), Q(Qin)
{
    S_perp = S_perp0 * pow( ((double) A)/((double) A0), 2.0/3.0); // Sperp is defined wrt Au nucleus
    //and scaled accordingly with the number of nucleons
    Qs = Qs0 * pow( ((double) A)/((double) A0), 1.0/6.0); // Qs is 1 GeV for the cold nucleus,
    // otherwise scaled accordign to A
    s = sqrtS*sqrtS; // s is W^2; introduced to make the code run faster
    prefactor = alpha_em * alpha_s * sum_charge2;

    // Begin: Initialization of the 2d bilinear interpolating functions, see the latest GSL manual.
    T=interp2d_bilinear;
    xa = gsl_interp_accel_alloc();
    ya = gsl_interp_accel_alloc();
    interp_xG = interp2d_spline_alloc(T, Pt_size, Y_size);
    interp_xH = interp2d_spline_alloc(T, Pt_size, Y_size);

    Y_arr = new double[Y_size];
    Pt_arr = new double[Pt_size];
    xG_arr = new double[Pt_size*Y_size];
    xH_arr = new double[Pt_size*Y_size];

    load_data();// Loads the external file

    interp2d_spline_init(interp_xG, Pt_arr, Y_arr, xG_arr, Pt_size, Y_size);
    interp2d_spline_init(interp_xH, Pt_arr, Y_arr, xH_arr, Pt_size, Y_size);
    // End: Initialization of the 2d bilinear interpolating functions, see the latest GSL manual.
}

TMD::~TMD()
{
	delete(Y_arr);
	delete(Pt_arr);
	delete(xG_arr);
	delete(xH_arr);
	interp2d_spline_free(interp_xG);
	interp2d_spline_free(interp_xH);
	gsl_interp_accel_free(xa);
	gsl_interp_accel_free(ya);
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
    if(qt>min_qt_to_Pt_ratio*Pt) return true; // This constraints ensures that only event with q⊥<P⊥*min_qt_to_Pt_ratio  are generated
    if(qt>qt_max) return true;
    if(x<.000129349) return true; //JIMWLK smallest x

    return false;
}


double TMD::get_Xsection(double Pt, double qt, double z, double phi, int pol)
// Returns X-section for given polarization
{
    double x = get_x(Pt, qt, z);
    if(Out_of_Kinematic_Constrains(Pt, qt, x)) return 0.0;

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
// Returns X-section's as components of the vector, the zeros component is te transverse polarization.
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
    ifstream data ("misc.dat"); //External file. With the results from JIMWLK.


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
/*
The class generating random events that correspond to the specific x-section's for given W^2 and Q^2.
 */
{

    TMD *Xsection;
    //Begin: Random number generators
    std::uniform_real_distribution<> *z_sample;
    std::uniform_real_distribution<> *Pt_sample;
    std::uniform_real_distribution<> *qt_sample;
    std::uniform_real_distribution<> *phi_sample;
    std::uniform_real_distribution<> *r_sample;
    //End: Random number generators

    static double f_minimize(const gsl_vector * x, void * params); //Generates the random numbers we need
    // the maximum value of the distribution; this function is a helper function, returns the x-section
    double z_min, z_max, x0, Xsmax;
public:
    DiJetEvent(TMD* Xs); // The constructor; a pointer to a TMD object has to be provided
		~DiJetEvent(); // The destructor
    vector<double>  operator() (int pol); // Returns a vector of all random kinematic variables, e.g. q⊥, P⊥
    vector<double> k1k2f(vector<double>  params); // Computes k1 and k2, from Eq 3 of 1508.04438

};

vector<double> DiJetEvent::k1k2f(vector<double>  params)
// Computes k1 and k2, from Eq 3 of 1508.04438
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
/* Returns a vector of all random kinematic variables, e.g. q⊥, P⊥
 for a given polarization
 */
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

    return ( - (*(instance->Xsection))(Pt, qt, z, 0.0,2) ); //Returns negative x-section to be minimized
}

DiJetEvent::DiJetEvent(TMD* Xs): Xsection(Xs)
{
    z_min =  z_reg;
    z_max = 1.0-z_reg;

    x0 = Xsection->get_x0();


    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point for the minimization of the nexatife x-section*/
    x = gsl_vector_alloc (3);
    gsl_vector_set (x, 0, qt_min/min_qt_to_Pt_ratio);
    gsl_vector_set (x, 1, qt_min);
    gsl_vector_set (x, 2, 0.5);

    /* Set initial step sizes to 0.1 */
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

    Xsmax = - 2.0*s->fval ; //With Xsmax we can generate the events.

    //cerr << Xsmax << "\n";
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);

    //Begin: Initialization of random number generators
    z_sample = new uniform_real_distribution<> ( z_min, z_max );
    phi_sample = new uniform_real_distribution<> ( 0.0, 2.0*M_PI );
    Pt_sample = new uniform_real_distribution<> ( qt_min, qt_max );
    qt_sample = new uniform_real_distribution<> ( qt_min, qt_max );
    r_sample = new uniform_real_distribution<> ( 0.0, 1.0 );
    //End: Initialization of random number generators
}

DiJetEvent::~DiJetEvent()
{
	delete(z_sample);
	delete(phi_sample);
	delete(Pt_sample);
	delete(qt_sample);
	delete(r_sample);
}


struct tmd_parameters
/* To perform GSL integration of the X-section we need to pass the parameters.
 It can be conveniently done using this container;
 */
{
    TMD* gen;
    int* pol;
    double* Pt;
    double* qt;
    double* z;
};


class DIS
/*
 This is DIS part of the generator; it randomly finds W^2 and Q^2 weighted by the TMD x-section
 and photonic flows.
 */
{
    std::uniform_real_distribution<> *Q2_sample;
    std::uniform_real_distribution<> *W2_sample;
    std::uniform_real_distribution<> *x_sample;
    std::uniform_real_distribution<> *r_sample;


    double sqrtS; // "real" √s
    double Q2_min, Q2_max; // minimum Q^2 and maximum Q^2
    double x0; // minimum Bjorken x
    int  A; // The atomic number
    double S;

    //Arrays for the interpolations and integration
    double* Xs_L;
    double* Xs_T;
    double* vec_Q2;
    double* vec_W2;

    //Integrated X-sections (in the defined kinematical regime)
    double IntXS_T, IntXS_L;

    size_t ind_Q2, ind_W2;

    //Interpolating functions for the transverse and longitudonal x-sections, see Eqs. (???)
    interp2d_spline* interp_Xs_L;
    interp2d_spline* interp_Xs_T;

    gsl_interp_accel *xa, *ya;
    const interp2d_type* T;

    TMD* generator;
    DiJetEvent* DJ;

    double integrated_Xs(double Q, double W, int pol);
public:
    DIS(double sqrtSin, int Ain); // The constructor. Parameters √s and A.
    ~DIS(); // The destructor 
		vector<double> operator() (void);
};



double u_Int_qt(double qt, void* params)
/*
 Helper function to integrate the x-sections
 */
{
    double Pt = *(((struct tmd_parameters *) params )->Pt);
    double z = *(((struct tmd_parameters *) params )->z);
    int pol = *(((struct tmd_parameters *) params )->pol);
    TMD* gen  = ((struct tmd_parameters *) params )->gen;

    vector<double> X_section = gen->get_Xsection_components(Pt, qt, z, M_PI/4.0); //So that cos(2 phi) = 0  --  the angular integration results in this

    return X_section.at(pol);
}


double u_Int_Pt(double Pt, void* params)
/*
 Helper function to integrate the x-sections
 */
{
    double result, error;
    gsl_function F;

    F.function = &u_Int_qt;

    tmd_parameters  pass = *(struct tmd_parameters *) params;
    pass.Pt = &Pt;
    F.params = &pass;

    gsl_integration_workspace * gsl_int_ws = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag(&F, qt_min, Pt*min_qt_to_Pt_ratio, 0, epsrel, 1000, 1, gsl_int_ws, &result, &error); //the upper integration limit is Pt
    gsl_integration_workspace_free (gsl_int_ws);

    return result;
}


double u_Int_z(double z, void* params)
/*
 Helper function to integrate the x-sections
 */
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
/*
 Integrated X_section including the photonic fluxes for given Q, W and the polarization. This is required
 to randomly assign the polarization to the photon.
 */
{
    TMD generator = TMD(W, Q, A); //create the TMD class

    double result, error;
    gsl_function F;

    F.function = &u_Int_z;

    tmd_parameters  pass;
    pass.gen = &generator;
    pass.pol = &pol;
    F.params = &pass;

    gsl_integration_workspace * gsl_int_ws = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag(&F, z_reg, 1.0-z_reg, 0, epsrel, 1000, 1, gsl_int_ws, &result, &error);
    gsl_integration_workspace_free (gsl_int_ws);

    return result;
}

DIS::DIS(double sqrtSin, int Ain):sqrtS(sqrtSin),A(Ain)
{
    S = sqrtS*sqrtS;

    ind_Q2 = 10; // The numerical integration of the x-sections cannot be done for all Q and W.
    ind_W2 = 10; // Thus it will be only computed at ind_Q2xind_W2 points and then extrapolated between those.

    Xs_L = new double[ind_Q2*ind_W2]; //The matrices to be populated with the integrated x-section
    Xs_T = new double[ind_Q2*ind_W2];

    vec_Q2 = new double[ind_Q2];
    vec_W2 = new double[ind_W2];

    double Q2_max = (S-M*M)*X0/(1.0-X0); // See Eqs. (???)
    Q2_min = 4;

    double Q2_step = (Q2_max-Q2_min)/(ind_Q2-1);


    double W2_max = S;
    double W2_min = M*M+Q2_min*(1.0-X0)/X0;
    double W2_step = (W2_max-W2_min)/(ind_W2-1);

    cout << "# generating interpolating cache\n";

    IntXS_T =0.0;
    IntXS_L =0.0;

    for(int i=0; i<ind_Q2; i++)
    {
        double Q2 = Q2_min + Q2_step*i;
        vec_Q2[i] =  Q2;
        for(int j=0; j<ind_W2; j++)
        {
            double W2 = W2_min + W2_step*j;
            vec_W2[j] =  W2;
            double Q = sqrt(Q2);
            double W = sqrt(W2);
            Xs_T[i+j*ind_Q2] = Flux_T(Q2, W2, S) * integrated_Xs(Q, W, 0);
            Xs_L[i+j*ind_Q2] = Flux_L(Q2, W2, S) * integrated_Xs(Q, W, 1);

            IntXS_T+=Xs_T[i+j*ind_Q2]*W2_step*Q2_step;
            IntXS_L+=Xs_L[i+j*ind_Q2]*W2_step*Q2_step;
        }
    }
    cout << "# done generating interpolating cache\n";
    cout << "# sqrt(s) A  integrated x_sections\n";
    cout << sqrtS << " " << A << " " << IntXS_T << " " <<  IntXS_L << "\n"<<flush;

    T=interp2d_bilinear;
    xa = gsl_interp_accel_alloc();
    ya = gsl_interp_accel_alloc();
    interp_Xs_L = interp2d_spline_alloc(T, ind_Q2, ind_W2);
    interp_Xs_T = interp2d_spline_alloc(T, ind_Q2, ind_W2);

    interp2d_spline_init(interp_Xs_L, vec_Q2, vec_W2, Xs_L, ind_Q2, ind_W2);
    interp2d_spline_init(interp_Xs_T, vec_Q2, vec_W2, Xs_T, ind_Q2, ind_W2);

    //Initialization of the random number generators
    Q2_sample = new uniform_real_distribution<> ( Q2_min, Q2_max );
    r_sample = new uniform_real_distribution<> ( 0, 1 );
}


DIS::~DIS()
{
	delete(Q2_sample);
	delete(r_sample);
  gsl_interp_accel_free(xa);
  gsl_interp_accel_free(ya);
	interp2d_spline_free(interp_Xs_L);
	interp2d_spline_free(interp_Xs_T);
  delete(Xs_L);
  delete(Xs_T);
  delete(vec_Q2);
  delete(vec_W2);
}

vector<double> DIS::operator() (void)
// Returns an event
{

    double Q2, x, W2, r, Xs_L, Xs_T, ratio;
    double Q, W;
    bool longit;
    int pol;
    TMD*  generator;
    DiJetEvent* DJ;
    vector<double> event;
    vector<double> k1k2; 

    double W2_max = S;
    do
    {
        Q2 = (*Q2_sample)(*rng);
        
				double W2_min = M*M+Q2_min*(1.0-X0)/X0;

        W2_sample = new uniform_real_distribution<> ( W2_min, W2_max );
        
				W2 = (*W2_sample)(*rng);
        r = (*r_sample)(*rng);

        Xs_L =  interp2d_spline_eval(interp_Xs_L,  Q2, W2,  xa, ya);
        Xs_T =  interp2d_spline_eval(interp_Xs_T,  Q2, W2,  xa, ya);

        ratio  = Xs_L/( Xs_L + Xs_T );

        longit = false;
        if(r<ratio) longit = true;

        //Step 5 from the notes

        Q = sqrt(Q2);
        W = sqrt(W2);

        TMD*  generator = new TMD(W,Q,A);
        DJ = new DiJetEvent (generator);

        pol = 0; //transverse
        if (longit) pol=1; //longit

        event = (*DJ)(pol);
    		if(event.size()>1) k1k2 = DJ->k1k2f(event);
    		delete(W2_sample);
    		delete(generator);
    		delete(DJ);
    }
    while(event.size()<1);

    vector<double> output = event;
    
		output.push_back(W);
    output.push_back(Q);
    output.push_back(pol);
    
		for(int i=0; i<k1k2.size(); i++)
    {
        output.push_back(k1k2.at(i));
    }

    return output;
}

int main(int argc, char** argv)
{
    gsl_ieee_env_setup();

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

}
