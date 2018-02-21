#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>
//#include <acado/matrix_vector/matrix.hpp>
//#include <cmath>
//#include <stdlib.h>
//#include <iostream>

int main() {

USING_NAMESPACE_ACADO

    // Constants
    const double     g = 9.81;  // the gravitational constant [m/s^2]
    const double     pi = 3.141592; // the constant PI

    const double     a = 1; // weight factor of cost function
    const double     b = 1; // weight factor of velocity change
    const double     c = 1; // weight factor of bank angle rate change

    const double     alpha = 0.9999; // weight factor of received power vs. fuel consumption
    const double     beta = 1000; // fitting constant of the power error function

    const double     v_a_min = 8; // min air speed [m/s]
    const double     v_a_max = 18; // max air speed [m/s]

    const double     phi_min = -0.7; // min roll [rad]
    const double     phi_max = 0.7; // max roll [rad]

    const double     dphi_min = -1.4; // min roll angle rate [rad/s]
    const double     dphi_max = 1.4; // max roll angle rate [rad/s]

    const double     r_c = 1; // min distance between UAVs

    const double     h0 = 0; // vessel negative altitude [m]

    const double     h1 = -100; // UAV 1 negative altitude [m]

    const double     hN = -20; // ground station negative altitude [m]
    const double     xN = 0; // ground station x [m]
    const double     yN = 0; // ground station y [m]

    const double     P_d = 0.25;//100*100/(100*100); // Power at ~100 m

    const double     v_w = 3;
    const double     psi_w = 0.7;

    const double     G0 = 10;
    const double     G1 = 10;
    const double     GN = 10;

    const double V_batt           = 14.8;    // Volts
    const double V_init           = 14.8;    
    const double W_start_glide    = 50;      // Weight in N
    const double rho              = 1.225;   // kg/m3
    const double S                = 0.75;    // AC surface main wing
    const double eta_p_bldc       = 0.5;     // Propulsion efficiency ESC + motor + cable + discharge losses
    
    //
    DifferentialState x1, y1, v_a1, psi1, phi1;
    DifferentialState x0, y0;

    Control u_phi1, u_v1;
    
    DifferentialEquation f;
    OCP ocp(0.0, 30, 90);

	Expression d01, d1N, D01, D10, D1N, DN1, P_f01, P_f1N, P01N_min, E1, F_u1, F_t, F1, J1, L1;
    Expression E_tmp, P_cons;

    E_tmp << W_start_glide / (S * rho * v_a1*v_a1);
    P_cons << V_batt * S * rho * v_a1*v_a1*v_a1 * (2032*E_tmp*E_tmp - 230*E_tmp + 89) / (20000*V_init*eta_p_bldc);

    d01 << sqrt((h1-h0)*(h1-h0) + (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    d1N << sqrt((hN-h1)*(hN-h1) + (xN-x1)*(xN-x1) + (yN-y1)*(yN-y1));

    D01 << G0 * (10 - acos(-(h1-h0) / d01)*acos(-(h1-h0) / d01));
    D10 << G1 * (10 - acos(((x1-x0)*(-sin(phi1)*sin(psi1)) + (y1-y0)*(sin(phi1)*cos(psi1)) + (h1-h0)*(-cos(phi1))) / d01)*acos(((x1-x0)*(-sin(phi1)*sin(psi1)) + (y1-y0)*(sin(phi1)*cos(psi1)) + (h1-h0)*(-cos(phi1))) / d01));
    D1N << G1 * (10 - acos(((x1-xN)*(-sin(phi1)*sin(psi1)) + (y1-yN)*(sin(phi1)*cos(psi1)) + (h1-hN)*(-cos(phi1))) / d1N)*acos(((x1-xN)*(-sin(phi1)*sin(psi1)) + (y1-yN)*(sin(phi1)*cos(psi1)) + (h1-hN)*(-cos(phi1))) / d1N));
    DN1 << GN * (10 - acos(-(hN-h1) / d1N)*acos(-(hN-h1) / d1N));
    
    P_f01 << D10 * D01 / (d01*d01); // with constants removed
    P_f1N << DN1 * D1N / (d1N*d1N); // with constants removed

    P01N_min << 0.5 * (P_f01 + P_f1N - sqrt((P_f01-P_f1N)*(P_f01-P_f1N)));
    E1 << (P_d - P01N_min) * ( pi/2 - atan(beta*(P01N_min-P_d)) ) / pi;
    //E1 << P_d - P01N_min;
    //E1 << d01 + d1N;
    //F_u1 << 0.1 * (v_a1 - 12) * (v_a1 - 12);
    F_u1 << P_cons;
    F_t << F_u1;
    F1 << F_u1*F_t;
    J1 << alpha*E1 + (1-alpha)*F1;
    L1 << a*J1 + b*u_v1*u_v1 + c*u_phi1*u_phi1;

	ocp.minimizeLagrangeTerm(L1);
    
    f << dot(x1) == v_a1*cos(psi1) + v_w*cos(psi_w);
    f << dot(y1) == v_a1*sin(psi1) + v_w*sin(psi_w);
    f << dot(psi1) == g*tan(phi1) / sqrt( (v_a1*cos(psi1) + v_w*cos(psi_w))*(v_a1*cos(psi1) + v_w*cos(psi_w)) + (v_a1*sin(psi1) + v_w*sin(psi_w))*(v_a1*sin(psi1) + v_w*sin(psi_w)) );
    f << dot(v_a1) == u_v1;
    f << dot(phi1) == u_phi1;

    f << dot(x0) == 2;
    f << dot(y0) == 0;
    
    //
    ocp.subjectTo(f);
    
    ocp.subjectTo(AT_START, x1 == -10);
    ocp.subjectTo(AT_START, y1 == -10);
    ocp.subjectTo(AT_START, psi1 == -1.0);
    ocp.subjectTo(AT_START, phi1 == 0);
    ocp.subjectTo(AT_START, x0 == 100);
    ocp.subjectTo(AT_START, y0 == 100);
    ocp.subjectTo(AT_START, v_a1 == 12);

    ocp.subjectTo(v_a_min <= v_a1 <= v_a_max);
    ocp.subjectTo(phi_min <= phi1 <= phi_max);
    ocp.subjectTo(dphi_min <= u_phi1 <= dphi_max);
    ocp.subjectTo(-0.2 <= u_v1 <= 0.2);

    OptimizationAlgorithm algo(ocp);

// algo.set(PARETO_FRONT_GENERATION,PFG_NORMAL_BOUNDARY_INTERSECTION');
// algo.set(PARETO_FRONT_DISCRETIZATION,41);
    algo.set(HESSIAN_APPROXIMATION, EXACT_HESSIAN);
    //algo.set(QP_SOLVER,                   QP_NONE);
// algo.set( HOTSTART_QP,                 NO             	);
// algo.set( LEVENBERG_MARQUARDT, 		 1e-10				);
// algo.set( CG_HARDCODE_CONSTRAINT_VALUES,YES);
// algo.set( SPARSE_QP_SOLUTION,          FULL_CONDENSING   );
    //algo.set(INTEGRATOR_TYPE,             INT_IRK_GL2);
    //algo.set(NUM_INTEGRATOR_STEPS,        50);
    algo.set(MAX_NUM_ITERATIONS, 300);
    //algo.set(INTEGRATOR_TYPE      , INT_RK78);
   //algo.set( INTEGRATOR_TOLERANCE , 1e-5);
   //algo.set( DISCRETIZATION_TYPE  , SINGLE_SHOOTING );
   algo.set( KKT_TOLERANCE        , 5e-2);
//	  algo.set(KKT_TOLERANCE,1e-12);

    algo.solve();
    algo.getDifferentialStates("states.txt");
    algo.getControls("controls.txt");


	return 0;