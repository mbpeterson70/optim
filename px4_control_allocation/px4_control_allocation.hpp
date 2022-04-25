#include "px4_matrix/matrix/math.hpp"
#include "px4_optim/bfgs.hpp"

#define CA_ROTOR_RIGHT 0
#define CA_ROTOR_LEFT 1
#define CA_ROTOR_REAR 2
#define CA_SERVO_RIGHT 3
#define CA_SERVO_LEFT 4
#define CA_ELEVON_RIGHT 5
#define CA_ELEVON_LEFT 6

// VTOL Vehicle parameters
#define VTOL_RHO 1.2682
#define VTOL_NCELLS 3
#define VTOL_V_MAX (3.7 * VTOL_NCELLS)
#define VTOL_S_WING 0.2589
#define VTOL_B 1.4224
#define VTOL_C 0.3305
#define VTOL_C_ELL_DELTA_A 0.018
#define VTOL_C_M_DELTA_E (-0.05)
// #define VTOL_K (matrix::Matrix<float, NUM_AXIS, NUM_AXIS>
// Identity(VTOL_NUM_AXIS, VTOL_NUM_AXIS))
#define VTOL_SERVO_MAX (115.0 * M_PI / 180)

// Rotor positions
#define VTOL_Q0_X 0.12
#define VTOL_Q0_Y 0.2
#define VTOL_Q0_Z 0.0
#define VTOL_Q1_X 0.12
#define VTOL_Q1_Y (-0.2)
#define VTOL_Q1_Z 0.0
#define VTOL_Q2_X (-0.24)
#define VTOL_Q2_Y 0.0
#define VTOL_Q2_Z 0.0

// Rear rotor paramaters
#define C_Q0_REAR 0.0216
#define C_Q1_REAR 0.0292
#define C_Q2_REAR (-0.0368)
#define C_T0_REAR 0.2097
#define C_T1_REAR 0.0505
#define C_T2_REAR (-0.1921)
#define D_PROP_REAR (5.5*(0.0254))
#define KV_REAR 1550.0
#define KQ_REAR ((1. / KV_REAR) * 60. / (2. * M_PI))
#define R_MOTOR_REAR 0.4
#define I0_REAR 0.6

// Front rotor parameters
#define C_Q0_FRONT 0.0088
#define C_Q1_FRONT 0.0129
#define C_Q2_FRONT (-0.0216)
#define C_T0_FRONT 0.1167
#define C_T1_FRONT 0.0144
#define C_T2_FRONT (-0.1480)
#define D_PROP_FRONT 7.0*(0.0254)
#define KV_FRONT 1450.0
#define KQ_FRONT (1. / KV_FRONT) * 60. / (2. * M_PI)
#define R_MOTOR_FRONT 0.3
#define I0_FRONT 0.83

using std::pow;

struct OptFunData {
    matrix::Vector<float, VTOL_NUM_AXIS> thrustTorqueDesired;
    matrix::Vector3f vBody;
    float Gamma;  
    float VaRear;
};

/**
 * Returns the thrust and torque of a given rotor.
 * thrustTorqueDers is an array that contains at array indices:
 * - 1: thrust
 * - 2: torque
 * - 3: thrust derivative
 * - 4: torque derivative
 * Va: airspeed coming into rotor
 * rotorNum: indicates which rotor's thrust and torque should be found
 */
void rotorThrustTorque(float *thrustTorqueDers, float delta, float Va, uint8_t rotorNum) {
    float C_Q0;
    float C_Q1;
    float C_Q2;
    float C_T0;
    float C_T1;
    float C_T2;
    float D_prop;
    float KQ;
    float R_motor;
    float i0;
    
    // Set parameters for correct rotor
    if (rotorNum == CA_ROTOR_REAR) {
        C_Q0 = C_Q0_REAR;
        C_Q1 = C_Q1_REAR;
        C_Q2 = C_Q2_REAR;
        C_T0 = C_T0_REAR;
        C_T1 = C_T1_REAR;
        C_T2 = C_T2_REAR;
        D_prop = D_PROP_REAR;
        KQ = KQ_REAR;
        R_motor = R_MOTOR_REAR;
        i0 = I0_REAR;
    } else {
        C_Q0 = C_Q0_FRONT;
        C_Q1 = C_Q1_FRONT;
        C_Q2 = C_Q2_FRONT;
        C_T0 = C_T0_FRONT;
        C_T1 = C_T1_FRONT;
        C_T2 = C_T2_FRONT;
        D_prop = D_PROP_FRONT;
        KQ = KQ_FRONT;
        R_motor = R_MOTOR_FRONT;
        i0 = I0_FRONT;
    }

    // map delta_t throttle command(0 to 1) into motor input voltage
    float V_in = VTOL_V_MAX * delta;
    float V_in_der = VTOL_V_MAX;
    // Quadratic formula to solve for motor speed
    float a = C_Q0 * VTOL_RHO * pow(D_prop, 5)
        / (pow((2.f * M_PI), 2));
    float b = (C_Q1 * VTOL_RHO * pow(D_prop, 4)
        / (2. * M_PI)) * Va + pow(KQ, 2) / R_motor;
    float c = C_Q2 * VTOL_RHO * pow(D_prop, 3)
        * pow(Va, 2) - (KQ / R_motor) * V_in + KQ * i0;
    float c_der = (KQ / R_motor) * V_in_der;
    // Consider only positive root
    float Omega_op = (-b + std::sqrt(pow(b, 2) - 4*a*c)) / (2.*a);
    float Omega_op_der = c_der / std::sqrt(pow(b, 2) - 4*a*c);
    // compute advance ratio
    float J_op = 2.f * M_PI * Va / (Omega_op * D_prop);
    float J_op_der = -2.f * M_PI * Va * Omega_op_der / (pow(Omega_op, 2) * D_prop);
    // compute non-dimensionalized coefficients of thrust and torque
    float C_T = C_T2 * pow(J_op, 2) + C_T1 * J_op + C_T0;
    float C_Q = C_Q2 * pow(J_op, 2) + C_Q1 * J_op + C_Q0;
    float C_T_der = 2 * C_T2 * J_op * J_op_der + C_T1 * J_op_der;
    float C_Q_der = 2 * C_Q2 * J_op * J_op_der + C_Q1 * J_op_der;
    // add thrust and torque due to propeller
    float n = Omega_op / (2 * M_PI);
    thrustTorqueDers[0] = VTOL_RHO * pow(n, 2) * pow(D_prop, 4) * C_T; // thrust value
    thrustTorqueDers[1] = VTOL_RHO * pow(n, 2) * pow(D_prop, 5) * C_Q; // torque value
    thrustTorqueDers[2] = VTOL_RHO * Omega_op * Omega_op_der * pow(D_prop, 4) * C_T / (2 * pow(M_PI, 2)) +
            VTOL_RHO * pow(Omega_op, 2) * pow(D_prop, 4) * C_T_der / pow(2 * M_PI, 2); // thrust derivative
    thrustTorqueDers[3] = VTOL_RHO * Omega_op * Omega_op_der * pow(D_prop, 5) * C_Q / (2 * pow(M_PI, 2)) +
            VTOL_RHO * pow(Omega_op, 2) * pow(D_prop, 5) * C_Q_der / pow(2 * M_PI, 2); // torque derivative
    // Negates torque moment for left rotor to account for its CW direction of rotation
    if (rotorNum == CA_ROTOR_LEFT) {
        thrustTorqueDers[1] *= -1;
        thrustTorqueDers[3] *= -1;
    }
}

/**
 * Calculates the thrust and torque achieved by a given set of rotor/servo/elevon setpoints
 * thrustTorqueAchieved: Will contain the thrust and torque achieved
 * thrust: thrust array (length 3) of each of the rotors
 * torque: torque array (length 3) of each of the rotors
 * x: Actuator setpoints
 * Gamma: Term used to determine elevon effectiveness
 */
void calcThrustTorqueAchieved(matrix::Vector<float, VTOL_NUM_AXIS> *thrustTorqueAchieved, 
    float *thrust, float *torque, matrix::Vector<float, VTOL_NUM_ACTUATORS> x, float Gamma) {
    float x_0 = std::cos(x(CA_SERVO_RIGHT));
    float z_0 = std::sin(x(CA_SERVO_RIGHT));
    float x_1 = std::cos(x(CA_SERVO_LEFT));
    float z_1 = std::sin(x(CA_SERVO_LEFT));

    float T_x =  (thrust[0] * x_0) + (thrust[1] * x_1);
    float T_z = -(thrust[0] * z_0) - (thrust[1] * z_1) - thrust[2];
    float Tau_x =   - (torque[0] * x_0) - (thrust[0] * VTOL_Q0_Y * z_0)
                    - (torque[1] * x_1) - (thrust[1] * VTOL_Q1_Y * z_1)
                    - (thrust[2] * VTOL_Q2_Y)
                    - Gamma * VTOL_B * VTOL_C_ELL_DELTA_A * x(CA_ELEVON_RIGHT)
                    + Gamma * VTOL_B * VTOL_C_ELL_DELTA_A * x(CA_ELEVON_LEFT);
    float Tau_y =   thrust[0] * (VTOL_Q0_Z * x_0 + VTOL_Q0_X * z_0) +
                    thrust[1] * (VTOL_Q1_Z * x_1 + VTOL_Q1_X * z_1) +
                    thrust[2] * VTOL_Q2_X +
                    Gamma * VTOL_C * VTOL_C_M_DELTA_E * x(CA_ELEVON_RIGHT) +
                    Gamma * VTOL_C * VTOL_C_M_DELTA_E * x(CA_ELEVON_LEFT);
    float Tau_z =   - (thrust[0] * VTOL_Q0_Y * x_0) + (torque[0] * z_0)
                    - (thrust[1] * VTOL_Q1_Y * x_1) + (torque[1] * z_1)
                    + torque[2];

    (*thrustTorqueAchieved)(0) = T_x;
    (*thrustTorqueAchieved)(1) = T_z;
    (*thrustTorqueAchieved)(2) = Tau_x;
    (*thrustTorqueAchieved)(3) = Tau_y;
    (*thrustTorqueAchieved)(4) = Tau_z;
}

/**
 * Calculates the thrust and torque achieved by a given set of rotor/servo/elevon setpoints
 * thrustTorqueAchieved: Will contain the thrust and torque achieved
 * thrust: thrust array (length 3) of each of the rotors
 * torque: torque array (length 3) of each of the rotors
 * thrustDer: thrust derivative array (length 3) of each of the rotors
 * torqueDer: torque derivative array (length 3) of each of the rotors
 * x: Actuator setpoints
 * Gamma: Term used to determine elevon effectiveness
 */
void calcThrustTorqueAchievedDer(
    matrix::Matrix<float, VTOL_NUM_ACTUATORS, VTOL_NUM_AXIS> *thrustTorqueAchievedDer, 
    float *thrust, float *torque, float *thrustDer, float *torqueDer, 
    matrix::Vector<float, VTOL_NUM_ACTUATORS> x, float Gamma) {
    float x_0 = std::cos(x(CA_SERVO_RIGHT));
    float z_0 = std::sin(x(CA_SERVO_RIGHT));
    float x_1 = std::cos(x(CA_SERVO_LEFT));
    float z_1 = std::sin(x(CA_SERVO_LEFT));

    matrix::Vector<float, VTOL_NUM_ACTUATORS> T_x_der;
    matrix::Vector<float, VTOL_NUM_ACTUATORS> T_z_der;
    matrix::Vector<float, VTOL_NUM_ACTUATORS> Tau_x_der;
    matrix::Vector<float, VTOL_NUM_ACTUATORS> Tau_y_der;
    matrix::Vector<float, VTOL_NUM_ACTUATORS> Tau_z_der;

    T_x_der(0) = thrustDer[0] * x_0;
    T_x_der(1) = thrustDer[1] * x_1;
    T_x_der(2) = 0;
    T_x_der(3) = -thrust[0] * z_0;
    T_x_der(4) = -thrust[1] * z_1;
    T_x_der(5) = 0;
    T_x_der(6) = 0;

    T_z_der(0) = -thrustDer[0] * z_0;
    T_z_der(1) = -thrustDer[1] * z_1;
    T_z_der(2) = -thrustDer[2];
    T_z_der(3) = -thrust[0] * x_0;
    T_z_der(4) = -thrust[1] * x_1;
    T_z_der(5) = 0;
    T_z_der(6) = 0;

    Tau_x_der(0) = -(torqueDer[0] * x_0) - (thrustDer[0] * VTOL_Q0_Y * z_0);
    Tau_x_der(1) = -(torqueDer[1] * x_1) - (thrustDer[1] * VTOL_Q1_Y * z_1);
    Tau_x_der(2) = - (thrustDer[2] * VTOL_Q2_Y);
    Tau_x_der(3) = torque[0] * z_0 - thrust[0] * VTOL_Q0_Y * x_0;
    Tau_x_der(4) = torque[1] * z_1 - thrust[1] * VTOL_Q1_Y * x_1;
    Tau_x_der(5) = - Gamma * VTOL_B * VTOL_C_ELL_DELTA_A;
    Tau_x_der(6) = + Gamma * VTOL_B * VTOL_C_ELL_DELTA_A;

    Tau_y_der(0) = thrustDer[0] * (VTOL_Q0_Z * x_0 + VTOL_Q0_X * z_0);
    Tau_y_der(1) = thrustDer[1] * (VTOL_Q1_Z * x_1 + VTOL_Q1_X * z_1);
    Tau_y_der(2) = thrustDer[2] * VTOL_Q2_X;
    Tau_y_der(3) = thrust[0] * (-VTOL_Q0_Z * z_0 + VTOL_Q0_X * x_0);
    Tau_y_der(4) = thrust[1] * (-VTOL_Q1_Z * z_1 + VTOL_Q1_X * x_1);
    Tau_y_der(5) = Gamma * VTOL_C * VTOL_C_M_DELTA_E;
    Tau_y_der(6) = Gamma * VTOL_C * VTOL_C_M_DELTA_E;

    Tau_z_der(0) = -(thrustDer[0] * VTOL_Q0_Y * x_0) + (torqueDer[0] * z_0);
    Tau_z_der(1) = -(thrustDer[1] * VTOL_Q1_Y * x_1) + (torqueDer[1] * z_1);
    Tau_z_der(2) = torqueDer[2];
    Tau_z_der(3) = thrust[0] * VTOL_Q0_Y * z_0 + torque[0] * x_0;
    Tau_z_der(4) = thrust[1] * VTOL_Q1_Y * z_1 + torque[1] * x_1;
    Tau_z_der(5) = 0;
    Tau_z_der(6) = 0;

    thrustTorqueAchievedDer->setCol(0, T_x_der);
    thrustTorqueAchievedDer->setCol(1, T_z_der);
    thrustTorqueAchievedDer->setCol(2, Tau_x_der);
    thrustTorqueAchievedDer->setCol(3, Tau_y_der);
    thrustTorqueAchievedDer->setCol(4, Tau_z_der);
}

/**
 * Optimization function that is used to calculate the norm of the difference between desired thrust/torque
 * and achieved thrust/torque, as well as that function's gradient.
 * vals_inp: actuator setpoints, the input to the function
 * grad_out: The gradient of the function at the vals_inp point
 * opt_data: Extra arguments to the function (Gamma, VaRear, vBody, and thrustTorqueDesired)
 */
float nonlinearCtrlOptFun(const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp, matrix::Vector<float, VTOL_NUM_ACTUATORS> *grad_out, void *opt_data) {
    float x_0 = std::cos(vals_inp(CA_SERVO_RIGHT));
    float z_0 = std::sin(vals_inp(CA_SERVO_RIGHT));
    float x_1 = std::cos(vals_inp(CA_SERVO_LEFT));
    float z_1 = std::sin(vals_inp(CA_SERVO_LEFT));
    matrix::Matrix<float, VTOL_NUM_AXIS, VTOL_NUM_AXIS> K;
    K.setIdentity();

    OptFunData* optFunData = reinterpret_cast<OptFunData*>(opt_data);
    float VaRight = (matrix::Vector3f(x_0, 0.0, -z_0).transpose() * optFunData->vBody)(0, 0);
    float VaLeft = (matrix::Vector3f(x_1, 0.0, -z_1).transpose() * optFunData->vBody)(0, 0);

    float thrustTorqueDersRotor0[4];
    float thrustTorqueDersRotor1[4];
    float thrustTorqueDersRotor2[4];

    rotorThrustTorque(thrustTorqueDersRotor0, vals_inp(0), VaRight, 0);
    rotorThrustTorque(thrustTorqueDersRotor1, vals_inp(1), VaLeft, 1);
    rotorThrustTorque(thrustTorqueDersRotor2, vals_inp(2), optFunData->VaRear, 2);

    float thrust[3] = {thrustTorqueDersRotor0[0], thrustTorqueDersRotor1[0], thrustTorqueDersRotor2[0]};
    float torque[3] = {thrustTorqueDersRotor0[1], thrustTorqueDersRotor1[1], thrustTorqueDersRotor2[1]};
    float thrustDer[3] = {thrustTorqueDersRotor0[2], thrustTorqueDersRotor1[2], thrustTorqueDersRotor2[2]};
    float torqueDer[3] = {thrustTorqueDersRotor0[3], thrustTorqueDersRotor1[3], thrustTorqueDersRotor2[3]};

    matrix::Vector<float, VTOL_NUM_AXIS> thrustTorqueAchieved;
    calcThrustTorqueAchieved(&thrustTorqueAchieved, thrust, torque, vals_inp, optFunData->Gamma);
    matrix::Vector<float, VTOL_NUM_AXIS> thrustTorqueDiff = optFunData->thrustTorqueDesired - thrustTorqueAchieved;
    float normDiff = (0.5f * thrustTorqueDiff.transpose() * K * thrustTorqueDiff)(0, 0);
    
    matrix::Matrix<float, VTOL_NUM_ACTUATORS, VTOL_NUM_AXIS> thrustTorqueDer;
    calcThrustTorqueAchievedDer(&thrustTorqueDer, thrust, torque, thrustDer, torqueDer,
        vals_inp, optFunData->Gamma);
    matrix::Vector<float, VTOL_NUM_ACTUATORS> normDiffDer = -1.f * thrustTorqueDer * K * thrustTorqueDiff;

    if (grad_out)
        *grad_out = normDiffDer;

    return normDiff;
}

/**
 * Runs the nonlinear optimization
 * thrustTorqueDesired: desired thrust/torque
 * x0: Inital guess
 * vBody: velocity in the body frame (ideally should be in the airframe)
 * airspeed: vehicle airspeed
 * iterMax: maximum number of optimization iterations
 */
matrix::Vector<float, VTOL_NUM_ACTUATORS> computeNonlinearOpt(
    matrix::Vector<float, VTOL_NUM_AXIS> thrustTorqueDesired, 
    matrix::Vector<float, VTOL_NUM_ACTUATORS> x0, matrix::Vector3f vBody, float airspeed, size_t iterMax) {
    
    bfgs_settings_t settings;
    float lowerBoundsf[VTOL_NUM_ACTUATORS] = {0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0};
    float upperBoundsf[VTOL_NUM_ACTUATORS] = {1.0, 1.0, 1.0, VTOL_SERVO_MAX, VTOL_SERVO_MAX, 1.0, 1.0};
    matrix::Vector<float, VTOL_NUM_ACTUATORS> lowerBounds(lowerBoundsf);
    matrix::Vector<float, VTOL_NUM_ACTUATORS> upperBounds(upperBoundsf);
    settings.lower_bounds = lowerBounds;
    settings.upper_bounds = upperBounds;
    settings.iter_max = iterMax;
    
    OptFunData optFunData;
    optFunData.Gamma = .5f * VTOL_RHO * pow(airspeed, 2) * VTOL_S_WING;
    optFunData.VaRear = (matrix::Vector3f(0.0, 0.0, -1.0).transpose() * vBody)(0, 0);
    optFunData.thrustTorqueDesired = thrustTorqueDesired;
    optFunData.vBody = vBody;

    bool success = ctrlalloc_bfgs(x0, nonlinearCtrlOptFun, &optFunData, settings);

    return x0;
}
