/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file ControlAllocationNonlinearOptimization.hpp
 *
 * Nonlinear optimizer based off of a quadratic thrust/torque rotor model
 *
 * @author Mason Peterson <mbpeterson70@gmail.com>
 */

#include "ControlAllocationNonlinearOptimization.hpp"
#include "bfgs_optimization.hpp"

using std::pow;


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
    float a = C_Q0 * VTOL_RHO * pow(D_prop, 5.f)
        / (pow((2.f * (float) M_PI), 2.f));
    float b = (C_Q1 * VTOL_RHO * pow(D_prop, 4.f)
        / (2.f * (float) M_PI)) * Va + pow(KQ, 2.f) / R_motor;
    float c = C_Q2 * VTOL_RHO * pow(D_prop, 3.f)
        * pow(Va, 2.f) - (KQ / R_motor) * V_in + KQ * i0;
    float c_der = (KQ / R_motor) * V_in_der;
    // Consider only positive root
    float Omega_op = (-b + std::sqrt(pow(b, 2.f) - 4.f*a*c)) / (2.f*a);
    float Omega_op_der = c_der / std::sqrt(pow(b, 2.f) - 4.f*a*c);
    // compute advance ratio
    float J_op = 2.f * (float) M_PI * Va / (Omega_op * D_prop);
    float J_op_der = -2.f * (float) M_PI * Va * Omega_op_der / (pow(Omega_op, 2.f) * D_prop);
    // compute non-dimensionalized coefficients of thrust and torque
    float C_T = C_T2 * pow(J_op, 2.f) + C_T1 * J_op + C_T0;
    float C_Q = C_Q2 * pow(J_op, 2.f) + C_Q1 * J_op + C_Q0;
    float C_T_der = 2 * C_T2 * J_op * J_op_der + C_T1 * J_op_der;
    float C_Q_der = 2 * C_Q2 * J_op * J_op_der + C_Q1 * J_op_der;
    // add thrust and torque due to propeller
    float n = Omega_op / (2 * (float) M_PI);
    thrustTorqueDers[0] = VTOL_RHO * pow(n, 2.f) * pow(D_prop, 4.f) * C_T; // thrust value
    thrustTorqueDers[1] = VTOL_RHO * pow(n, 2.f) * pow(D_prop, 5.f) * C_Q; // torque value
    thrustTorqueDers[2] = VTOL_RHO * Omega_op * Omega_op_der * pow(D_prop, 4.f) * C_T / (2 * pow((float) M_PI, 2.f)) +
            VTOL_RHO * pow(Omega_op, 2.f) * pow(D_prop, 4.f) * C_T_der / pow(2 * (float) M_PI, 2.f); // thrust derivative
    thrustTorqueDers[3] = VTOL_RHO * Omega_op * Omega_op_der * pow(D_prop, 5.f) * C_Q / (2 * pow((float) M_PI, 2.f)) +
            VTOL_RHO * pow(Omega_op, 2.f) * pow(D_prop, 5.f) * C_Q_der / pow(2 * (float) M_PI, 2.f); // torque derivative
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
void calcThrustTorqueAchieved(matrix::Vector<float, VTOL_NUM_AXES> *thrustTorqueAchieved,
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
    matrix::Matrix<float, VTOL_NUM_ACTUATORS, VTOL_NUM_AXES> *thrustTorqueAchievedDer,
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
    matrix::Matrix<float, VTOL_NUM_AXES, VTOL_NUM_AXES> K;
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

    matrix::Vector<float, VTOL_NUM_AXES> thrustTorqueAchieved;
    calcThrustTorqueAchieved(&thrustTorqueAchieved, thrust, torque, vals_inp, optFunData->Gamma);
    matrix::Vector<float, VTOL_NUM_AXES> thrustTorqueDiff = optFunData->thrustTorqueDesired - thrustTorqueAchieved;
    float normDiff = (0.5f * thrustTorqueDiff.transpose() * K * thrustTorqueDiff)(0, 0);

    matrix::Matrix<float, VTOL_NUM_ACTUATORS, VTOL_NUM_AXES> thrustTorqueDer;
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
    matrix::Vector<float, VTOL_NUM_AXES> thrustTorqueDesired,
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
    optFunData.Gamma = .5f * VTOL_RHO * pow(airspeed, 2.f) * VTOL_S_WING;
    optFunData.VaRear = (matrix::Vector3f(0.0, 0.0, -1.0).transpose() * vBody)(0, 0);
    optFunData.thrustTorqueDesired = thrustTorqueDesired;
    optFunData.vBody = vBody;

    ctrlalloc_bfgs(x0, nonlinearCtrlOptFun, &optFunData, settings);

    return x0;
}


// void
// ControlAllocationNonlinearOptimization::setEffectivenessMatrix(
// 	const matrix::Matrix<float, ControlAllocation::NUM_AXES, ControlAllocation::NUM_ACTUATORS> &effectiveness,
// 	const matrix::Vector<float, ControlAllocation::NUM_ACTUATORS> &actuator_trim, int num_actuators)
// {
// 	// does nothing
// 	// ControlAllocation::setEffectivenessMatrix(effectiveness, actuator_trim, num_actuators);
// 	// _mix_update_needed = true;
// }

void
ControlAllocationNonlinearOptimization::allocate()
{
	// PX4_INFO("actuator_sp: %.4f %.4f %.4f %.4f %.4f %.4f %.4f",
	// 	(double) _actuator_sp(0), (double) _actuator_sp(1), (double) _actuator_sp(2), (double) _actuator_sp(3),
	// 	(double) _actuator_sp(4), (double) _actuator_sp(5), (double) _actuator_sp(6));

	matrix::Vector<float, VTOL_NUM_AXES> thrustTorqueDesired;
	matrix::Vector<float, VTOL_NUM_ACTUATORS> solution;
	matrix::Vector3f vBody(_airspeed, 0.f, 0.f);

	thrustTorqueDesired(0) = _control_sp(3);
	thrustTorqueDesired(1) = _control_sp(5);
	thrustTorqueDesired(2) = _control_sp(0);
	thrustTorqueDesired(3) = _control_sp(1);
	thrustTorqueDesired(4) = _control_sp(2);

	bool set_init_actuators =
		(_actuator_sp(0) <= 1E-30f) &&
		(_actuator_sp(1) <= 1E-30f) &&
		(_actuator_sp(2) <= 1E-30f) &&
		(_actuator_sp(3) <= 1E-30f) &&
		(_actuator_sp(4) <= 1E-30f) &&
		(_actuator_sp(6) <= 1E-30f && _actuator_sp(6) >= -1E-30f) &&
		(_actuator_sp(7) <= 1E-30f && _actuator_sp(7) >= -1E-30f);
	if (!set_init_actuators) {
		solution(0) = _actuator_sp(0);
		solution(1) = _actuator_sp(1);
		solution(2) = _actuator_sp(2);
		solution(3) = _actuator_sp(3);
		solution(4) = _actuator_sp(4);
		solution(5) = _actuator_sp(6);
		solution(6) = _actuator_sp(7);
	} else {
		solution(0) = 0.5f;
		solution(1) = 0.5f;
		solution(2) = 0.5f;
		solution(3) = 1.0f;
		solution(4) = 1.0f;
		solution(5) = 0.1f;
		solution(6) = 0.1f;
	}

	solution = computeNonlinearOpt(thrustTorqueDesired, solution, vBody, _airspeed, ITER_MAX);

	// PX4_INFO("actuator_setpoint: %.4f %.4f %.4f %.4f %.4f %.4f %.4f",
	// 	(double) solution(0), (double) solution(1), (double) solution(2), (double) solution(3),
	// 	(double) solution(4), (double) solution(5), (double) solution(6));

	_actuator_sp(0) = solution(0);
	_actuator_sp(1) = solution(1);
	_actuator_sp(2) = solution(2);
	_actuator_sp(3) = solution(3);
	_actuator_sp(4) = solution(4);
	_actuator_sp(5) = 0.f;
	_actuator_sp(6) = solution(5);
	_actuator_sp(7) = solution(6);

	_control_allocated = _control_sp;

}
