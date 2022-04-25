#ifndef PX4_CONTROL_ALLOCATION_TEST
#define PX4_CONTROL_ALLOCATION_TEST

#include <iostream>
#include <chrono>
#include <ctime>  
#include <cstdlib>
#include "px4_matrix/matrix/math.hpp"
#include "PX4_Autopilot_impl/ControlAllocationNonlinearOptimization.hpp"
#include "bfgs_optimization.hpp"

// void runThrustTorqueTest() {
//     float thrustTorqueDers0[4];
//     float thrustTorqueDers1[4];
//     float thrustTorqueDers2[4];
//     rotorThrustTorque(thrustTorqueDers0, .5f, 2.0, 0);
//     rotorThrustTorque(thrustTorqueDers1, .4f, 2.5, 1);
//     rotorThrustTorque(thrustTorqueDers2, .3f, 3.0, 2);


//     std::cout << "Front rotor right: throttle = 0.5, airspeed = 2.0" << std::endl;
//     std::cout << "thrust: " << thrustTorqueDers0[0] << std::endl;
//     std::cout << "torque: " << thrustTorqueDers0[1] << std::endl;
//     std::cout << "thrust_der: " << thrustTorqueDers0[2] << std::endl;
//     std::cout << "torque_der: " << thrustTorqueDers0[3] << std::endl << std::endl;

//     std::cout << "Front rotor left: throttle = 0.4, airspeed = 2.5" << std::endl;
//     std::cout << "thrust: " << thrustTorqueDers1[0] << std::endl;
//     std::cout << "torque: " << thrustTorqueDers1[1] << std::endl;
//     std::cout << "thrust_der: " << thrustTorqueDers1[2] << std::endl;
//     std::cout << "torque_der: " << thrustTorqueDers1[3] << std::endl << std::endl;

//     std::cout << "Back rotor: throttle = 0.3, airspeed = 3.0" << std::endl;
//     std::cout << "thrust: " << thrustTorqueDers2[0] << std::endl;
//     std::cout << "torque: " << thrustTorqueDers2[1] << std::endl;
//     std::cout << "thrust_der: " << thrustTorqueDers2[2] << std::endl;
//     std::cout << "torque_der: " << thrustTorqueDers2[3] << std::endl;
// }

// // void setupThrustTorqueDer(float *thrust, float *torque, float *thrust_der, float *torque_der) {
// //     thrust = {1.2, 1.5, 1.0};
// //     torque = {0.019, 0.02, 0.018};
// //     thrust_der = {5.51, 5.21, 2.98};
// //     torque_der = {0.079, 0.069, 0.061};
// // }

// void runCalcThrustTorqueAchievedTest() {
//     matrix::Vector<float, VTOL_NUM_AXES> thrustTorqueAchieved;
//     float thrust[3] = {1.2, 1.5, 1.0};
//     float torque[3] = {0.019, 0.02, 0.018};
//     float actuatorsf[VTOL_NUM_ACTUATORS] = {0.0, 0.0, 0.0, 0.3, 0.9, 0.3, 0.6};
//     matrix::Vector<float, VTOL_NUM_ACTUATORS> actuators(actuatorsf);

//     matrix::Vector3f vBody(5.0, 0.2, -0.2);
//     float airspeed = vBody.norm();
//     float Gamma = .5 * VTOL_RHO * pow(airspeed, 2) * VTOL_S_WING;
//     std::cout << "Gamma: " << Gamma << std::endl;

//     calcThrustTorqueAchieved(&thrustTorqueAchieved, thrust, torque, actuators, Gamma);
//     std::cout << "Thrust/torque achieved:" << std::endl;
//     for (uint8_t i = 0; i < VTOL_NUM_AXES; ++i)
//         std::cout << thrustTorqueAchieved(i) << std::endl;
// }

// void runCalcThrustTorqueAchievedDerTest() {
//     matrix::Matrix<float, VTOL_NUM_ACTUATORS, VTOL_NUM_AXES> thrustTorqueAchievedDer;
//     float thrust[3] = {1.2, 1.5, 1.0};
//     float torque[3] = {0.019, 0.02, 0.018};
//     float thrust_der[3] = {5.51, 5.21, 2.98};
//     float torque_der[3] = {0.079, 0.069, 0.061};
//     float actuatorsf[VTOL_NUM_ACTUATORS] = {0.0, 0.0, 0.0, 0.3, 0.9, 0.3, 0.6};
//     matrix::Vector<float, VTOL_NUM_ACTUATORS> actuators(actuatorsf);
    
//     matrix::Vector3f vBody(5.0, 0.2, -0.2);
//     float airspeed = vBody.norm();
//     float Gamma = .5 * VTOL_RHO * pow(airspeed, 2) * VTOL_S_WING;
//     std::cout << "Gamma: " << Gamma << std::endl;

//     calcThrustTorqueAchievedDer(&thrustTorqueAchievedDer, thrust, torque, thrust_der, torque_der, actuators, Gamma);
//     std::cout << "Thrust/torque achieved derivative:" << std::endl;
//     for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
//         for (uint8_t j = 0; j < VTOL_NUM_AXES; ++j)
//             std::cout << thrustTorqueAchievedDer(i, j) << " ";
//         std::cout << std::endl;
//     }

// }

// // Compare with python function:
// // nonlinear_ctrl_optimization(
// //     [0.5, 0.52, 0.48, 0.79, 0.81, 0.2, 0.1], [1.2, -5.3, 0.12, 0.14, -0.1], 
// //     np.array([[5.0], [0.2], [-0.2]]), 0.2, 4.11735, np.eye(5))
// void runNonlinearCtrlOptFunTest() {    
//     matrix::Vector3f vBody(5.0, 0.2, -0.2);
//     float airspeed = vBody.norm();
//     float Gamma = .5 * VTOL_RHO * pow(airspeed, 2) * VTOL_S_WING;
//     std::cout << "Gamma: " << Gamma << std::endl;
    
//     float thrustTorqueDesiredf[VTOL_NUM_AXES] = {1.2, -5.3, 0.12, 0.14, -0.1};
//     matrix::Vector<float, VTOL_NUM_AXES> thrustTorqueDesired(thrustTorqueDesiredf);

//     OptFunData data;
//     data.thrustTorqueDesired = thrustTorqueDesired;
//     data.vBody = vBody;
//     data.Gamma = Gamma;
//     data.VaRear = 0.2;

//     float vals_inpf[VTOL_NUM_ACTUATORS] = {0.5, 0.52, 0.48, 0.79, 0.81, 0.2, 0.1};
//     matrix::Vector<float, VTOL_NUM_ACTUATORS> vals_inp(vals_inpf);
//     matrix::Vector<float, VTOL_NUM_ACTUATORS> grad_out;

//     nonlinearCtrlOptFun(vals_inp, &grad_out, &data);

//     std::cout << "norm: ";
//     std::cout << nonlinearCtrlOptFun(vals_inp, &grad_out, &data) << std::endl;
//     std::cout << "grad_out:" << std::endl;
//     for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i)
//         std::cout << grad_out(i) << std::endl;
// }

void printVector(matrix::Vector<float, VTOL_NUM_ACTUATORS> v) {
    for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i)
        std::cout << v(i) << std::endl;
}

void runFullTest(char *argv[]) {
    // prevSolution << 0.7556, 0.7582, 0.9235, 1.5349, 1.6065, 0.0, 0.0;
    float prevSolutionf[VTOL_NUM_ACTUATORS] = {0.5, 0.5, 0.5, VTOL_SERVO_MAX / 2, VTOL_SERVO_MAX / 2, 0.0, 0.0};
    // float prevSolutionf[VTOL_NUM_ACTUATORS] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    float thrustTorqueDesiredf[VTOL_NUM_AXES] = {0.0, -8.0, 0.0, 0.0, 0.0};
    matrix::Vector<float, VTOL_NUM_AXES> thrustTorqueDesired(thrustTorqueDesiredf);
    matrix::Vector<float, VTOL_NUM_ACTUATORS>  prevSolution(prevSolutionf);
    matrix::Vector3f vBody(5.0, 0.2, -0.2); 
    float airspeed = vBody.norm();

    if (std::string(argv[1]) == "speed") {    
        std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

        for (uint8_t i = 0; i < 250; ++i) {
            prevSolution = computeNonlinearOpt(thrustTorqueDesired, prevSolution, vBody, airspeed, atoi(argv[2]));
            thrustTorqueDesired(i % VTOL_NUM_AXES) += .02;
        }
        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        std::cout << "average time per optimization run: " << 1000 * elapsed_seconds.count() / 250 << "ms\n";
        std::cout << "speed: " << 250 / elapsed_seconds.count() << " Hz" << std::endl;
    } else if (std::string(argv[1]) == "accuracy") {
        // use: ca._compute_nonlinear_optimization([1.0, -5.0, 0.2, 0.2, 0.2], np.array([[5.0], [0.2], [-0.2]]), 5.00799)
        float thrustTorqueDesiredf[VTOL_NUM_AXES] = {1.0, -5.0, 0.2, 0.2, 0.2};
        matrix::Vector<float, VTOL_NUM_AXES> thrustTorqueDesired(thrustTorqueDesiredf);
        matrix::Vector<float, VTOL_NUM_ACTUATORS> solution = computeNonlinearOpt(thrustTorqueDesired, prevSolution, vBody, airspeed, atoi(argv[2]));
        OptFunData optFunData;
        optFunData.Gamma = .5f * VTOL_RHO * pow(airspeed, 2) * VTOL_S_WING;
        optFunData.VaRear = (matrix::Vector3f(0.0, 0.0, -1.0).transpose() * vBody)(0, 0);
        optFunData.thrustTorqueDesired = thrustTorqueDesired;
        optFunData.vBody = vBody;
        
        std::cout << "airspeed: " << airspeed << std::endl;
        std::cout << "Solution:\n";
        printVector(solution);
        std::cout << "Norm:\n";
        std::cout << nonlinearCtrlOptFun(solution, nullptr, &optFunData) << std::endl << std::endl;
    }
}

/**
 * simple test to test speed or accuracy of omptimization
 * Command line options: <test mode> <iter max>
 * test mode should either be "speed" or "accuracy"
 * iter max should be an integer
 */
int main(int argc, char *argv[]) {

    runFullTest(argv);

    return 0;
}

#endif