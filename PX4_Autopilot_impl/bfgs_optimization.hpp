#ifndef BFGS_OPTIMIZATION
#define BFGS_OPTIMIZATION

#include "px4_matrix/matrix/math.hpp"
#include <functional>
#include <iostream>

#define VTOL_NUM_ACTUATORS 7
#define VTOL_NUM_AXES 5
#define EPSILON_FLOAT 1.175494E-38

struct bfgs_settings_t {
    matrix::Vector<float,VTOL_NUM_ACTUATORS> lower_bounds;
    matrix::Vector<float,VTOL_NUM_ACTUATORS> upper_bounds;
    uint8_t iter_max;
    float opt_fn_value; // returned by algorithm
    size_t opt_iter;
    float opt_error_value;
};



matrix::Vector<float, VTOL_NUM_ACTUATORS> inv_transform(
    matrix::Vector<float, VTOL_NUM_ACTUATORS> vals_inp,
    matrix::Vector<float, VTOL_NUM_ACTUATORS> lower_bounds,
    matrix::Vector<float, VTOL_NUM_ACTUATORS> upper_bounds) {

        matrix::Vector<float, VTOL_NUM_ACTUATORS> vals_inv_trans;
        float eps_flt = EPSILON_FLOAT;

        for (size_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
            if (!std::isfinite(vals_inp(i))) {
                if (std::isnan(vals_inp(i))) {
                    vals_inv_trans(i) = (upper_bounds(i) - lower_bounds(i)) / 2.0f;
                }
                else if (vals_inp(i) < 0.0f) {
                    vals_inv_trans(i) = lower_bounds(i) + eps_flt;
                } else {
                    vals_inv_trans(i) = upper_bounds(i) - eps_flt;
                }
            } else {
                vals_inv_trans(i) = ( lower_bounds(i) - eps_flt + (upper_bounds(i) + eps_flt)*std::exp(vals_inp(i)) ) \
                                / ( 1.0f + std::exp(vals_inp(i)) );

                if (!std::isfinite(vals_inv_trans(i))) {
                    vals_inv_trans(i) = upper_bounds(i) - eps_flt;
                }
            }
        }
        return vals_inv_trans;
    }

inline
float
mt_sup_norm(const float a,
            const float b,
            const float c)
{
    return std::max( std::max(std::abs(a), std::abs(b)), std::abs(c) );
}

inline
size_t
mt_step(
    float& st_best,
    float& f_best,
    float& d_best,
    float& st_other,
    float& f_other,
    float& d_other,
    float& step,
    float& f_step,
    float& d_step,
    bool& bracket,
    float step_min,
    float step_max)
{
    bool bound = false;
    size_t info = 0;
    float sgnd = d_step*(d_best / std::abs(d_best));

    float theta,s,gamma, p,q,r, step_c,step_q,step_f;

    if (f_step > f_best) {
        info = 1;
        bound = true;

        theta = 3*(f_best - f_step)/(step - st_best) + d_best + d_step;
        s = mt_sup_norm(theta,d_best,d_step); // sup norm

        gamma = s*std::sqrt(std::pow(theta/s,2.f) - (d_best/s)*(d_step/s));
        if (step < st_best) {
            gamma = -gamma;
        }

        p = (gamma - d_best) + theta;
        q = ((gamma - d_best) + gamma) + d_step;
        r = p/q;

        step_c = st_best + r*(step - st_best);
        step_q = st_best + ((d_best / ((f_best - f_step)/(step - st_best) + d_best)) / 2.0f)*(step - st_best);

        if (std::abs(step_c - st_best) < std::abs(step_q - st_best)) {
            step_f = step_c;
        } else {
            step_f = step_c + (step_q - step_c)/2;
        }

        bracket = true;
    } else if (sgnd < 0.0f) {
        info = 2;
        bound = false;

        theta = 3*(f_best - f_step)/(step - st_best) + d_best + d_step;
        s = mt_sup_norm(theta,d_best,d_step); // sup norm

        gamma = s * std::sqrt(std::pow(theta/s,2.f) - (d_best/s)*(d_step/s));
        if (step > st_best) {
            gamma = -gamma;
        }

        p = (gamma - d_step) + theta;
        q = ((gamma - d_step) + gamma) + d_best;
        r = p/q;

        step_c = step + r*(st_best - step);
        step_q = step + (d_step/(d_step-d_best))*(st_best - step);

        if (std::abs(step_c-step) > std::abs(step_q-step)) {
            step_f = step_c;
        } else {
            step_f = step_q;
        }

        bracket = true;
    } else if (std::abs(d_step) < std::abs(d_best)) {
        info = 3;
        bound = true;

        theta = 3*(f_best - f_step)/(step - st_best) + d_best + d_step;
        s = mt_sup_norm(theta,d_best,d_step); // sup norm

        gamma = s*std::sqrt(std::max(0.0f,std::pow(theta/s,2.f) - (d_best/s)*(d_step/s)));
        if (step > st_best) {
            gamma = -gamma;
        }

        p = (gamma - d_step) + theta;
        q = (gamma + (d_best - d_step)) + gamma;
        r = p/q;

        if (r < 0.0f && (gamma > 1E-30f || gamma < -1E-30f)) {
            step_c = step + r*(st_best - step);
        } else if (step > st_best) {
            step_c = step_max;
        } else {
            step_c = step_min;
        }

        step_q = step + (d_step/(d_step-d_best))*(st_best - step);

        if (bracket) {
            if (std::abs(step-step_c) < std::abs(step-step_q)) {
                step_f = step_c;
            } else {
                step_f = step_q;
            }
        } else {
            if (std::abs(step-step_c) > std::abs(step-step_q)) {
                step_f = step_c;
            } else {
                step_f = step_q;
            }
        }
    } else {
        info = 4;
        bound = false;

        if (bracket) {
            theta = 3*(f_step - f_other)/(st_other - step) + d_other + d_step;
            s = mt_sup_norm(theta,d_other,d_step);

            gamma = s*std::sqrt(std::pow(theta/s,2.f) - (d_other/s)*(d_step/s));
            if (step > st_other) {
                gamma = -gamma;
            }

            p = (gamma - d_step) + theta;
            q = ((gamma - d_step) + gamma) + d_other;
            r = p/q;

            step_c = step + r*(st_other - step);
            step_f = step_c;
        }  else if (step > st_best) {
            step_f = step_max;
        } else {
            step_f = step_min;
        }
    }

    /*
     * Update the interval of uncertainty.
     */

    if (f_step > f_best) {
        st_other = step;
        f_other = f_step;
        d_other = d_step;
    } else {
        if (sgnd < 0.0f) {
            st_other = st_best;
            f_other = f_best;
            d_other = d_best;
        }

        st_best = step;
        f_best = f_step;
        d_best = d_step;
    }

    /*
     * Compute the new step and safeguard it.
     */

    step_f = std::max(step_min, std::min(step_max,step_f));
    step = step_f;

    if (bracket && bound) {
        if (st_other > st_best) {
            step = std::min(st_best + 0.66f*(st_other - st_best), step);
        } else {
            step = std::max(st_best + 0.66f*(st_other - st_best), step);
        }
    }

    //

    return info;
}

float line_search_mt(
    float step,
    matrix::Vector<float, VTOL_NUM_ACTUATORS>& x,
    matrix::Vector<float, VTOL_NUM_ACTUATORS>& grad,
    const matrix::Vector<float, VTOL_NUM_ACTUATORS>& direc,
    const float* wolfe_cons_1_inp,
    const float* wolfe_cons_2_inp,
    std::function<float (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp, matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out, void* opt_data)> opt_objfn,
    void* opt_data)
{
    const size_t iter_max = 100;

    const float step_min = 0.0f;
    const float step_max = 10.0f;
    const float xtol = 1E-04;

    // Wolfe parameters
    const float wolfe_cons_1 = (wolfe_cons_1_inp) ? *wolfe_cons_1_inp : 1E-03f; // tolerance on the Armijo sufficient decrease condition; sometimes labelled 'mu'.
    const float wolfe_cons_2 = (wolfe_cons_2_inp) ? *wolfe_cons_2_inp : 0.90f;  // tolerance on the curvature condition; sometimes labelled 'eta'.

    //

    size_t info = 0, infoc = 1;
    const float extrap_delta = 4; // 'delta' on page 20

    matrix::Vector<float, VTOL_NUM_ACTUATORS> x_0 = x;

    float f_step = opt_objfn(x,&grad,opt_data); // q(0)

    float dgrad_init = grad.dot(direc);

    if (dgrad_init >= 0.0f) {
        return step;
    }

    float dgrad = dgrad_init;

    //

    size_t iter = 0;

    bool bracket = false, stage_1 = true;

    float f_init = f_step, dgrad_test = wolfe_cons_1*dgrad_init;
    float width = step_max - step_min, width_old = 2*width;

    float st_best = 0.0f, f_best = f_init, dgrad_best = dgrad_init;
    float st_other = 0.0f, f_other = f_init, dgrad_other = dgrad_init;

    std::cout << "x2: ";
    for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
        std::cout << x(i) << " ";
    }
    std::cout << std::endl;

    while (1) {
        std::cout << "x3: ";
        for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
            std::cout << x(i) << " ";
        }
        std::cout << std::endl;

        ++iter;

        float st_min, st_max;

        if (bracket) {
            st_min = std::min(st_best,st_other);
            st_max = std::max(st_best,st_other);
        } else {
            st_min = st_best;
            st_max = step + extrap_delta*(step - st_best);
        }

        step = std::min(std::max(step,step_min),step_max);

        if ( (bracket && (step <= st_min || step >= st_max)) \
                || iter >= iter_max-1 || infoc == 0 || (bracket && st_max-st_min <= xtol*st_max) ) {
            step = st_best;
        }

        //
        std::cout << "x4: ";
        for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
            std::cout << x(i) << " ";
        }
        std::cout << std::endl;
        std::cout << "direc: ";
        for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
            std::cout << direc(i) << " ";
        }
        std::cout << std::endl;
        x = x_0 + step * direc;
        std::cout << "x4.25: ";
        for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
            std::cout << x(i) << " ";
        }
        std::cout << std::endl;
        f_step = opt_objfn(x,&grad,opt_data);
        std::cout << "x4.5: ";
        for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
            std::cout << x(i) << " ";
        }
        std::cout << std::endl;

        dgrad = grad.dot(direc);
        float armijo_check_val = f_init + step*dgrad_test;

        // check stop conditions

        if ((bracket && (step <= st_min || step >= st_max)) || infoc == 0) {
            info = 6;
        }
        if ((step <= step_max + 1E-30f && step >= step_max - 1E-30f) && f_step <= armijo_check_val && dgrad <= dgrad_test) {
            info = 5;
        }
        if ((step <= step_min + 1E-30f && step >= step_min - 1E-30f) && (f_step > armijo_check_val || dgrad >= dgrad_test)) {
            std::cout << "step1: " << step << std::endl;
            std::cout << "step_min1: " << step_min << std::endl;
            info = 4;
        }
        if (iter >= iter_max) {
            info = 3;
        }
        if (bracket && st_max-st_min <= xtol*st_max) {
            info = 2;
        }

        if (f_step <= armijo_check_val && std::abs(dgrad) <= wolfe_cons_2*(-dgrad_init))
        {   // strong Wolfe conditions
            std::cout << "f_step: " << f_step << std::endl;
            std::cout << "armijo_check_val: " << armijo_check_val << std::endl;
            std::cout << "dgrad: " << dgrad << std::endl;
            std::cout << "wolfe_cons_2: " << wolfe_cons_2 << std::endl;
            std::cout << "dgrad_init: " << dgrad_init << std::endl;
            info = 1;
        }

        if (info != 0) {
            return step;
        }

        std::cout << "info: " << info << std::endl;
        std::cout << "step: " << step << std::endl;
        std::cout << "step_min: " << step_min << std::endl;
        //

        if (stage_1 && f_step <= armijo_check_val && dgrad >= std::min(wolfe_cons_1,wolfe_cons_2)*dgrad_init) {
            stage_1 = false;
        }

        if (stage_1 && f_step <= f_best && f_step > armijo_check_val) {
            float f_mod  = f_step - step*dgrad_test;
            float f_best_mod = f_best - st_best*dgrad_test;
            float f_other_mod = f_other - st_other*dgrad_test;

            float dgrad_mod  = dgrad - dgrad_test;
            float dgrad_best_mod = dgrad_best - dgrad_test;
            float dgrad_other_mod = dgrad_other - dgrad_test;

            infoc = mt_step(st_best,f_best_mod,dgrad_best_mod,st_other,f_other_mod,dgrad_other_mod,step,f_mod,dgrad_mod,bracket,st_min,st_max);

            //

            f_best = f_best_mod + st_best*dgrad_test;
            f_other = f_other_mod + st_other*dgrad_test;

            dgrad_best = dgrad_best_mod + dgrad_test;
            dgrad_other = dgrad_other_mod + dgrad_test;
        } else {
            infoc = mt_step(st_best,f_best,dgrad_best,st_other,f_other,dgrad_other,step,f_step,dgrad,bracket,st_min,st_max);
        }

        //

        if (bracket) {
            if (std::abs(st_other - st_best) >= 0.66f*width_old) {
                step = st_best + 0.5f*(st_other - st_best);
            }

            width_old = width;
            width = std::abs(st_other - st_best);
        }

        std::cout << "x5: ";
        for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
            std::cout << x(i) << " ";
        }
        std::cout << std::endl;
    }

    //

    return step;
}


inline
void
error_reporting(matrix::Vector<float, VTOL_NUM_ACTUATORS>& out_vals,
                const matrix::Vector<float, VTOL_NUM_ACTUATORS>& x_p,
                std::function<float (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp, matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out, void* opt_data)> opt_objfn,
                void* opt_data,
                bool& success,
                const float err,
                const float err_tol,
                const size_t iter,
                const size_t iter_max,
                const int conv_failure_switch,
                bfgs_settings_t* settings_inp)
{
    success = false;

    if (conv_failure_switch == 0) {
        out_vals = x_p;

        if (err <= err_tol && iter <= iter_max) {
            success = true;
        }
    }

    if (settings_inp) {
        settings_inp->opt_fn_value     = opt_objfn(x_p,nullptr,opt_data);
        settings_inp->opt_iter         = iter;
        settings_inp->opt_error_value  = err;
    }
}

inline
void
error_reporting(matrix::Vector<float, VTOL_NUM_ACTUATORS>& out_vals,
                const matrix::Vector<float, VTOL_NUM_ACTUATORS>& x_p,
                std::function<float (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp, matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out, matrix::Matrix<float, VTOL_NUM_ACTUATORS, VTOL_NUM_ACTUATORS>* hess_out, void* opt_data)> opt_objfn,
                void* opt_data,
                bool& success,
                const float err,
                const float err_tol,
                const size_t iter,
                const size_t iter_max,
                const int conv_failure_switch,
                bfgs_settings_t* settings_inp)
{
    std::function<float (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp, matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out, void* opt_data)> lam_objfn \
    = [opt_objfn] (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp2, matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out2, void* opt_data2)
    -> float
    {
        return opt_objfn(vals_inp2,grad_out2,nullptr,opt_data2);
    };

    //

    error_reporting(out_vals, x_p, lam_objfn, opt_data, success, err, err_tol, iter, iter_max, conv_failure_switch, settings_inp);
}


inline
matrix::Vector<float, VTOL_NUM_ACTUATORS> cwiseDivide(
    matrix::Vector<float, VTOL_NUM_ACTUATORS> v1,
    matrix::Vector<float, VTOL_NUM_ACTUATORS> v2)
{
    matrix::Vector<float, VTOL_NUM_ACTUATORS> r;
    for (size_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
        r(i) = v1(i) / v2(i);
    }
    return r;
}

inline
matrix::Matrix<float, VTOL_NUM_ACTUATORS, 1> getMatrix(
    matrix::Vector<float, VTOL_NUM_ACTUATORS> v)
{
    matrix::Matrix<float, VTOL_NUM_ACTUATORS, 1> r;
    for (size_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
        r(i, 0) = v(i);
    }
    return r;
}

bool bfgs_impl(
    matrix::Vector<float, VTOL_NUM_ACTUATORS>& init_out_vals,
    std::function<float (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp,
        matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out,
        void* opt_data)> opt_objfn,
    void* opt_data,
    bfgs_settings_t* settings_inp)
{
    // notation: 'p' stands for '+1'.

    bool success = false;

    const size_t n_vals = VTOL_NUM_ACTUATORS;

    //
    // BFGS settings

    bfgs_settings_t settings;

    if (settings_inp) {
        settings = *settings_inp;
    }

    const size_t iter_max = settings.iter_max;
    const float grad_err_tol = 1E-08;
    const float rel_sol_change_tol = 1E-14;

    const float wolfe_cons_1 = 1E-03; // line search tuning parameter
    const float wolfe_cons_2 = 0.90; // line search tuning parameter

    const matrix::Vector<float, VTOL_NUM_ACTUATORS> lower_bounds = settings.lower_bounds;
    const matrix::Vector<float,VTOL_NUM_ACTUATORS> upper_bounds = settings.upper_bounds;

    float eps_flt = EPSILON_FLOAT;

    // lambda function for box constraints

    std::function<float (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp,
        matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out,
        void* box_data)> box_objfn \
        = [opt_objfn, lower_bounds, upper_bounds] \
        (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp2,
        matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out2, void* opt_data2) \
        -> float
    {

        matrix::Vector<float, n_vals> vals_inv_trans = inv_transform(vals_inp2, lower_bounds, upper_bounds);
        float eps_flt2 = EPSILON_FLOAT;
        float ret;

        if (grad_out2) {
            matrix::Vector<float, VTOL_NUM_ACTUATORS> grad_obj = *grad_out2;

            ret = opt_objfn(vals_inv_trans,&grad_obj,opt_data2);

            matrix::SquareMatrix<float, n_vals> jacob_matrix;
            jacob_matrix.setIdentity();

            for (size_t i = 0; i < n_vals; ++i) {
                jacob_matrix(i,i) = std::exp(vals_inp2(i)) * (2*eps_flt2 + upper_bounds(i) - lower_bounds(i)) \
                                / (std::exp(2 * vals_inp2(i)) + 2*std::exp(vals_inp2(i)) + 1);
            }
            matrix::Vector<float, n_vals> jacob_vec = jacob_matrix.diag();

            for (size_t i = 0; i < n_vals; ++i)
                (*grad_out2)(i) = jacob_vec(i) * grad_obj(i);
        } else {
            ret = opt_objfn(vals_inv_trans,nullptr,opt_data2);
        }

        return ret;
    };

    // initialization

    matrix::Vector<float, VTOL_NUM_ACTUATORS> x = init_out_vals;

    // we don't check for finite values
    // if (! OPTIM_MATOPS_IS_FINITE(x) ) {
    //     printf("bfgs error: non-finite initial value(s).\n");
    //     return false;
    // }

    std::cout << "x0: ";
    for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
        std::cout << x(i) << " ";
    }
    std::cout << std::endl;

    for (size_t i = 0; i < n_vals; ++i) {
        x(i) = std::log(x(i) - lower_bounds(i) + eps_flt) - std::log(upper_bounds(i) - x(i) + eps_flt);
    }

    std::cout << "x1: ";
    for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
        std::cout << x(i) << " ";
    }
    std::cout << std::endl;


    matrix::SquareMatrix<float, n_vals> I_mat;
    I_mat.setIdentity();

    matrix::SquareMatrix<float, n_vals> W = I_mat;  // initial approx. to (inverse) Hessian
    matrix::Vector<float, n_vals> grad;             // gradient vector
    matrix::Vector<float, n_vals> d;                // direction vector
    matrix::Vector<float, n_vals> s;
    matrix::Vector<float, n_vals> y;
    d.setZero();
    s.setZero();
    y.setZero();

    box_objfn(x, &grad, opt_data);

    std::cout << "grad: ";
    for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
        std::cout << grad(i) << " ";
    }
    std::cout << std::endl;

    float grad_err = grad.norm();

    // OPTIM_BFGS_TRACE(-1, grad_err, 0.0f, x, d, grad, s, y, W);

    if (grad_err <= grad_err_tol) {
        return true;
    }

    // if ||gradient(initial values)|| > tolerance, continue

    d = - W*grad; // direction

    matrix::Vector<float, VTOL_NUM_ACTUATORS> x_p = x, grad_p = grad;

    std::cout << "x_p0: ";
    for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
        std::cout << x_p(i) << " ";
    }
    std::cout << std::endl;

    line_search_mt(1.0, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, box_objfn, opt_data);

    std::cout << "x_p1: ";
    for (uint8_t i = 0; i < VTOL_NUM_ACTUATORS; ++i) {
        std::cout << x_p(i) << " ";
    }
    std::cout << std::endl;

    s = x_p - x;
    y = grad_p - grad;

    // update approx. inverse Hessian (W)

    float W_denom_term = y.dot(s);
    matrix::SquareMatrix<float, n_vals> W_term_1;

    if (W_denom_term > 1E-10f) {
        // checking whether the curvature condition holds: y's > 0
        W_term_1 = I_mat - getMatrix(s) * (y.transpose()) / W_denom_term;

        // perform rank-1 update of inverse Hessian approximation
        W = W_term_1 * W * (W_term_1.transpose()) + getMatrix(s) * (s.transpose()) / W_denom_term;
    } else {
        W = W * 0.1f;
    }

    grad = grad_p;

    grad_err = grad_p.norm();
    float rel_sol_change = cwiseDivide(s, (x.abs() + 1.0e-08)).l1norm();

    if (grad_err <= grad_err_tol) {
        init_out_vals = x_p;
        return true;
    }

    // begin loop

    size_t iter = 0;

    while (grad_err > grad_err_tol && rel_sol_change > rel_sol_change_tol && iter < iter_max) {
        ++iter;

        //

        d = - W*grad;

        line_search_mt(1.0, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, box_objfn, opt_data);

        //

        s = x_p - x;
        y = grad_p - grad;

        W_denom_term = y.dot(s);

        if (W_denom_term > 1E-10f) {
            // checking the curvature condition y.s > 0
            W_term_1 = I_mat - getMatrix(s) * y.transpose() / W_denom_term;

            W = W_term_1 * W * W_term_1.transpose() + getMatrix(s) * s.transpose() / W_denom_term;
        }

        //

        grad_err = grad_p.norm();
        rel_sol_change = cwiseDivide(s, (x.abs() + 1.0e-08)).l1norm();

        x = x_p;
        grad = grad_p;

    }

    //

    x_p = inv_transform(x_p, lower_bounds, upper_bounds);

    int conv_failure_switch = 0;
    error_reporting(init_out_vals, x_p, opt_objfn, opt_data,
                    success, grad_err, grad_err_tol, iter, iter_max,
                    conv_failure_switch, settings_inp);

    //

    return success;
}

bool ctrlalloc_bfgs(
    matrix::Vector<float, VTOL_NUM_ACTUATORS>& init_out_vals,
    std::function<float (const matrix::Vector<float, VTOL_NUM_ACTUATORS>& vals_inp,
        matrix::Vector<float, VTOL_NUM_ACTUATORS>* grad_out,
        void* opt_data)> opt_objfn,
    void* opt_data,
    bfgs_settings_t& settings)
{
    return bfgs_impl(init_out_vals, opt_objfn, opt_data, &settings);
}

#endif