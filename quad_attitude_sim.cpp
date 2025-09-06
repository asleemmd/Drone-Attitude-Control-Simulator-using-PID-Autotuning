// quadcopter_attitude_sim_smooth_autotune.cpp
// Quadcopter attitude simulator with smooth PID autotuning
// Compile: g++ -std=c++17 -O2 quadcopter_attitude_sim_smooth_autotune.cpp -o quad_sim

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <tuple>
using namespace std;

constexpr double PI = 3.14159265358979323846;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double DEG2RAD = PI / 180.0;

// ---------------- Vector 3D ----------------
struct Vec3 {
    double x=0, y=0, z=0;
    Vec3() = default;
    Vec3(double a, double b, double c): x(a), y(b), z(c) {}
    Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x, y+o.y, z+o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x, y-o.y, z-o.z); }
    Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s, y/s, z/s); }
    static Vec3 cross(const Vec3& a, const Vec3& b) {
        return Vec3(a.y*b.z - a.z*b.y,
                    a.z*b.x - a.x*b.z,
                    a.x*b.y - a.y*b.x);
    }
};

// ---------------- Quaternion ----------------
struct Quaternion {
    double w=1, x=0, y=0, z=0;
    Quaternion() = default;
    Quaternion(double ww, double xx, double yy, double zz): w(ww), x(xx), y(yy), z(zz) {}
    Quaternion operator+(const Quaternion& o) const { return Quaternion(w+o.w, x+o.x, y+o.y, z+o.z); }
    Quaternion operator*(double s) const { return Quaternion(w*s, x*s, y*s, z*s); }
    Quaternion operator*(const Quaternion& o) const {
        return Quaternion(
            w*o.w - x*o.x - y*o.y - z*o.z,
            w*o.x + x*o.w + y*o.z - z*o.y,
            w*o.y - x*o.z + y*o.w + z*o.x,
            w*o.z + x*o.y - y*o.x + z*o.w
        );
    }
    void normalize() {
        double n = sqrt(w*w + x*x + y*y + z*z);
        if (n < 1e-12) { w=1; x=y=z=0; return; }
        w/=n; x/=n; y/=n; z/=n;
    }
    Vec3 toEuler() const {
        double sinr_cosp = 2.0 * (w*x + y*z);
        double cosr_cosp = 1.0 - 2.0 * (x*x + y*y);
        double roll = atan2(sinr_cosp, cosr_cosp);

        double sinp = 2.0 * (w*y - z*x);
        double pitch;
        if (fabs(sinp) >= 1.0) pitch = copysign(PI/2, sinp);
        else pitch = asin(sinp);

        double siny_cosp = 2.0 * (w*z + x*y);
        double cosy_cosp = 1.0 - 2.0 * (y*y + z*z);
        double yaw = atan2(siny_cosp, cosy_cosp);

        return Vec3(roll, pitch, yaw);
    }
};

// ---------------- PID Controller ----------------
struct PID {
    double kp, ki, kd;
    double integral = 0;
    double prev_error = 0;
    double out_min, out_max;
    PID() = default;
    PID(double P, double I, double D, double omin=-1e9, double omax=1e9)
        : kp(P), ki(I), kd(D), out_min(omin), out_max(omax) {}
    double update(double error, double dt) {
        integral += error * dt;
        double deriv = (dt>0) ? (error - prev_error)/dt : 0.0;
        double out = kp*error + ki*integral + kd*deriv;
        if (out > out_max) { integral -= error*dt; out = out_max; }
        if (out < out_min) { integral -= error*dt; out = out_min; }
        prev_error = error;
        return out;
    }
};

// ---------------- State ----------------
struct State {
    Quaternion q;
    Vec3 omega;
};

State derivative(const State& s, const Vec3& tau, const Vec3& Idiag) {
    Quaternion omega_q(0, s.omega.x, s.omega.y, s.omega.z);
    Quaternion qdot = s.q * omega_q * 0.5;

    Vec3 Iomega(s.omega.x*Idiag.x, s.omega.y*Idiag.y, s.omega.z*Idiag.z);
    Vec3 cross = Vec3::cross(s.omega, Iomega);
    Vec3 omegadot((tau.x - cross.x)/Idiag.x,
                  (tau.y - cross.y)/Idiag.y,
                  (tau.z - cross.z)/Idiag.z);

    State sdot;
    sdot.q = qdot;
    sdot.omega = omegadot;
    return sdot;
}

State add_scaled(const State& s, const State& k, double scale) {
    State r;
    r.q = s.q + k.q * scale;
    r.omega = s.omega + k.omega * scale;
    return r;
}

State rk4_step(const State& s, const Vec3& tau, double dt, const Vec3& Idiag) {
    State k1 = derivative(s, tau, Idiag);
    State k2 = derivative(add_scaled(s, k1, dt*0.5), tau, Idiag);
    State k3 = derivative(add_scaled(s, k2, dt*0.5), tau, Idiag);
    State k4 = derivative(add_scaled(s, k3, dt), tau, Idiag);

    Quaternion q_new = s.q + (k1.q + k2.q*2.0 + k3.q*2.0 + k4.q) * (dt/6.0);
    Vec3 omega_new   = s.omega + (k1.omega + k2.omega*2.0 + k3.omega*2.0 + k4.omega) * (dt/6.0);

    State r;
    r.q = q_new; r.q.normalize();
    r.omega = omega_new;
    return r;
}

// ---------------- Autotune functions ----------------
double autotune_PID_Ku(double dt, const Vec3& Idiag, const string& axis) {
    double Ku = 1.0;
    State s; s.q = Quaternion(1,0,0,0); s.omega = Vec3(0,0,0);
    double setpoint_rad = (axis=="roll") ? 10*DEG2RAD : (axis=="pitch") ? -6*DEG2RAD : 30*DEG2RAD;
    PID pid(Ku,0,0,-1,1);
    bool found = false;

    while(!found && Ku<50.0) {
        pid.kp = Ku; pid.integral = 0; pid.prev_error=0;
        s.q = Quaternion(1,0,0,0); s.omega = Vec3(0,0,0);
        double max_val=-1e9, min_val=1e9;
        int steps = (int)(5.0/dt);
        for(int i=0;i<=steps;i++){
            Vec3 euler = s.q.toEuler();
            double error = setpoint_rad - ((axis=="roll")?euler.x:(axis=="pitch")?euler.y:euler.z);
            Vec3 tau_vec(0,0,0);
            double tau = pid.update(error, dt);
            if(axis=="roll") tau_vec.x=tau;
            if(axis=="pitch") tau_vec.y=tau;
            if(axis=="yaw") tau_vec.z=tau;
            s = rk4_step(s, tau_vec, dt, Idiag);
            double val = (axis=="roll")?s.q.toEuler().x:(axis=="pitch")?s.q.toEuler().y:s.q.toEuler().z;
            if(val>max_val) max_val=val;
            if(val<min_val) min_val=val;
        }
        if(max_val > 1.5*setpoint_rad || min_val < -1.5*setpoint_rad) found=true;
        else Ku+=1.0;
    }
    return Ku;
}

double estimate_Pu(const Vec3& Idiag, double dt, double Ku, const string& axis) {
    State s; s.q = Quaternion(1,0,0,0); s.omega = Vec3(0,0,0);
    double setpoint_rad = (axis=="roll") ? 10*DEG2RAD : (axis=="pitch") ? -6*DEG2RAD : 30*DEG2RAD;
    PID pid(Ku,0,0,-1,1);
    int steps = (int)(5.0/dt);
    vector<double> vals;
    for(int i=0;i<=steps;i++){
        Vec3 euler = s.q.toEuler();
        double error = setpoint_rad - ((axis=="roll")?euler.x:(axis=="pitch")?euler.y:euler.z);
        Vec3 tau_vec(0,0,0);
        double tau = pid.update(error, dt);
        if(axis=="roll") tau_vec.x=tau;
        if(axis=="pitch") tau_vec.y=tau;
        if(axis=="yaw") tau_vec.z=tau;
        s = rk4_step(s, tau_vec, dt, Idiag);
        double val = (axis=="roll")?s.q.toEuler().x:(axis=="pitch")?s.q.toEuler().y:s.q.toEuler().z;
        vals.push_back(val);
    }
    vector<int> peak_idx;
    for(size_t i=1;i<vals.size()-1;i++)
        if(vals[i]>vals[i-1] && vals[i]>vals[i+1]) peak_idx.push_back(i);
    return (peak_idx.size()>=2) ? (peak_idx[1]-peak_idx[0])*dt : 0.5;
}

// Smooth gain computation
// ---------------- Smooth gain computation (tweaked for damping) ----------------
tuple<double,double,double> compute_smooth_gains(double Ku, double Pu){
    double Kp = 1*Ku;           // reduced from 0.4 for less overshoot
    double Ki = 0.2*Kp/Pu;        // slightly smaller integral
    double Kd = 2.2*Kp*Pu/8.0;    // slightly larger derivative for damping
    return make_tuple(Kp, Ki, Kd);
}


// ---------------- Main ----------------
int main() {
    double sim_time = 10.0, dt = 0.002;
    Vec3 Idiag(0.005,0.005,0.009);

    cout << "Starting smooth PID autotuning..." << endl;

    // Autotune each axis
    double Ku_roll = autotune_PID_Ku(dt, Idiag, "roll");
    double Pu_roll = estimate_Pu(Idiag, dt, Ku_roll, "roll");
    double Ku_pitch = autotune_PID_Ku(dt, Idiag, "pitch");
    double Pu_pitch = estimate_Pu(Idiag, dt, Ku_pitch, "pitch");
    double Ku_yaw = autotune_PID_Ku(dt, Idiag, "yaw");
    double Pu_yaw = estimate_Pu(Idiag, dt, Ku_yaw, "yaw");

    double Kp_roll, Ki_roll, Kd_roll;
    tie(Kp_roll, Ki_roll, Kd_roll) = compute_smooth_gains(Ku_roll, Pu_roll);
    double Kp_pitch, Ki_pitch, Kd_pitch;
    tie(Kp_pitch, Ki_pitch, Kd_pitch) = compute_smooth_gains(Ku_pitch, Pu_pitch);
    double Kp_yaw, Ki_yaw, Kd_yaw;
    tie(Kp_yaw, Ki_yaw, Kd_yaw) = compute_smooth_gains(Ku_yaw, Pu_yaw);

    cout << "Smoothed autotuned PID gains:" << endl;
    cout << "Roll  | Kp="<<Kp_roll<<" Ki="<<Ki_roll<<" Kd="<<Kd_roll<<"\n";
    cout << "Pitch | Kp="<<Kp_pitch<<" Ki="<<Ki_pitch<<" Kd="<<Kd_pitch<<"\n";
    cout << "Yaw   | Kp="<<Kp_yaw<<" Ki="<<Ki_yaw<<" Kd="<<Kd_yaw<<"\n";

    // Initialize state and PID
    State s; s.q = Quaternion(1,0,0,0); s.omega = Vec3(0,0,0);
    PID pid_roll(Kp_roll, Ki_roll, Kd_roll, -0.6,0.6);
    PID pid_pitch(Kp_pitch, Ki_pitch, Kd_pitch, -0.6,0.6);
    PID pid_yaw(Kp_yaw, Ki_yaw, Kd_yaw, -0.3,0.3);

    ofstream csv("sim_results.csv");
    csv << "time,roll_deg,pitch_deg,yaw_deg,p_deg,q_deg,r_deg,tau_x,tau_y,tau_z\n";
    csv << fixed << setprecision(6);

    cout << "Running simulation with smooth PID..." << endl;
    int steps = int(sim_time/dt);
    for(int i=0;i<=steps;i++){
        double t = i*dt;
        Vec3 euler = s.q.toEuler();
        double roll_deg = euler.x*RAD2DEG;
        double pitch_deg = euler.y*RAD2DEG;
        double yaw_deg = euler.z*RAD2DEG;

        double roll_sp_deg  = 10.0;
        double pitch_sp_deg = -6.0;
        double yaw_sp_deg   = 30.0;

        double err_roll  = roll_sp_deg*DEG2RAD - euler.x;
        double err_pitch = pitch_sp_deg*DEG2RAD - euler.y;
        double err_yaw   = yaw_sp_deg*DEG2RAD - euler.z;

        double tau_x = pid_roll.update(err_roll, dt);
        double tau_y = pid_pitch.update(err_pitch, dt);
        double tau_z = pid_yaw.update(err_yaw, dt);

        s = rk4_step(s, Vec3(tau_x, tau_y, tau_z), dt, Idiag);

        csv << t << "," << roll_deg << "," << pitch_deg << "," << yaw_deg << ","
            << s.omega.x*RAD2DEG << "," << s.omega.y*RAD2DEG << "," << s.omega.z*RAD2DEG << ","
            << tau_x << "," << tau_y << "," << tau_z << "\n";
    }
    csv.close();
    cout << "Simulation complete. Results in sim_results.csv" << endl;
    return 0;
}
