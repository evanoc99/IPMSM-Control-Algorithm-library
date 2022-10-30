// Author: Evan O'Connor
// This is a control loop designed for an IPMSM it takes a current reference in and produces an output to be used in a modulator
// This control algorithm uses a sychronous reference frame, with angular position given by a resolver. The algorithm implements standard current control as well as more advanced techniques which include MTPA and field weakening
// I have tried to design this code to be modular so with only a basic understanding of control theory and the hardware all you will need to do is change the parameters and include this library as part of your project and run the methods in the correct order.
// This code was produced as part of my Capstone
// https://github.com/evanoc99

#include <math.h>

#define __SQRT2				(1.414213562)
#define __SQRT3				(1.732050808)
#define __PI				3.141592653

#define CL_FREQ				20000										//Frequency at which the current loop runs
#define FW_FREQ				20000										//Frequency at which the Field weakening loop runs
#define R_Freq				20000										//Sampling Frequency of Resolver

//parameters for motor and resolver
#define VDC					(float)(300.0/2.0)							// DC bus voltage of supply
#define deltat				(float)(1.0/(CL_FREQ))						// sampling time
#define Td					(float)((3.0/4.0)*(2.0*deltat)
#define wc					(float)(((__PI/2.0)-(45.0*(__PI/180.0)))/Td)
#define Ld					(float)(0.000291) 							// d-axis inductance
#define Lq					(float)(0.000402) 							// q-axis inductance
#define Kp_d				(float)((wc*Ld)/VDC)
#define Kp_q				(float)((wc*Lq)/VDC)
#define tau_i				(float)(10.0/wc)	

#define M_pole_pairs		(float)(4.0)								// motor pole pairs
#define R_pole_pairs		(float)(2.0)								// Resolver pole pairs

#define deltat_r			(float)(1.0/(R_FREQ))						// sampling time
#define wn_r				(float)(0.06*R_FREQ)						// response time of resolver. higher constant is faster but at the cost of stability
#define epsilon_r			(float)(0.707)								// dampening of resolver
#define amplitude_r			(float)(1.0)								// resolver output amplitude
#define KP_R				(float)(2.0*epsilon_r*wn_r/amplitude_r)
#define tau_R				(float)((wn_r*wn_r)/amplitude_r)

//field weakening parameters
#define deltat_fw			(float)(1.0/(FW_FREQ))						// loop time delta
#define V_max				(float)(0.90)								//sets the threshold at which field weakening starts to take affect. In this case field weakening tries to hold the voltage at 90% of the maximum
#define	I_d_max_offset		(float)(-60)								//limits the field weakening controller to a maximum current offset from baseline
#define Kp_fw				(float)(100*(wc*Lq)/VDC)					// This parmeter will need to be changed to suite the mechanical model of the motor (disregard the equation)
#define tau_fw				(float)(100*(10.0/wc))						// This parmeter will need to be changed to suite the mechanical model of the motor (disregard the equation)



class IPMSM
{
    private:
	// IPMSM Variables
	float Is_ref = 0;
	float I_q_ref = 0;
	float I_d_ref = 0;
	float I_alpha = 0;
	float I_beta = 0;
	float I_q_meas = 0;
	float I_d_meas = 0;
	float error_d = 0;
	float error_q = 0;
	float u_d = 0;
	float u_q = 0;
	float error_d_prev = 0;
	float error_q_prev = 0;
	float u_d_prev = 0;
	float u_q_prev = 0;
	float u_alpha = 0;
	float u_beta = 0;
	float u_a = 0;
	float u_b = 0;
	float u_c = 0;
	
	// Resolver variables
	float PLL_q = 0;
	float PLL_w = 0;
	float omega_PLL = 0;
	float theta = 0;
	float PLL_q_prev = 0;
	float PLL_w_prev = 0;
	float theta_e = 0;

	//FW variables
	float FW_error = 0;
	float FW_mag = 0;
	float FW_mag_prev = 0;
	float FW_error_prev = 0;
	float I_q_max = 0;

    protected:

    public:
	void IPMSM::ResolverInitialization(const float &, const float &);			//Use at the start of the program to initialize resolver
	void IPMSM::ResolverCalc(const float &, const float &);					//Use this to calculate the angular position when a new sampling instance has been received from the resolver
	void IPMSM::ResolverCalc();													//Use this to calculate the angular position when there is no new sampling instance for the resolver.
	void IPMSM::setIref(const float & );										//Sets the current reference for the system
	void IPMSM::MTPAcalc();														//Calculates the Maximum torque per ampere for this specific motor (the polynomial in the method can be changed to suite your motor)
	void IPMSM::FWcalc();														//Calculates fieldweakening for motor
	void IPMSM::currentcontrol(const float &, const float &, const float &);
	float IPMSM::GetControlSignalA();
	float IPMSM::GetControlSignalB();
	float IPMSM::GetControlSignalC();

};

void IPMSM::ResolverInitialization(const float & R_cos, const float & R_sin) //input the scaled resolver output
{
	theta_e = atan2(R_cos,R_sin)*M_pole_pairs/R_pole_pairs;
}

void IPMSM::ResolverCalc(const float & R_cos, const float & R_sin) //input the scaled resolver output
{
	// SRF PLL control loop

	PLL_q = (float)(-sin(theta)*R_cos + cos(theta)*R_sin);
	PLL_w = (float)(PLL_w_prev + KP_R*(PLL_q - PLL_q_prev) + (KP_R/tau_R)*PLL_q*deltat_r);
	omega_PLL = (float)(PLL_w*M_pole_pairs/R_pole_pairs);
	theta = (theta + (PLL_w*deltat))%(2.0*__PI*R_pole_pairs/M_pole_pairs);

	PLL_q_prev = PLL_q;
	PLL_w_prev = PLL_w;

	theta_e = (theta*M_pole_pairs/R_pole_pairs);
}

void IPMSM::ResolverCalc() //
{
	// SRF PLL control loop
	theta = (theta + (PLL_w*deltat))%(2.0*__PI*R_pole_pairs/M_pole_pairs);
	theta_e = (theta*M_pole_pairs/R_pole_pairs);
}

//calculate resolver control output before current control output
void IPMSM::setIref(const float & I)
{
	Is_ref = I;
}

void IPMSM::MTPAcalc()
{
    //MTPA
	if(Is_ref == 0)
	{
		I_q_ref = 0;
		I_d_ref = 0;
	}
	else if (Is_ref < 0)
	{
		I_q_ref = -(0.9785*-Is_ref + 0.3205);
		I_d_ref = -0.0032*Is_ref*Is_ref - 0.0082*-Is_ref + 0.2445;
	}
	else if (Is_ref > 0)
	{
		I_q_ref = 0.9785*Is_ref + 0.3205;
		I_d_ref = -0.0032*Is_ref*Is_ref - 0.0082*Is_ref + 0.2445;
	}
}

void IPMSM::FWcalc()
{
	// field weakening
	FW_error = (float)(V_max - sqrt((u_d*u_d) + (u_q*u_q)));

	FW_mag = (float)(FW_mag_prev + Kp_fw*(FW_error - FW_error_prev) + (Kp_fw/tau_fw)*FW_error*deltat_fw); //needs to be tuned for your motor
	
	//set limits for field weakening controller
	if (FW_mag > 0)		// constraining below zero ensures we only weaken the field and do not strengthen
	{
		FW_mag = 0;
	}
	else if (FW_mag < I_d_max_offset)
	{
		FW_mag = I_d_max_offset;
	}

	FW_mag_prev = FW_mag;
	FW_error_prev = FW_error;
	
	I_d_ref = I_d_ref + FW_mag;

	// I_d current limiter  
	if (I_d_ref < -abs(Is_ref))
	{
		I_d_ref = -abs(Is_ref);
	}
	else if (I_d_ref > abs(Is_ref))
	{
		I_d_ref = abs(Is_ref);
	}

	//I_q current limiter
	I_q_max = sqrt(abs((Is_ref*Is_ref) - (I_d_ref*I_d_ref)));

	if((I_q_ref < 0) && (abs(I_q_ref) > I_q_max))
	{
		I_q_ref = -I_q_max;
	}
	else if((I_q_ref > 0) && (abs(I_q_ref) > I_q_max))
	{
		I_q_ref = I_q_max;
	}
}

// windup and limitter not implemented
void IPMSM::currentcontrol(const float & I_a, const float & I_b, const float & I_c) //input current of phase A,B,C respectively in amperes (A)
{
	I_alpha = (float)(2.0*I_a/3.0 - I_b/3.0 - I_c/3.0);
	I_beta = (float)(__SQRT3*I_b/3.0 - __SQRT3*I_c/3.0);

	I_d_meas = (float)(cos(theta_e)*I_alpha + sin(theta_e)*I_beta);
	I_q_meas = (float)(-sin(theta_e)*I_alpha + cos(theta_e)*I_beta);

	error_d = (float)(I_d_ref -I_d_meas); 
	error_q = (float)(I_q_ref - I_q_meas);

	u_d = (float)(u_d_prev + Kp_d*(error_d - error_d_prev) + (Kp_d/tau_i)*error_d*deltat); //
	u_q = (float)(u_q_prev + Kp_q*(error_q - error_q_prev) + (Kp_q/tau_i)*error_q*deltat); // 

	error_d_prev = (float)error_d;
	error_q_prev = (float)error_q;
	u_d_prev = (float)u_d;
	u_q_prev = (float)u_q;

	u_alpha = (float)(cos(theta_e)*u_d - sin(theta_e)*u_q);
	u_beta = (float)(sin(theta_e)*u_d + cos(theta_e)*u_q);

	u_a = (float)u_alpha;
	u_b = (float)((-1.0/2.0)*u_alpha + (__SQRT3/2.0)*u_beta);
	u_c = -u_a - u_b;
}

//Control signal is used for input into modulator
float IPMSM::GetControlSignalA()
{
	return u_a;
}

float IPMSM::GetControlSignalB()
{
	return u_b;
}

float IPMSM::GetControlSignalC()
{
	return u_c;
}
