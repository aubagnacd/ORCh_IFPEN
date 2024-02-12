
#include "flamemodel.h"


//---AutoIgnition---

AutoIgnition::AutoIgnition(double T_fuel, double T_oxidizer, double Pressure, double delta_t, double max_t, double Phi, string Y_fuel, string X_fuel, string Y_oxidizer, string X_oxidizer) //Constructeur
   :m_T_fuel(T_fuel), m_T_oxidizer(T_oxidizer), m_Pressure(Pressure), m_delta_t(delta_t), m_max_t(max_t), m_Phi(Phi), m_Y_fuel(Y_fuel), m_X_fuel(X_fuel), m_Y_oxidizer(Y_oxidizer), m_X_oxidizer(X_oxidizer)
{}

AutoIgnition::~AutoIgnition() //Destructeur
{}


//---Premixed---

PremixedFlames::PremixedFlames(double T_fuel, double T_oxidizer, double Pressure, double Phi, string Yf_str, string Xf_str, string Yo_str, string Xo_str, string Initial_Flame, string Final_Flame) //Constructeur
   :m_T_fuel(T_fuel), m_T_oxidizer(T_oxidizer), m_Pressure(Pressure), m_Phi(Phi), m_Yf_str(Yf_str), m_Xf_str(Xf_str), m_Yo_str(Yo_str), m_Xo_str(Xo_str), m_Initial_Flame(Initial_Flame), m_Final_Flame(Final_Flame)
{}

PremixedFlames::~PremixedFlames() //Destructeur
{}



//---MultipleInlet---

MultipleInlet::MultipleInlet(bool equil, double Temperature, double Pressure, double flowRate, string X_Species, string Y_Species, bool liquid, double DropletsDiameter, double Tau_vj, double density_liquid, double EvaporationLatentHeat) //Constructeur
   :m_equil(equil), m_Temperature(Temperature), m_Pressure(Pressure), m_flowRate(flowRate), m_X_Species(X_Species), m_Y_Species(Y_Species), m_liquid(liquid), m_DropletsDiameter(DropletsDiameter), m_Tau_vj(Tau_vj), m_density_liquid(density_liquid), m_EvaporationLatentHeat(EvaporationLatentHeat)   
{}

MultipleInlet::~MultipleInlet() //Destructeur
{}