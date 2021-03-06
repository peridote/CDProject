#ifndef __Simulation_h__
#define __Simulation_h__

#include "Common/Common.h"
#include "SimulationModel.h"
#include "ParameterObject.h"
#include "TimeStep.h"


namespace PBD
{
	enum class SimulationMethods { PBD = 0, XPBD, IBDS, TEST, NumSimulationMethods };

	/** \brief Class to manage the current simulation time and the time step size. 
	* This class is a singleton.
	*/
	class Simulation : public GenParam::ParameterObject
	{
	public:
		static int GRAVITATION;
		static int SIMULATION_METHOD;

		static int ENUM_SIMULATION_PBD;
		static int ENUM_SIMULATION_XPBD;
		static int ENUM_SIMULATION_IBDS;

	protected:
		SimulationModel *m_model;
		SimulationMethods m_simulationMethod;
		TimeStep *m_timeStep;
		TimeStep* m_timeStep2;
		Vector3r m_gravitation;
		std::function<void()> m_simulationMethodChanged;

		virtual void initParameters();
		

	private:
		static Simulation *current;

	public:
		static Simulation* m_simul1;
		static Simulation* m_simul2;

	public:
		Simulation ();
		~Simulation ();

		void init();
		void reset();

		// Singleton
		static Simulation* getCurrent ();
		// switch current between m_simul1 and m_simul2
		static void switchCurrent();
		static void setCurrent (Simulation* tm);
		static bool hasCurrent();

		SimulationModel *getModel() { return m_model; }
		void setModel(SimulationModel *model) { m_model = model; }

		int getSimulationMethod() const { return static_cast<int>(m_simulationMethod); }
		void setSimulationMethod(const int val);

		void setSimulationMethodChangedCallback(std::function<void()> const& callBackFct);

		TimeStep *getTimeStep() { return m_timeStep; }
		TimeStep* getTimeStep2() { return m_timeStep2; }
		void setTimeStep(TimeStep *ts) { m_timeStep = ts; }
	};
}

#endif
