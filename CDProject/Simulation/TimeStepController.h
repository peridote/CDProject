#ifndef __TIMESTEPCONTROLLER_h__
#define __TIMESTEPCONTROLLER_h__

#include "Common/Common.h"
#include "TimeStep.h"
#include "SimulationModel.h"
#include "CollisionDetection.h"
#include "Simulation.h"

namespace PBD
{
	class TimeStepController : public TimeStep
	{
	public: 		
// 		static int SOLVER_ITERATIONS;
// 		static int SOLVER_ITERATIONS_V;
		static int NUM_SUB_STEPS;
		static int MAX_ITERATIONS;
		static int MAX_ITERATIONS_V;
		static int VELOCITY_UPDATE_METHOD;

		static int ENUM_VUPDATE_FIRST_ORDER;
		static int ENUM_VUPDATE_SECOND_ORDER;

	protected:
		int m_velocityUpdateMethod;
		unsigned int m_iterations;
		unsigned int m_iterationsV;
		unsigned int m_subSteps;
		unsigned int m_maxIterations;
		unsigned int m_maxIterationsV;
		Simulation* m_simulation;

		virtual void initParameters();
		
		void positionConstraintProjection(SimulationModel &model);
		void velocityConstraintProjection(SimulationModel &model);
		void step_liu(SimulationModel& model);


	public:
		TimeStepController(Simulation* simulation);
		virtual ~TimeStepController(void);

		virtual void step(SimulationModel &model);
		virtual void steps();
		virtual void reset();
	};
}

#endif
