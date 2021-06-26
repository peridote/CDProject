#ifndef _TIMEMANAGER_H
#define _TIMEMANAGER_H

#include "Common/Common.h"

namespace PBD
{
	class TimeManager
	{
	public:
		Real time;
		static TimeManager *current;
		static TimeManager* current2;
		Real h;

	public:
		TimeManager ();
		~TimeManager ();

		// Singleton
		static TimeManager* getCurrent ();
		static TimeManager* getCurrent2();
		static void setCurrent (TimeManager* tm);
		static bool hasCurrent();

		Real getTime();
		void setTime(Real t);
		Real getTimeStepSize();
		void setTimeStepSize(Real tss);
	};
}

#endif
