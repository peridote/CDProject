#include "TimeManager.h"

using namespace PBD;

TimeManager* TimeManager::current = 0;
TimeManager* TimeManager::current2 = 0;

TimeManager::TimeManager () 
{
	time = 0;
	h = static_cast<Real>(0.005);
}

TimeManager::~TimeManager () 
{
	current = 0;
	current2 = 0;
}

TimeManager* TimeManager::getCurrent ()
{
	if (current == 0)
	{
		current = new TimeManager ();
	}
	return current;
}

TimeManager* TimeManager::getCurrent2()
{
	if (current2 == 0)
	{
		current2 = new TimeManager();
	}
	return current2;
}


void TimeManager::setCurrent (TimeManager* tm)
{
	current = tm;
}

bool TimeManager::hasCurrent()
{
	return (current != 0);
}

Real TimeManager::getTime()
{
	return time;
}

void TimeManager::setTime(Real t)
{
	time = t;
}

Real TimeManager::getTimeStepSize()
{
	return h;
}

void TimeManager::setTimeStepSize(Real tss)
{
	h = tss;
}
