#ifndef TIMING_H
#define TIMING_H

#include <iostream>
#include <tinyxml2.h>

namespace PHYSICSMOTION {
//timing
extern void saveTimingSetting();
extern void loadTimingSetting();
extern void enableTiming();
extern void disableTiming();
extern std::string TGETNAME(const std::string& path,const std::string& pathFrm,tinyxml2::XMLElement& pt);
extern void TFRMRESET(const std::string& pathFrm,tinyxml2::XMLElement& pt);
extern void TFRMADVANCE(const std::string& pathFrm,tinyxml2::XMLElement& pt);

extern void TBEG();
extern void TBEG(const std::string& name);
extern void TENDT(const std::string& path,tinyxml2::XMLElement& pt);
extern void TENDT(const std::string& path,const std::string& pathFrm,tinyxml2::XMLElement& pt);
extern void TEND(std::ostream& os);
extern void TEND();
extern double TENDV();
}

#endif
