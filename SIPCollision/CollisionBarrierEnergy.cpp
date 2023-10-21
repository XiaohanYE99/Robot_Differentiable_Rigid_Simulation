#include "DistanceFunction.h"
#include "CollisionBarrierEnergy.h"
#include <Utils/DebugGradient.h>
#include <Utils/Utils.h>
#include <Environment/MeshExact.h>

namespace PHYSICSMOTION {
template <typename T,typename PFunc>
CollisionBarrierEnergy<T,PFunc>::CollisionBarrierEnergy
(const ArticulatedBody& body,
 const Vec& controlPoints,const ThetaTrajectory<T>& tt,
 std::vector<std::pair<EntityId<T>,EntityId<T>>>& triPairs,CCSeparatingPlanes<T>& CCPlanes,
 std::unordered_map<T,GradInfo>& gradInfo,T d0,T x0,T coef,bool JTJApprox,bool implicit)
  :TrajectorySIPEnergy<T>(body,controlPoints,tt,coef),
   _TTPairs(triPairs),_CCPlanes(CCPlanes),_gradInfo(gradInfo),
   _JTJApprox(JTJApprox),_constructGradInfo(false),
   _implicit(implicit),_output(false),_d0(d0) {
  _pTT._x0=(double)x0;
  _pCC._x0=(double)x0/2;
}
template <typename T,typename PFunc>
const std::vector<std::pair<EntityId<T>,EntityId<T>>>& CollisionBarrierEnergy<T,PFunc>::getTriPairs() const {
  return _TTPairs;
}
template <typename T,typename PFunc>
const std::vector<std::pair<EntityId<T>,EntityId<T>>>& CollisionBarrierEnergy<T,PFunc>::getFailedTT() const {
  return _failedTT;
}
template <typename T,typename PFunc>
const std::vector<std::pair<EntityId<T>,EntityId<T>>>& CollisionBarrierEnergy<T,PFunc>::getFailedCC() const {
  return _failedCC;
}
template <typename T,typename PFunc>
const CCSeparatingPlanes<T>& CollisionBarrierEnergy<T,PFunc>::getCCPlanes() const {
  return _CCPlanes;
}
template <typename T,typename PFunc>
bool CollisionBarrierEnergy<T,PFunc>::eval(EFunc* E,GFunc* G,HFunc* H) {
  bool succ=true;
  if(_constructGradInfo) {
    _gradInfo.clear();
  }
  //stage 1: TT pairs
  _failedTT.clear();
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_TTPairs.size(); i++) {
    if(!succ)
      continue;
    //make sure time stamp is the same
    T timeAvg=_TTPairs[i].first.getTimeAvg();
    if(!_TTPairs[i].second.isObstacle()) {
      T timeAvgSecond=_TTPairs[i].second.getTimeAvg();
      if(timeAvg!=timeAvgSecond) {
        _TTPairs[i].first.print();
        _TTPairs[i].second.print();
      }
      ASSERT_MSG(timeAvgSecond==timeAvg,"TimeAvg for each pair should be the same!")
    }
    //make sure gradInfo is computed
    OMP_CRITICAL_
    if(_constructGradInfo && _gradInfo.find(timeAvg)==_gradInfo.end()) {
      _gradInfo[timeAvg]=GradInfo(_body,_thetaTrajectory.getPoint(_controlPoints,timeAvg));
    }
    ASSERT_MSG(_gradInfo.find(timeAvg)!=_gradInfo.end(),"Cannot find gradInfo at timeAvg!")
    //evaluate energy in parallel
    ASSERT_MSG(!_TTPairs[i].first.isRobotConvexHull(),"CollisionBarrierEnergy._TTPairs should only contain TTPairs!")
    if(!evalTT(_TTPairs[i],E,G,H,_gradInfo[timeAvg])) {
      OMP_CRITICAL_
      _failedTT.push_back(_TTPairs[i]);
      succ=false;
    }
  }
  if(!succ)
    return false;
  //stage 2: CC pairs
  _CCPlanes.getCCPairs(_CCPairs,_planeHandles);
  _failedCC.clear();
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_CCPairs.size(); i++) {
    if(!succ)
      continue;
    //make sure time stamp is the same
    T timeAvg=_CCPairs[i].first.getTimeAvg();
    if(!_CCPairs[i].second.isObstacle()) {
      T timeAvgSecond=_CCPairs[i].second.getTimeAvg();
      if(timeAvg!=timeAvgSecond) {
        _CCPairs[i].first.print();
        _CCPairs[i].second.print();
      }
      ASSERT_MSG(timeAvgSecond==timeAvg,"TimeAvg for each pair should be the same!")
    }
    //make sure gradInfo is computed
    OMP_CRITICAL_
    if(_constructGradInfo && _gradInfo.find(timeAvg)==_gradInfo.end()) {
      _gradInfo[timeAvg]=GradInfo(_body,_thetaTrajectory.getPoint(_controlPoints,timeAvg));
    }
    ASSERT_MSG(_gradInfo.find(timeAvg)!=_gradInfo.end(),"Cannot find gradInfo at timeAvg!")
    //evaluate energy in parallel
    ASSERT_MSG(_CCPairs[i].first.isRobotConvexHull(),"CollisionBarrierEnergy._CCPairs should only contain CCPairs!")
    if(!evalCC(_CCPairs[i],*_planeHandles[i],E,G,H,_gradInfo[timeAvg])) {
      OMP_CRITICAL_
      _failedCC.push_back(_CCPairs[i]);
      succ=false;
    }
  }
  if(!succ)
    return false;
  //assemble
  Vec GTheta=Vec::Zero(_body.nrDOF());
  Mat3XT tmpDTG=Mat3XT::Zero(3,_body.nrJ()*4);
  for(auto& infoPair:_gradInfo) {
    if(G) {
      GTheta.setZero();
      infoPair.second._info.DTG(_body,mapM(tmpDTG=infoPair.second._DTG),mapV(GTheta));
    }
    if(H) {
      infoPair.second._info.toolB(_body,mapM(tmpDTG=infoPair.second._DTG),[&](int row,int col,T val) {
        SIPEnergy<T>::parallelAdd(infoPair.second._HTheta(row,col),val);
      });
    }
    _thetaTrajectory.assembleEnergy(infoPair.first,G,H,&GTheta,&(infoPair.second._HTheta));
  }
  return succ;
}
template <typename T,typename PFunc>
bool CollisionBarrierEnergy<T,PFunc>::evalCC(const std::pair<EntityId<T>,EntityId<T>>& pair,CCSeparatingPlane<T>& plane,EFunc* E,bool DG,bool DH,GradInfo& gradInfo) {
  T energyVal,timeDiff=pair.first.getTimeDiff();
  if(timeDiff==0)
    timeDiff=1; //this implies you are optimizing a static pose, we do not use timeDiff here
  ASSERT_MSG(pair.second.isRobotConvexHull() || pair.second.isObstacleConvexHull(),
             "Collision between convex hull and triangle not supported!")
  const GJKPolytope<T>& p1=pair.first.getPolytope(gradInfo);
  const GJKPolytope<T>& p2=pair.second.getPolytope(gradInfo);
  CCBarrierEnergyType cc(p1,p2,_pCC,_d0,&gradInfo,_coef*timeDiff*2,_implicit);
  //initialize
  if(!_implicit) {
    if(!plane._initialized) {
      if(!cc.initialize(&(plane._x),&_body))
        return false;
      plane._initialized=true;
    }
    cc.initialize(plane._x);
  }
  //evaluate
  cc.setOutput(_output);
  if(!cc.eval(&energyVal,&_body,(DG||DH)?&gradInfo:NULL,NULL,NULL) || !isfinite(energyVal)) {
    std::cout<<"CC Eval Error"<<std::endl;
    return false;
  }
  if(E)
    (*E)(energyVal);
  return true;
}
template <typename T,typename PFunc>
bool CollisionBarrierEnergy<T,PFunc>::evalTT(const std::pair<EntityId<T>,EntityId<T>>& pair,EFunc* E,bool DG,bool DH,GradInfo& gradInfo) const {
  int cnt=0;
  T energyVal;
  Vec12T Grad;
  Mat12T Hessian;
  bool isObs=pair.second.isObstacle();
  T timeDiff=pair.first.getTimeDiff(),mollifierCoef;
  if(timeDiff==0) {
    timeDiff=1; //this implies you are optimizing a static pose, we do not use timeDiff here
    //if user optimizes a static pose, we need discrete collision check for an entire triangle
    Vec3T va[3]= {pair.first.globalV(gradInfo,0), pair.first.globalV(gradInfo,1), pair.first.globalV(gradInfo,2)};
    Vec3T vb[3]= {pair.second.globalV(gradInfo,0),pair.second.globalV(gradInfo,1),pair.second.globalV(gradInfo,2)};
    if(triangleTriangleIntersect(va,vb))
      return false;
  }
  ASSERT(pair.first.isRobotTriangle())
  ASSERT(pair.second.isRobotTriangle() || pair.second.isObstacleTriangle())
  //triangle/vertex
  for(int i=0; i<3; ++i) {
    if(!pair.second.checkBss(i)) {
      cnt++;
      continue;
    }
    VTBarrierEnergy<T,PFunc>(pair.first.globalV(gradInfo,0),
                             pair.first.globalV(gradInfo,1),
                             pair.first.globalV(gradInfo,2),
                             pair.second.globalV(gradInfo,i),
                             _pTT,_d0,_coef*timeDiff*2).eval(&energyVal,DG?&Grad:NULL,DH?&Hessian:NULL);
    if(energyVal==0) {
      cnt++;
      continue;
    }
    if(!isfinite(energyVal))
      return false;
    if(E)
      (*E)(energyVal);
    if(DG) {
      for(int r=0; r<4; r++) {
        if(r<3)
          SIPEnergy<T>::template parallelAdd<3,-1>(gradInfo._DTG,0,pair.first._jid*4,computeDTG(Grad.template segment<3>(3+r*3),pair.first.localV(r)));
        else if(r==3 && !isObs)
          SIPEnergy<T>::template parallelAdd<3,-1>(gradInfo._DTG,0,pair.second._jid*4,computeDTG(Grad.template segment<3>(0),pair.second.localV(i)));
        if(DH) {
          for(int c=0; c<4; c++) {
            if(isObs && (r==3 || c==3))
              continue;
            gradInfo._info.toolALR(r<3?pair.first._jid:pair.second._jid,
                                   c<3?pair.first._jid:pair.second._jid,
                                   _body,Hessian.template block<3,3>(r<3?3+r*3:0,c<3?3+c*3:0),
                                   r<3?pair.first.localV(r):pair.second.localV(i),
                                   c<3?pair.first.localV(c):pair.second.localV(i),
            [&](int offr,int offc,T val) {
              SIPEnergy<T>::parallelAdd(gradInfo._HTheta(offr,offc),val);
            });
          }
        }
      }
    }
    cnt++;
  }
  //vertex/triangle
  for(int i=0; i<3; ++i) {
    if(!pair.first.checkBss(i)) {
      cnt++;
      continue;
    }
    VTBarrierEnergy<T,PFunc>(pair.second.globalV(gradInfo,0),
                             pair.second.globalV(gradInfo,1),
                             pair.second.globalV(gradInfo,2),
                             pair.first.globalV(gradInfo,i),
                             _pTT,_d0,_coef*timeDiff*2).eval(&energyVal,DG?&Grad:NULL,DH?&Hessian:NULL);
    if(energyVal==0) {
      cnt++;
      continue;
    }
    if(!isfinite(energyVal))
      return false;
    if(E)
      (*E)(energyVal);
    if(DG) {
      for(int r=0; r<4; r++) {
        if(r<3 && !isObs)
          SIPEnergy<T>::template parallelAdd<3,-1>(gradInfo._DTG,0,pair.second._jid*4,computeDTG(Grad.template segment<3>(3+r*3),pair.second.localV(r)));
        else if(r==3)
          SIPEnergy<T>::template parallelAdd<3,-1>(gradInfo._DTG,0,pair.first._jid*4,computeDTG(Grad.template segment<3>(0),pair.first.localV(i)));
        if(DH) {
          for(int c=0; c<4; c++) {
            if(isObs && (r<3 || c<3))
              continue;
            gradInfo._info.toolALR(r<3?pair.second._jid:pair.first._jid,
                                   c<3?pair.second._jid:pair.first._jid,
                                   _body,Hessian.template block<3,3>(r<3?3+r*3:0,c<3?3+c*3:0),
                                   r<3?pair.second.localV(r):pair.first.localV(i),
                                   c<3?pair.second.localV(c):pair.first.localV(i),
            [&](int offr,int offc,T val) {
              SIPEnergy<T>::parallelAdd(gradInfo._HTheta(offr,offc),val);
            });
          }
        }
      }
    }
    cnt++;
  }
  //edge/edge
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      if(!(pair.first.checkBss(3+i) && pair.second.checkBss(3+j))) {
        cnt++;
        continue;
      }
      mollifierCoef=std::max<T>((pair.first.localV(i)-pair.first.localV((i+1)%3)).squaredNorm()*
                                (pair.second.localV(j)-pair.second.localV((j+1)%3)).squaredNorm(),1)*1e-2;
      EEBarrierEnergy<T,PFunc>(pair.first.globalV(gradInfo,i),
                               pair.first.globalV(gradInfo,(i+1)%3),
                               pair.second.globalV(gradInfo,j),
                               pair.second.globalV(gradInfo,(j+1)%3),
                               _JTJApprox,_pTT,_d0,_coef*timeDiff*2,mollifierCoef).eval(&energyVal,DG?&Grad:NULL,DH?&Hessian:NULL);
      if(energyVal==0) {
        cnt++;
        continue;
      }
      if(!isfinite(energyVal))
        return false;
      if(E)
        (*E)(energyVal);
      if(DG) {
        for(int r=0; r<4; r++) {
          if(r<2)
            SIPEnergy<T>::template parallelAdd<3,-1>(gradInfo._DTG,0,pair.first._jid*4,computeDTG(Grad.template segment<3>(r*3),pair.first.localV((i+r)%3)));
          else if(r>=2 && !isObs)
            SIPEnergy<T>::template parallelAdd<3,-1>(gradInfo._DTG,0,pair.second._jid*4,computeDTG(Grad.template segment<3>(r*3),pair.second.localV((j+r-2)%3)));
          if(DH) {
            for(int c=0; c<4; c++) {
              if(isObs && (r>=2 || c>=2))
                continue;
              gradInfo._info.toolALR(r<2?pair.first._jid:pair.second._jid,
                                     c<2?pair.first._jid:pair.second._jid,
                                     _body,Hessian.template block<3,3>(r*3,c*3),
                                     r<2?pair.first.localV((i+r)%3):pair.second.localV((j+r-2)%3),
                                     c<2?pair.first.localV((i+c)%3):pair.second.localV((j+c-2)%3),
              [&](int offr,int offc,T val) {
                SIPEnergy<T>::parallelAdd(gradInfo._HTheta(offr,offc),val);
              });
            }
          }
        }
      }
      cnt++;
    }
  }
  return true;
}
template <typename T,typename PFunc>
void CollisionBarrierEnergy<T,PFunc>::randomizeCCPairs(std::shared_ptr<MeshExact> obstacle) {
  std::shared_ptr<GJKPolytope<T>> obs(new GJKPolytope<T>(obstacle));
  _TTPairs.clear();
  _CCPlanes.clear();
  //obstacle
  for(int i=0; i<_body.nrJ(); i++) {
    //robot
    EntityId<T> r;
    r._timeFrom=fmod(rand(),_thetaTrajectory.getNumSegment());
    r._timeTo=fmod(rand(),_thetaTrajectory.getNumSegment());
    sort2(r._timeFrom,r._timeTo);
    r._link=std::dynamic_pointer_cast<MeshExact>(_body.joint(i)._mesh);
    if(!r._link)
      continue;
    r._jid=i;
    r._tid=-1;
    //obstacle
    EntityId<T> o;
    o._timeFrom=r._timeFrom;
    o._timeTo=r._timeTo;
    o._obs=obs;
    if(!o._obs)
      continue;
    o._jid=-1;
    o._tid=-1;
    _CCPlanes.insertPlane(std::make_pair(r,o));
  }
  //self
  for(int i=0; i<_body.nrJ(); i++)
    for(int j=i+1; j<_body.nrJ(); j++) {
      if(_body.joint(i)._parent==j || _body.joint(j)._parent==i)
        continue;
      //robot
      EntityId<T> r;
      r._timeFrom=fmod(rand(),_thetaTrajectory.getNumSegment());
      r._timeTo=fmod(rand(),_thetaTrajectory.getNumSegment());
      sort2(r._timeFrom,r._timeTo);
      r._link=std::dynamic_pointer_cast<MeshExact>(_body.joint(i)._mesh);
      if(!r._link)
        continue;
      r._jid=i;
      r._tid=-1;
      //robot
      EntityId<T> o;
      o._timeFrom=r._timeFrom;
      o._timeTo=r._timeTo;
      o._link=std::dynamic_pointer_cast<MeshExact>(_body.joint(j)._mesh);
      if(!o._link)
        continue;
      o._jid=j;
      o._tid=-1;
      _CCPlanes.insertPlane(std::make_pair(r,o));
    }
}
template <typename T,typename PFunc>
void CollisionBarrierEnergy<T,PFunc>::randomizeTTPairs(std::shared_ptr<MeshExact> obstacle) {
  std::shared_ptr<GJKPolytope<T>> obs(new GJKPolytope<T>(obstacle));
  _TTPairs.clear();
  _CCPlanes.clear();
  //obstacle
  for(int i=0; i<_body.nrJ(); i++) {
    //robot
    EntityId<T> r;
    r._timeFrom=fmod(rand(),_thetaTrajectory.getNumSegment());
    r._timeTo=fmod(rand(),_thetaTrajectory.getNumSegment());
    sort2(r._timeFrom,r._timeTo);
    r._link=std::dynamic_pointer_cast<MeshExact>(_body.joint(i)._mesh);
    if(!r._link)
      continue;
    r._jid=i;
    r._tid=rand()%r._link->iss().size();
    //obstacle
    EntityId<T> o;
    o._timeFrom=r._timeFrom;
    o._timeTo=r._timeTo;
    o._obs=obs;
    if(!o._obs)
      continue;
    o._jid=-1;
    o._tid=rand()%obstacle->iss().size();
    _TTPairs.emplace_back(r,o);
  }
  //self
  for(int i=0; i<_body.nrJ(); i++)
    for(int j=i+1; j<_body.nrJ(); j++) {
      if(_body.joint(i)._parent==j || _body.joint(j)._parent==i)
        continue;
      //robot
      EntityId<T> r;
      r._timeFrom=fmod(rand(),_thetaTrajectory.getNumSegment());
      r._timeTo=fmod(rand(),_thetaTrajectory.getNumSegment());
      sort2(r._timeFrom,r._timeTo);
      r._link=std::dynamic_pointer_cast<MeshExact>(_body.joint(i)._mesh);
      if(!r._link)
        continue;
      r._jid=i;
      r._tid=rand()%r._link->iss().size();
      //robot
      EntityId<T> o;
      o._timeFrom=r._timeFrom;
      o._timeTo=r._timeTo;
      o._link=std::dynamic_pointer_cast<MeshExact>(_body.joint(j)._mesh);
      if(!o._link)
        continue;
      o._jid=j;
      o._tid=rand()%o._link->iss().size();
      _TTPairs.emplace_back(r,o);
    }
}
template <typename T,typename PFunc>
void CollisionBarrierEnergy<T,PFunc>::setConstructGradInfo(bool construct) {
  _constructGradInfo=construct;
}
template <typename T,typename PFunc>
void CollisionBarrierEnergy<T,PFunc>::setOutput(bool output) {
  _output=output;
}
template <typename T,typename PFunc>
bool CollisionBarrierEnergy<T,PFunc>::implicit() const {
  return _implicit;
}
template <typename T,typename PFunc>
bool CollisionBarrierEnergy<T,PFunc>::updateCCPlanes() {
  if(_implicit)
    return true;
  bool succ=true;
  _CCPlanes.getCCPairs(_CCPairs,_planeHandles);
  OMP_PARALLEL_FOR_
  for(int i=0; i<(int)_CCPairs.size(); i++) {
    if(!succ)
      continue;
    //make sure time stamp is the same
    T timeAvg=_CCPairs[i].first.getTimeAvg();
    if(!_CCPairs[i].second.isObstacle()) {
      T timeAvgSecond=_CCPairs[i].second.getTimeAvg();
      if(timeAvg!=timeAvgSecond) {
        _CCPairs[i].first.print();
        _CCPairs[i].second.print();
      }
      ASSERT_MSG(timeAvgSecond==timeAvg,"TimeAvg for each pair should be the same!")
    }
    //update CCPlane only for CCEnergy
    ASSERT_MSG(_CCPairs[i].first.isRobotConvexHull(),"CollisionBarrierEnergy._CCPairs should only contain CCPairs!")
    const std::pair<EntityId<T>,EntityId<T>>& pair=_CCPairs[i];
    //make sure gradInfo is computed
    OMP_CRITICAL_
    if(_constructGradInfo && _gradInfo.find(timeAvg)==_gradInfo.end())
      _gradInfo[timeAvg]=GradInfo(_body,_thetaTrajectory.getPoint(_controlPoints,timeAvg));
    ASSERT_MSG(_gradInfo.find(timeAvg)!=_gradInfo.end(),"Cannot find gradInfo at timeAvg!")
    //construct CCBarrierEnergy
    T timeDiff=pair.first.getTimeDiff();
    ASSERT_MSG(pair.second.isRobotConvexHull() || pair.second.isObstacleConvexHull(),
               "Collision between convex hull and triangle not supported!")
    const GJKPolytope<T>& p1=pair.first.getPolytope(_gradInfo[timeAvg]);
    const GJKPolytope<T>& p2=pair.second.getPolytope(_gradInfo[timeAvg]);
    CCBarrierEnergyType cc(p1,p2,_pCC,_d0,&_gradInfo[timeAvg],_coef*timeDiff*2,_implicit);
    //initialize
    if(!_planeHandles[i]->_initialized) {
      if(!cc.initialize(&(_planeHandles[i]->_x),&_body))
        succ=false;
      _planeHandles[i]->_initialized=true;
    }
    //update
    if(succ) {
      cc.initialize(_planeHandles[i]->_x);
      if(!cc.update(&(_planeHandles[i]->_x)))
        succ=false;
    }
  }
  return succ;
}
//instance
template class CollisionBarrierEnergy<FLOAT,Px>;
template class CollisionBarrierEnergy<FLOAT,Logx>;
template class CollisionBarrierEnergy<FLOAT,CLogx>;
template class CollisionBarrierEnergy<FLOAT,Cubicx>;
template class CollisionBarrierEnergy<FLOAT,InvQuadraticx>;
}
