#ifndef ARTICULATED_LOADER_H
#define ARTICULATED_LOADER_H

#include "ArticulatedBody.h"
#include <Environment/MeshExact.h>
#include <filesystem>
#include <tinyxml2.h>

namespace PHYSICSMOTION {
class ArticulatedLoader {
 public:
  typedef double T;
  DECL_MAT_VEC_MAP_TYPES_T
  DECL_MAP_FUNCS
  struct JointInfo {
    std::string _child;
    std::string _name;
    std::string _type;
    Mat4T _T;
    Vec3T _axis;
    Vec2T _lmt;
    T _ctrl;
    T _damp;
    //mimic
    std::string _mimic;
    T _multiplier;
    T _offset;
  };
  typedef std::unordered_map<std::string,Joint> Links;
  typedef std::unordered_map<std::string,std::vector<JointInfo>> Relations;
  static Mat3T RPY2Mat(const Vec3T& v);
  static void makeConvex(tinyxml2::XMLElement& pt);
  static void readLinks(Links& links,const tinyxml2::XMLElement& pt,const std::filesystem::path& path,bool visualMesh);
  static bool readRelations(Relations& relations,const tinyxml2::XMLElement& pt);
  static std::string findRelationByChild(const std::string& name,const Relations& relations,JointInfo& info);
  static std::string findRelationByName(const std::string& name,const Relations& relations,JointInfo& info);
  static bool buildBody(ArticulatedBody& body,const std::string& name,const Links& links,const Relations& relations);
  static ArticulatedBody readURDF(const std::string& file,bool convex,bool visualMesh);
  //built-in bodies
  static int inferDim(int rootType);
  static tinyxml2::XMLElement* addRootJoint(tinyxml2::XMLElement& pt,int rootType);
  static tinyxml2::XMLElement* createBalls(tinyxml2::XMLElement& pt,int rootType,int nr);
  static tinyxml2::XMLElement* createBird(tinyxml2::XMLElement& pt,int rootType,T bodySz,T bodyLen,T neckSz,T neckLen,T footLen1,T footLen2,T footLen3,T footRad1,T footRad2,T footRad3,bool fixFoot,bool head);
  static tinyxml2::XMLElement* createBipedal(tinyxml2::XMLElement& pt,int rootType,T bodySz,T bodyLen,T footLen1,T footLen2,T footLen3,T footRad1,T footRad2,T footRad3,bool fixFoot);
  static void addBipedalHand(tinyxml2::XMLElement& pt,T bodySz,T bodyLen,T handLen1,T handLen2,T handRad1,T handRad2);
  static tinyxml2::XMLElement* createChain(tinyxml2::XMLElement& pt,int rootType,int nr,T l,T rad,T rot,T rot0,T t,T t0,int geomDim,T yOff=0,T ratioX=1,T ratioZ=1);
  static tinyxml2::XMLElement* createGrasper(tinyxml2::XMLElement& pt,T armRad);
  static tinyxml2::XMLElement* createSpider(tinyxml2::XMLElement& pt,int rootType,T bodySz,T footLen,T footRad,T angLegLift,T rangeLeg1,T rangeLeg2,bool ball);
  static tinyxml2::XMLElement* createArm(tinyxml2::XMLElement& pt,T armLen,T armRad);
  static tinyxml2::XMLElement* createBall(tinyxml2::XMLElement& pt,int rootType,T Rad);
  static tinyxml2::XMLElement* createBox(tinyxml2::XMLElement& pt,T x,T y,T z,T rad);
  static tinyxml2::XMLElement* createDummy(tinyxml2::XMLElement& pt);
  //build-in bodies as ArticulatedBody
  static ArticulatedBody createInitJoints(int convexhulls,int sz);
  static ArticulatedBody createBird(int root,T bodySz=0.25f,T bodyLen=0.5f,T neckSz=0.1f,T neckLen=0.6f,T footLen1=0.3f,T footLen2=0.4f,T footLen3=0.45f,T footRad1=0.12f,T footRad2=0.1f,T footRad3=0.2f,bool fixFoot=false,bool head=false);
  static ArticulatedBody createBipedal(int root,T bodySz=0.25f,T bodyLen=0.5f,T footLen1=0.4f,T footLen2=0.4f,T footLen3=0.3f,T footRad1=0.12f,T footRad2=0.1f,T footRad3=0.15f,bool fixFoot=true,bool withHand=false,T handLen1=0.4f,T handLen2=0.4f,T handRad1=0.1f,T handRad2=0.1f);
  static ArticulatedBody createChain(int root,int nr,T l=0.5f,T rad=0.1f,T rot=D2R(360),T rot0=0,T t=0,T t0=0,int geomDim=3,T yOff=0,T ratioX=1,T ratioZ=1);
  static ArticulatedBody createGrasper(T armRad);
  static ArticulatedBody createSpider(int root,T bodySz=0.2f,T footLen=0.2f*sqrt(2.0f)+0.16f,T footRad=0.08f,T angLegLift=D2R(10),T rangeLeg1=D2R(45),T rangeLeg2=D2R(60),bool ball=true);
  static ArticulatedBody createArm(T armLen,T armRad);
  static ArticulatedBody createBall(int rootType,T Rad);
  static ArticulatedBody createBox(T x,T y,T z,T rad);
  static ArticulatedBody createDummy();
};
}

#endif
