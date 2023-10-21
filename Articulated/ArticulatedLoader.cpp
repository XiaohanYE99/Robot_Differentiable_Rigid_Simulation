#include "PBDArticulatedGradientInfo.h"
#include "ArticulatedLoader.h"
#include "ArticulatedUtils.h"
#include <Environment/MeshExact.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/SphericalBBoxExact.h>
#include <Utils/RotationUtils.h>
#include <Utils/Utils.h>
#include <filesystem>

namespace PHYSICSMOTION {
typename ArticulatedLoader::Mat3T ArticulatedLoader::RPY2Mat(const Vec3T& rpy) {
  Mat3T RX=expWGradV<T,Vec3T>(Vec3T::UnitX()*rpy[0]);
  Mat3T RY=expWGradV<T,Vec3T>(Vec3T::UnitY()*rpy[1]);
  Mat3T RZ=expWGradV<T,Vec3T>(Vec3T::UnitZ()*rpy[2]);
  return RZ*RY*RX;
}
void ArticulatedLoader::makeConvex(tinyxml2::XMLElement& pt) {
  for(tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(std::string(v->Name())=="link")
      put<char>(*v,"convex",1);
    else makeConvex(*v);
  }
}
void ArticulatedLoader::readLinks(Links& links,const tinyxml2::XMLElement& pt,const std::filesystem::path& path,bool visualMesh) {
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(std::string(v->Name())=="gazebo")
      continue;
    if(std::string(v->Name())=="link") {
      const tinyxml2::XMLElement& link=*v;
      if(!hasAttribute(link,"<xmlattr>.name"))
        continue;
      //name
      Joint joint;
      joint._name=get<std::string>(link,"<xmlattr>.name");
      std::vector<std::shared_ptr<ShapeExact>> geomsMerged;
      std::vector<CompositeShapeExact::Mat3X4T> transMerged;
      //mesh
      for(const tinyxml2::XMLElement* g=link.FirstChildElement(); g; g=g->NextSiblingElement()) {
        if(visualMesh && std::string(g->Name())!="visual")
          continue;
        if(!visualMesh && std::string(g->Name())!="collision")
          continue;
        const tinyxml2::XMLElement& meshPt=*g;
        //T
        Mat3X4T t=Mat3X4T::Identity();
        ROT(t)=RPY2Mat(parsePtreeDef<Vec3T>(meshPt,"origin.<xmlattr>.rpy",Vec3T::Zero()));
        CTR(t)=parsePtreeDef<Vec3T>(meshPt,"origin.<xmlattr>.xyz",Vec3T::Zero());
        //geometry
        std::shared_ptr<ShapeExact> mesh;
        if(hasAttribute(meshPt,"geometry.mesh")) {
          if(!hasAttribute(meshPt,"geometry.mesh.<xmlattr>.filename"))
            continue;
          std::string meshPath=get<std::string>(meshPt,"geometry.mesh.<xmlattr>.filename");
          if(beginsWith(meshPath,"package://"))
            meshPath=meshPath.substr(std::string("package://").length());
          if(beginsWith(meshPath,"file://"))
            meshPath=meshPath.substr(std::string("file://").length());
          std::filesystem::path p=path/meshPath;
          if(get<char>(link,"convex",(char)0))
            mesh.reset(new ConvexHullExact(p.string()));
          else if(visualMesh)
            mesh.reset(new CompositeShapeExact(p.string(),false));  //visual mesh is usually of low quality and we choose to not construct BVH for it
          else mesh.reset(new MeshExact(p.string(),true));
          std::shared_ptr<MeshExact> meshTest=std::dynamic_pointer_cast<MeshExact>(mesh);
          if(mesh->empty()) {
            mesh=NULL;
            continue;
          }
          //scale
          Vec3T scale=parsePtreeDef<Vec3T>(meshPt,"geometry.mesh.<xmlattr>.scale",Vec3T::Ones());
          if(scale!=Vec3T::Ones()) {
            std::vector<Eigen::Matrix<double,3,1>> vss;
            std::vector<Eigen::Matrix<int,3,1>> iss;
            mesh->getMesh(vss,iss);
            for(Eigen::Matrix<double,3,1>& v:vss)
              v.array()*=scale.array();
            if(get<char>(link,"convex",(char)0))
              mesh.reset(new ConvexHullExact(vss));
            else mesh.reset(new MeshExact(vss,iss,!visualMesh));
          }
        } else if(hasAttribute(meshPt,"geometry.box")) {
          Vec3T size=parsePtree<Vec3T>(meshPt,"geometry.box.<xmlattr>.size");
          mesh.reset(new BBoxExact(size[0]/2,size[1]/2,size[2]/2));
        } else if(hasAttribute(meshPt,"geometry.sphere")) {
          T r=get<T>(meshPt,"geometry.sphere.<xmlattr>.radius");
          mesh.reset(new SphericalBBoxExact(r));
        } else if(hasAttribute(meshPt,"geometry.cylinder")) {
          T r=get<T>(meshPt,"geometry.cylinder.<xmlattr>.radius");
          T l=get<T>(meshPt,"geometry.cylinder.<xmlattr>.length");
          mesh.reset(new SphericalBBoxExact(0,0,l/2,r));
        }
        if(mesh) {
          geomsMerged.push_back(mesh);
          transMerged.push_back(t.template cast<CompositeShapeExact::T>());
        }
      }
      //assemble
      ArticulatedUtils::mergeMesh(joint,geomsMerged,transMerged);
      if(joint._mesh) {
        if(hasAttribute(link,"inertial.mass.<xmlattr>.value") && !get<char>(link,"convex",(char)0)) { //if convex is required, recompute mass
          std::cout << "Using mass in URDF file!" << std::endl;
          Mat3T inertia;
          inertia(0,0)=get<T>(link,"inertial.inertia.<xmlattr>.ixx");
          inertia(1,1)=get<T>(link,"inertial.inertia.<xmlattr>.iyy");
          inertia(2,2)=get<T>(link,"inertial.inertia.<xmlattr>.izz");
          inertia(0,1)=inertia(1,0)=get<T>(link,"inertial.inertia.<xmlattr>.ixy");
          inertia(0,2)=inertia(2,0)=get<T>(link,"inertial.inertia.<xmlattr>.ixz");
          inertia(1,2)=inertia(2,1)=get<T>(link,"inertial.inertia.<xmlattr>.iyz");
          Mat3T MCCT=Mat3T::Identity()*inertia.trace()/2-inertia;
          Eigen::SelfAdjointEigenSolver<Mat3T> eig(MCCT);
          if(eig.eigenvalues().minCoeff()<0) {
            std::cout << "Found negative MCCT tensor (link=" << joint._name << "): " << eig.eigenvalues().transpose() << std::endl;
          }
          //convert
          Mat3T R=RPY2Mat(parsePtreeDef<Vec3T>(link,"inertial.origin.<xmlattr>.rpy",Vec3T::Zero()));
          Vec3T t=parsePtreeDef<Vec3T>(link,"inertial.origin.<xmlattr>.xyz",Vec3T::Zero());
          //build 6x6 matrix in body frame
          T mass=get<T>(link,"inertial.mass.<xmlattr>.value");
          joint._M=mass;
          joint._MC=mass*t;
          joint._MCCT=R*MCCT*R.transpose()+mass*t*t.transpose();
        } else {
          std::cout << "Using recomputed mass!" << std::endl;
          joint.assemble(1);
        }
      } else {
        joint._M=0;
        joint._MC.setZero();
        joint._MCCT.setZero();
      }
      links[joint._name]=joint;
    } else readLinks(links,*v,path,visualMesh);
  }
}
bool ArticulatedLoader::readRelations(Relations& relations,const tinyxml2::XMLElement& pt) {
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(std::string(v->Name())=="gazebo")
      continue;
    if(std::string(v->Name())=="joint") {
      JointInfo info;
      const tinyxml2::XMLElement& joint=*v;
      //name
      if(!hasAttribute(joint,"<xmlattr>.name"))
        continue;
      info._name=get<std::string>(joint,"<xmlattr>.name");
      //type
      if(!hasAttribute(joint,"<xmlattr>.type"))
        continue;
      if(get<std::string>(joint,"<xmlattr>.type")!="revolute" && get<std::string>(joint,"<xmlattr>.type")!="fixed" &&
          get<std::string>(joint,"<xmlattr>.type")!="continuous" && get<std::string>(joint,"<xmlattr>.type")!="prismatic")
        return false;
      info._type=get<std::string>(joint,"<xmlattr>.type");
      //parent
      if(!hasAttribute(joint,"parent.<xmlattr>.link"))
        continue;
      //child
      if(!hasAttribute(joint,"child.<xmlattr>.link"))
        continue;
      info._child=get<std::string>(joint,"child.<xmlattr>.link");
      //param
      info._T.setIdentity();
      ROT(info._T)=RPY2Mat(parsePtreeDef<Vec3T>(joint,"origin.<xmlattr>.rpy",Vec3T::Zero()));
      CTR(info._T)=parsePtreeDef<Vec3T>(joint,"origin.<xmlattr>.xyz",Vec3T::Zero());
      info._axis=parsePtreeDef<Vec3T>(joint,"axis.<xmlattr>.xyz",Vec3T::Zero());
      info._ctrl=get<T>(joint,"limit.<xmlattr>.effort",0);
      info._damp=get<T>(joint,"dynamics.<xmlattr>.damping",0);
      //lower
      if(hasAttribute(joint,"limit.<xmlattr>.lower"))
        info._lmt[0]=get<T>(joint,"limit.<xmlattr>.lower");
      else info._lmt[0]=-std::numeric_limits<T>::infinity();
      //upper
      if(hasAttribute(joint,"limit.<xmlattr>.upper"))
        info._lmt[1]=get<T>(joint,"limit.<xmlattr>.upper");
      else info._lmt[1]=std::numeric_limits<T>::infinity();
      //mimic
      if(get<std::string>(joint,"<xmlattr>.type")=="revolute") {
        if(hasAttribute(joint,"mimic.<xmlattr>.joint")) {
          info._mimic=get<std::string>(joint,"mimic.<xmlattr>.joint");
          info._multiplier=get<T>(joint,"mimic.<xmlattr>.multiplier",0);
          info._offset=get<T>(joint,"mimic.<xmlattr>.offset",0);
        }
      }
      relations[get<std::string>(joint,"parent.<xmlattr>.link")].push_back(info);
    } else if(!readRelations(relations,*v))
      return false;
  }
  return true;
}
std::string ArticulatedLoader::findRelationByChild(const std::string& name,const Relations& relations,JointInfo& info) {
  for(const std::pair<const std::string,std::vector<JointInfo>>& v:relations) {
    for(const JointInfo& c:v.second) {
      if(c._child==name) {
        info=c;
        return v.first;
      }
    }
  }
  ASSERT_MSGV(false,"Cannot find parent for link: %s!",name.c_str())
  return name;
}
std::string ArticulatedLoader::findRelationByName(const std::string& name,const Relations& relations,JointInfo& info) {
  for(const std::pair<const std::string,std::vector<JointInfo>>& v:relations) {
    for(const JointInfo& c:v.second) {
      if(c._name==name) {
        info=c;
        return v.first;
      }
    }
  }
  ASSERT_MSGV(false,"Cannot find parent with name: %s!",name.c_str())
  return name;
}
bool ArticulatedLoader::buildBody(ArticulatedBody& body,const std::string& name,const Links& links,const Relations& relations) {
  int JID=-1;
  Mat3X4T invA2Z=Mat3X4T::Identity();
  if(body._joints.size()==0) {
    //this is root
    if(links.find(name)==links.end()) {
      std::cout << "Cannot find link: " << name << std::endl;
      return false;
    }
    body._joints.push_back(links.find(name)->second);
    Joint& joint=body._joints[JID=0];
    joint._parent=-1;
    joint._depth=0;
    joint._typeJoint=Joint::FIX_JOINT;
    joint._offDOF=0;
    joint._offDDT=0;
    joint._limits.resize(3,0);
    joint._control.resize(0);
    joint._damping.resize(0);
    joint._trans.setIdentity();
  } else {
    //this is non-root, then get info
    JointInfo info;
    int nrDOF=body.nrDOF();
    int nrDDT=body.nrDDT();
    const std::string& parentName=findRelationByChild(name,relations,info);
    //build
    if(links.find(name)==links.end()) {
      std::cout << "Cannot find link: " << name << std::endl;
      return false;
    }
    body._joints.push_back(links.find(name)->second);
    Joint& joint=body._joints[JID=body.nrJ()-1];
    joint._parent=body.jointId(parentName);
    joint._depth=body.joint(joint._parent)._depth+1;
    joint._mult=joint._offset=0;
    if(info._type=="revolute" || info._type=="continuous") {
      joint._typeJoint=Joint::HINGE_JOINT;
      if(!info._mimic.empty()) {
        joint._mult=info._multiplier;
        joint._offset=info._offset;
      }
    } else if(info._type=="fixed")
      joint._typeJoint=Joint::FIX_JOINT;
    else if(info._type=="prismatic")
      joint._typeJoint=Joint::TRANS_1D;
    else {
      ASSERT_MSGV(false,"Unknown type: %s!",info._type.c_str())
    }
    joint._offDOF=nrDOF;
    joint._offDDT=nrDDT;
    if(joint._typeJoint==Joint::FIX_JOINT) {
      joint._limits.setZero(3,0);
      joint._control.resize(0);
      joint._damping.resize(0);
    } else {
      joint._limits=Vec3T(info._lmt[0],info._lmt[1],0);
      if(!std::isfinite(info._lmt[0]) || !std::isfinite(info._lmt[1]))
        joint._limits(2,0)=std::numeric_limits<T>::infinity();
      joint._control.resize(1);
      joint._control[0]=info._ctrl;
      joint._damping.resize(1);
      joint._damping[0]=info._damp;
    }
    //Axis2Z
    Mat3X4T A2Z=Mat3X4T::Identity(),trans0=info._T.block<3,4>(0,0);
    if(joint._typeJoint==Joint::HINGE_JOINT)
      ROT(A2Z)=QuatT::FromTwoVectors(Vec3T::UnitZ(),info._axis).toRotationMatrix();
    if(joint._typeJoint==Joint::TRANS_1D)
      ROT(A2Z)=QuatT::FromTwoVectors(Vec3T::UnitX(),info._axis).toRotationMatrix();
    APPLY_TRANS(joint._trans,trans0,A2Z)
    //joint
    if(joint._mesh)
      joint.transformMesh(ROT(A2Z).transpose(),Vec3T::Zero());
    INV(invA2Z,A2Z)
  }
  //recurse
  if(relations.find(name)!=relations.end()) {
    for(const JointInfo& c:relations.find(name)->second) {
      if(!buildBody(body,c._child,links,relations))
        return false;
    }
  }
  //children
  std::set<int> ids=body.children(JID,true);
  for(int c:ids) {
    Joint& joint=body._joints[c];
    Mat3X4T trans0=joint._trans;
    APPLY_TRANS(joint._trans,invA2Z,trans0)
  }
  return true;
}
ArticulatedBody ArticulatedLoader::readURDF(const std::string& file,bool convex,bool visualMesh) {
  ASSERT_MSGV(exists(file),"Cannot find %s",file.c_str())
  tinyxml2::XMLDocument pt;
  pt.LoadFile(file.c_str());
  if(convex)
    makeConvex(*(pt.RootElement()));
  //read links
  Links links;
  readLinks(links,*(pt.RootElement()),std::filesystem::path(file).parent_path(),visualMesh);
  //read joints
  Relations relations;
  if(!readRelations(relations,*(pt.RootElement())))
    return ArticulatedBody();
  //merge
  std::unordered_set<std::string> rootFilter;
  for(const std::pair<const std::string,std::vector<JointInfo>>& v:relations) {
    rootFilter.insert(v.first);
  }
  for(const std::pair<const std::string,std::vector<JointInfo>>& v:relations) {
    for(const JointInfo& c:v.second) {
      rootFilter.erase(c._child);
    }
  }
  ArticulatedBody body;
  ASSERT(rootFilter.size()==1)
  if(!buildBody(body,*(rootFilter.begin()),links,relations))
    return ArticulatedBody();
  body.fillChildren();
  //mimic
  for(const std::pair<const std::string,std::vector<JointInfo>>& v:relations) {
    for(const JointInfo& c:v.second) {
      if(!c._mimic.empty()) {
        JointInfo info;
        Joint& joint=body.joint(body.jointId(c._child));
        findRelationByName(c._mimic,relations,info);
        joint._mimic=body.jointId(info._child);
      }
    }
  }
  return body;
}
//built-in bodies
int ArticulatedLoader::inferDim(int rootType) {
  int dim=2;
  Joint::loopAllJointTypes([&](Joint::JOINT_TYPE t) {
    if(t&rootType) {
      if(t==Joint::TRANS_3D || t==Joint::ROT_3D_EXP || t==Joint::ROT_3D_XYZ || t==Joint::BALL_JOINT)
        dim=3;
    }
  });
  return dim;
}
tinyxml2::XMLElement* ArticulatedLoader::addRootJoint(tinyxml2::XMLElement& pt,int rootType) {
  int nrRotJoint=0;
  int nrTransJoint=0;
  Joint::JOINT_TYPE rotJoint=Joint::NR_JOINT_TYPE;
  Joint::JOINT_TYPE transJoint=Joint::NR_JOINT_TYPE;
  Joint::loopAllJointTypes([&](Joint::JOINT_TYPE t) {
    if(t&rootType) {
      Joint J;
      J._typeJoint=t;
      if(J.RBegEnd()[1]-J.RBegEnd()[0]) {
        rotJoint=t;
        nrRotJoint++;
      } else if(J.CBegEnd()[1]-J.CBegEnd()[0]) {
        transJoint=t;
        nrTransJoint++;
      }
    }
  });
  ASSERT_MSG(nrRotJoint<=1 && nrTransJoint<=1,"You can only have at most 1 rotJoint and 1 transJoint!")
  tinyxml2::XMLElement* ret=NULL;
  if(nrTransJoint>0) {
    ret=addChild(ret?*ret:pt,"joint");
    put<std::string>(*ret,"type",Joint::typeToString(transJoint));
  }
  if(nrRotJoint>0) {
    ret=addChild(ret?*ret:pt,"joint");
    put<std::string>(*ret,"type",Joint::typeToString(rotJoint));
  }
  if(!ret) {
    ret=addChild(ret?*ret:pt,"joint");
    put<std::string>(*ret,"type",Joint::typeToString(Joint::FIX_JOINT));
  }
  return ret;
}
tinyxml2::XMLElement* ArticulatedLoader::createBird(tinyxml2::XMLElement& pt,int rootType,T bodySz,T bodyLen,T neckSz,T neckLen,T footLen1,T footLen2,T footLen3,T footRad1,T footRad2,T footRad3,bool fixFoot,bool head) {
  Mat3T R0=expWGradV<T,Vec3T>(Vec3T(0,0,D2R(-85.0f)));
  Mat3T R1=expWGradV<T,Vec3T>(Vec3T(0,0,D2R(-45.0f)));
  Mat3T R2=expWGradV<T,Vec3T>(Vec3T(0,0,D2R(30.0f)));
  Mat3T R3=expWGradV<T,Vec3T>(Vec3T(0,0,D2R(-10.0f)));
  tinyxml2::XMLElement* body=addRootJoint(pt,rootType);
  put<T>(*body,"geom1D.lenY",bodyLen);
  put<T>(*body,"geom1D.rad",bodySz);
  putPtree<Vec3T>(*body,"meshDirY",R0*Vec3T::UnitY());
  for(int i=0; i<2; i++) {
    //leg1
    T sgn=i==0?1:-1;
    tinyxml2::XMLElement* leg1=addChild(*body,"joint");
    put<std::string>(*leg1,"name",std::string(i==0?"R":"L")+"_UpperLimb");
    put<std::string>(*leg1,"type",Joint::typeToString(Joint::ROT_3D_XYZ));
    putPtree(*leg1,"limit.lower",Vec3T(D2R(-60.0f),D2R(-60.0f),D2R(-45.0f)));
    putPtree(*leg1,"limit.upper",Vec3T(D2R(60.0f),D2R(60.0f),D2R(45.0f)));
    put<T>(*leg1,"geom1D.lenY",footLen1);
    put<T>(*leg1,"geom1D.rad",footRad1);
    putPtree<Vec3T>(*leg1,"trans",Vec3T::UnitZ()*bodySz*2/3*sgn-Vec3T::UnitY()*bodySz/2+R0*Vec3T::UnitY()*bodySz);
    putPtree<Vec3T>(*leg1,"jointDirY",Vec3T(0,0,1));
    putPtree<Vec3T>(*leg1,"jointDirZ",R1*Vec3T(1,0,0)*sgn);
    putPtree<Vec3T>(*leg1,"meshDirY",R1*Vec3T(0,-1,0));

    //leg2
    tinyxml2::XMLElement* leg2=addChild(*leg1,"joint");
    put<std::string>(*leg2,"name",std::string(i==0?"R":"L")+"_LowerLimb");
    put<std::string>(*leg2,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree<Vec>(*leg2,"limit.lower",Vec::Constant(1,D2R(-90.0f)));
    putPtree<Vec>(*leg2,"limit.upper",Vec::Constant(1,D2R(90.0f)));
    //putPtree(*leg2,"limit.coef",Vec::Constant(nrDOF,1000));  //use global joint limit
    put<T>(*leg2,"geom1D.lenY",footLen2);
    put<T>(*leg2,"geom1D.rad",footRad2);
    put<T>(*leg2,"transY",1);
    putPtree<Vec3T>(*leg2,"jointDirZ",Vec3T(0,0,1));
    putPtree<Vec3T>(*leg2,"meshDirY",R2*Vec3T(0,-1,0));

    //leg3
    tinyxml2::XMLElement* leg3=addChild(*leg2,"joint");
    put<std::string>(*leg3,"name",std::string(i==0?"R":"L")+"_Foot");
    if(fixFoot)
      put<std::string>(*leg3,"type",Joint::typeToString(Joint::FIX_JOINT));
    else {
      put<std::string>(*leg3,"type",Joint::typeToString(Joint::BALL_JOINT));
      putPtree<Vec>(*leg3,"limit.lower",Vec2T(D2R(-45.0f),D2R(-45.0f)));
      putPtree<Vec>(*leg3,"limit.upper",Vec2T(D2R(45.0f),D2R(70.0f)));
    }
    //putPtree(*leg3,"limit.coef",Vec::Constant(nrDOF,1000));  //use global joint limit
    put<T>(*leg3,"geom3D.lenX",footRad3);
    put<T>(*leg3,"geom3D.lenY",footLen3);
    put<T>(*leg3,"geom3D.lenZ",footRad3/2);
    put<T>(*leg3,"geom3D.rad",0);
    put<T>(*leg3,"transY",1);
    putPtree<Vec3T>(*leg3,"trans",-Vec3T::UnitY()*footRad2);
    putPtree<Vec3T>(*leg3,"jointDirY",Vec3T(1,0,0));
    putPtree<Vec3T>(*leg3,"jointDirZ",Vec3T(0,0,1));
    putPtree<Vec3T>(*leg3,"meshDirX",Vec3T(0,0,1));
    putPtree<Vec3T>(*leg3,"meshDirY",Vec3T(1,0,0));
    putPtree<Vec3T>(*leg3,"meshTrans",-Vec3T::UnitY()*footLen3/2);
  }
  if(head) {
    //neck
    tinyxml2::XMLElement* neck=addChild(*body,"joint");
    put<std::string>(*neck,"name","neck");
    put<std::string>(*neck,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree<Vec>(*neck,"limit.lower",Vec::Constant(1,D2R(-90.0f)));
    putPtree<Vec>(*neck,"limit.upper",Vec::Constant(1,D2R(90.0f)));
    put<T>(*neck,"geom1D.lenY",neckLen);
    put<T>(*neck,"geom1D.rad",neckSz);
    putPtree<Vec3T>(*neck,"trans",Vec3T::UnitY()*bodySz/2+R0*Vec3T::UnitY()*(bodyLen+bodySz/2));
    putPtree<Vec3T>(*neck,"jointDirZ",Vec3T(0,0,1));
    putPtree<Vec3T>(*neck,"meshDirY",Vec3T(0,1,0));

    //head
    tinyxml2::XMLElement* head=addChild(*neck,"joint");
    put<std::string>(*head,"name","head");
    put<std::string>(*head,"type",Joint::typeToString(Joint::FIX_JOINT));
    put<T>(*head,"geom1D.lenY",neckSz*2);
    put<T>(*head,"geom1D.rad",neckSz*1.5f);
    put<T>(*head,"transY",1);
    putPtree<Vec3T>(*head,"meshDirY",R3*Vec3T(1,0,0));
  }
  return &pt;
}
tinyxml2::XMLElement* ArticulatedLoader::createBipedal(tinyxml2::XMLElement& pt,int rootType,T bodySz,T bodyLen,T footLen1,T footLen2,T footLen3,T footRad1,T footRad2,T footRad3,bool fixFoot) {
  tinyxml2::XMLElement* body=addRootJoint(pt,rootType);
  put<T>(*body,"geom1D.lenY",bodyLen);
  put<T>(*body,"geom1D.rad",bodySz);
  putPtree<Vec3T>(*body,"meshDirY",Vec3T::UnitY());
  for(int i=0; i<2; i++) {
    //leg1
    T sgn=i==0?1:-1;
    tinyxml2::XMLElement* leg1=addChild(*body,"joint");
    put<std::string>(*leg1,"name",std::string(i==0?"R":"L")+"_UpperLimb");
    put<std::string>(*leg1,"type",Joint::typeToString(Joint::ROT_3D_XYZ));
    putPtree(*leg1,"limit.lower",Vec3T(D2R(-45.0f),D2R(-45.0f),D2R(-60.0f)));
    putPtree(*leg1,"limit.upper",Vec3T(D2R(45.0f),D2R(150.0f),D2R(10.0f)));
    put<T>(*leg1,"geom1D.lenY",footLen1);
    put<T>(*leg1,"geom1D.rad",footRad1);
    putPtree<Vec3T>(*leg1,"trans",Vec3T::UnitZ()*bodySz*2/3*sgn-Vec3T::UnitY()*bodySz/2);
    putPtree<Vec3T>(*leg1,"jointDirY",Vec3T(0,0,1));
    putPtree<Vec3T>(*leg1,"jointDirZ",Vec3T(1,0,0)*sgn);
    putPtree<Vec3T>(*leg1,"meshDirY",Vec3T(0,-1,0));

    //leg2
    tinyxml2::XMLElement* leg2=addChild(*leg1,"joint");
    put<std::string>(*leg2,"name",std::string(i==0?"R":"L")+"_LowerLimb");
    put<std::string>(*leg2,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree<Vec>(*leg2,"limit.lower",Vec::Constant(1,D2R(-150.0f)));
    putPtree<Vec>(*leg2,"limit.upper",Vec::Constant(1,D2R(0.0f)));
    //putPtree(*leg2,"limit.coef",Vec::Constant(nrDOF,1000));  //use global joint limit
    put<T>(*leg2,"geom1D.lenY",footLen2);
    put<T>(*leg2,"geom1D.rad",footRad2);
    put<T>(*leg2,"transY",1);
    putPtree<Vec3T>(*leg2,"jointDirZ",Vec3T(0,0,1));
    putPtree<Vec3T>(*leg2,"meshDirY",Vec3T(0,-1,0));

    //leg3
    tinyxml2::XMLElement* leg3=addChild(*leg2,"joint");
    put<std::string>(*leg3,"name",std::string(i==0?"R":"L")+"_Foot");
    if(fixFoot)
      put<std::string>(*leg3,"type",Joint::typeToString(Joint::FIX_JOINT));
    else {
      put<std::string>(*leg3,"type",Joint::typeToString(Joint::BALL_JOINT));
      putPtree<Vec>(*leg3,"limit.lower",Vec2T(D2R(-45.0f),D2R(-45.0f)));
      putPtree<Vec>(*leg3,"limit.upper",Vec2T(D2R(45.0f),D2R(45.0f)));
    }
    //putPtree(*leg3,"limit.coef",Vec::Constant(nrDOF,1000));  //use global joint limit
    put<T>(*leg3,"geom3D.lenX",footRad3);
    put<T>(*leg3,"geom3D.lenY",footLen3);
    put<T>(*leg3,"geom3D.lenZ",footRad3/2);
    put<T>(*leg3,"geom3D.rad",0);
    put<T>(*leg3,"transY",1);
    putPtree<Vec3T>(*leg3,"trans",-Vec3T::UnitY()*footRad2);
    putPtree<Vec3T>(*leg3,"jointDirY",Vec3T(1,0,0));
    putPtree<Vec3T>(*leg3,"jointDirZ",Vec3T(0,0,1));
    putPtree<Vec3T>(*leg3,"meshDirX",Vec3T(0,0,1));
    putPtree<Vec3T>(*leg3,"meshDirY",Vec3T(1,0,0));
    putPtree<Vec3T>(*leg3,"meshTrans",-Vec3T::UnitY()*footRad2);
  }
  return body;
}
void ArticulatedLoader::addBipedalHand(tinyxml2::XMLElement& pt,T bodySz,T bodyLen,T handLen1,T handLen2,T handRad1,T handRad2) {
  tinyxml2::XMLElement* body=&pt;
  for(int i=0; i<2; i++) {
    //leg1
    T sgn=i==0?1:-1;
    tinyxml2::XMLElement* hand1=addChild(*body,"joint");
    put<std::string>(*hand1,"name",std::string(i==0?"R":"L")+"_UpperArm");
    put<std::string>(*hand1,"type",Joint::typeToString(Joint::BALL_JOINT));
    putPtree(*hand1,"limit.lower",Vec2T(D2R(-60.0f),D2R(-90.0f)));
    putPtree(*hand1,"limit.upper",Vec2T(D2R(90.0f),D2R(0.0f)));
    put<T>(*hand1,"geom1D.lenY",handLen1);
    put<T>(*hand1,"geom1D.rad",handRad1);
    putPtree<Vec3T>(*hand1,"trans",Vec3T::UnitZ()*(bodySz+handRad1)*sgn+Vec3T::UnitY()*bodyLen);
    putPtree<Vec3T>(*hand1,"jointDirY",Vec3T(0,0,1));
    putPtree<Vec3T>(*hand1,"jointDirZ",Vec3T(1,0,0)*sgn);
    putPtree<Vec3T>(*hand1,"meshDirY",-Vec3T::UnitY());

    //leg2
    tinyxml2::XMLElement* hand2=addChild(*hand1,"joint");
    put<std::string>(*hand2,"name",std::string(i==0?"R":"L")+"_LowerArm");
    put<std::string>(*hand2,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree<Vec>(*hand2,"limit.lower",Vec::Constant(1,D2R(0.0f)));
    putPtree<Vec>(*hand2,"limit.upper",Vec::Constant(1,D2R(150.0f)));
    //putPtree(*hand2,"limit.coef",Vec::Constant(nrDOF,1000));  //use global joint limit
    put<T>(*hand2,"geom1D.lenY",handLen2);
    put<T>(*hand2,"geom1D.rad",handRad2);
    put<T>(*hand2,"transY",1);
    putPtree<Vec3T>(*hand2,"jointDirZ",Vec3T(0,0,1));
    putPtree<Vec3T>(*hand2,"meshDirY",Vec3T(0,-1,0));
  }
}
tinyxml2::XMLElement* ArticulatedLoader::createChain(tinyxml2::XMLElement& pt,int rootType,int nr,T l,T rad,T rot,T rot0,T t,T t0,int geomDim,T yOff,T ratioX,T ratioZ) {
  int nrDOF;
  int dim=inferDim(rootType);
  tinyxml2::XMLElement* child=NULL;
  for(int i=0; i<nr; i++) {
    if(i == 0)
      child=addRootJoint(pt,rootType);
    else {
      nrDOF=dim == 2 ? 1 : 2;
      put<std::string>(*child,"type",dim == 2 ? "HINGE_JOINT" : "BALL_JOINT");
      putPtree<Vec>(*child,"limit.lower",Vec::Constant(nrDOF,-rot));
      putPtree<Vec>(*child,"limit.upper",Vec::Constant(nrDOF,rot));
      //putPtree<Vec>(*child,"limit.coef",Vec::Constant(nrDOF,1000));  //use global joint limit
      putPtree<Vec>(*child,"x0",Vec::Constant(nrDOF,rot0));
    }
    //geometry
    put<T>(*child,"geom"+std::to_string(geomDim)+"D.lenX",l*ratioX);
    put<T>(*child,"geom"+std::to_string(geomDim)+"D.lenY",l);
    put<T>(*child,"geom"+std::to_string(geomDim)+"D.lenZ",l*ratioZ);
    put<T>(*child,"geom"+std::to_string(geomDim)+"D.rad",rad);
    //rotation
    put<T>(*child,"transY",1+yOff);
    putPtree(*child,"meshDirX",Vec3T(1,0,0));
    putPtree(*child,"meshDirY",Vec3T(0,-1,0));
    putPtree(*child,"jointDirY",Vec3T(0,0,1));
    putPtree(*child,"jointDirZ",Vec3T(1,0,0));
    //translation
    if(t > 0 && i<nr-1) {
      child=addChild(*child,"joint");
      put<std::string>(*child,"type",Joint::typeToString(Joint::TRANS_1D));
      put<T>(*child,"limit.lower",-t);
      put<T>(*child,"limit.upper",t);
      //put<T>(*child,"limit.coef",1000);  //use global joint limit
      put<T>(*child,"transY",1);
      put<T>(*child,"x0",t0);
      putPtree(*child,"jointDirX",Vec3T(0,-1,0));
    }
    if(i<nr-1)
      child=addChild(*child,"joint");
  }
  return &pt;
}
tinyxml2::XMLElement* ArticulatedLoader::createSpider(tinyxml2::XMLElement& pt,int rootType,T bodySz,T footLen,T footRad,T angLegLift,T rangeLeg1,T rangeLeg2,bool ball) {
  int nrDOF;
  tinyxml2::XMLElement* body=addRootJoint(pt,rootType);
  put<T>(*body,"geom0D.rad",bodySz);
  for(int i=0; i<4; i++) {
    Mat3T R0=expWGradV<T,Vec3T>(Vec3T(0,M_PI/4+(T)i*M_PI/2,0));
    Mat3T R1=expWGradV<T,Vec3T>(Vec3T(0,0,angLegLift));
    Mat3T R2=expWGradV<T,Vec3T>(Vec3T(0,0,angLegLift-M_PI/2));

    //leg1
    nrDOF=ball?2:1;
    tinyxml2::XMLElement* leg1=addChild(*body,"joint");
    put<std::string>(*leg1,"name",std::string(i==0?"LF":i==1?"RF":i==2?"RB":"LB")+"_UpperLimb");
    put<std::string>(*leg1,"type",ball ? "BALL_JOINT" : "HINGE_JOINT");
    putPtree<Vec>(*leg1,"limit.lower",Vec::Constant(nrDOF,-rangeLeg1));
    putPtree<Vec>(*leg1,"limit.upper",Vec::Constant(nrDOF,rangeLeg1));
    //putPtree(*leg1,"limit.coef",Vec::Constant(nrDOF,1000));  //use global joint limit
    put<T>(*leg1,"geom1D.lenY",footLen);
    put<T>(*leg1,"geom1D.rad",footRad);
    putPtree<Vec3T>(*leg1,"trans",R0*Vec3T::UnitX()*bodySz);
    putPtree<Vec3T>(*leg1,"jointDirY",R0*Vec3T(0,0,1));
    putPtree<Vec3T>(*leg1,"jointDirZ",Vec3T(0,1,0));
    putPtree<Vec3T>(*leg1,"meshDirY",R0*R1*Vec3T::UnitX());

    //leg2
    nrDOF=1;
    tinyxml2::XMLElement* leg2=addChild(*leg1,"joint");
    put<std::string>(*leg2,"name",std::string(i==0?"LF":i==1?"RF":i==2?"RB":"LB")+"_LowerLimb");
    put<std::string>(*leg2,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree<Vec>(*leg2,"limit.lower",Vec::Constant(nrDOF,-rangeLeg2));
    putPtree<Vec>(*leg2,"limit.upper",Vec::Constant(nrDOF,rangeLeg2));
    //putPtree(*leg2,"limit.coef",Vec::Constant(nrDOF,1000));  //use global joint limit
    put<T>(*leg2,"geom1D.lenY",footLen);
    put<T>(*leg2,"geom1D.rad",footRad);
    put<T>(*leg2,"transY",1);
    putPtree<Vec3T>(*leg2,"jointDirZ",R0*Vec3T(0,0,1));
    putPtree<Vec3T>(*leg2,"meshDirY",R0*R2*Vec3T::UnitX());
  }
  return &pt;
}
tinyxml2::XMLElement* ArticulatedLoader::createDummy(tinyxml2::XMLElement&pt) {
  tinyxml2::XMLElement* root=addChild(pt,"joint");
  put<std::string>(*root,"type","FIX_JOINT");
  return &pt;
}
tinyxml2::XMLElement* ArticulatedLoader::createArm(tinyxml2::XMLElement& pt,T armLen,T armRad) {
  tinyxml2::XMLElement* root=addChild(pt,"joint");
  put<std::string>(*root,"type","FIX_JOINT");
  {
    //the table
    tinyxml2::XMLElement* child=addChild(*root,"joint");
    put<std::string>(*child,"type","FIX_JOINT");
    put<T>(*child,"geom3D.lenX",armRad*3);
    put<T>(*child,"geom3D.lenY",armRad);
    put<T>(*child,"geom3D.lenZ",armRad*3);
    put<T>(*child,"geom3D.rad",0);
    putPtree<Vec3T>(*child,"trans",-Vec3T::UnitY()*armRad*1.5f);
    putPtree<Vec3T>(*child,"meshDirX",Vec3T::UnitX());
    putPtree<Vec3T>(*child,"meshDirY",Vec3T::UnitY());

    tinyxml2::XMLElement* child2=addChild(*root,"joint");
    put<std::string>(*child2,"type","FIX_JOINT");
    put<T>(*child2,"geom3D.lenX",armLen*10);
    put<T>(*child2,"geom3D.lenY",armRad);
    put<T>(*child2,"geom3D.lenZ",armLen*10);
    put<T>(*child2,"geom3D.rad",0);
    putPtree<Vec3T>(*child2,"trans",-Vec3T::UnitY()*armRad*2.5f);
    putPtree<Vec3T>(*child2,"meshDirX",Vec3T::UnitX());
    putPtree<Vec3T>(*child2,"meshDirY",Vec3T::UnitY());
  }
  {
    //the arm
    tinyxml2::XMLElement* child=addChild(*root,"joint");
    put<std::string>(*child,"type","BALL_JOINT");
    putPtree(*child,"limit.lower",Vec2T(-D2R(180.0f),D2R(10.0f)));
    putPtree(*child,"limit.upper",Vec2T(D2R(180.0f),D2R(90.0f)));
    putPtree<Vec3T>(*child,"jointDirY",Vec3T(0,1,0));
    putPtree<Vec3T>(*child,"jointDirZ",Vec3T(0,0,1));
    put<T>(*child,"geom1D.lenY",armLen);
    put<T>(*child,"geom1D.rad",armRad);
    putPtree<Vec3T>(*child,"meshDirY",Vec3T::UnitY());

    tinyxml2::XMLElement* child2=addChild(*child,"joint");
    put<std::string>(*child2,"type","HINGE_JOINT");
    putPtree<Vec>(*child2,"limit.lower",Vec::Constant(1,D2R(0.0f)));
    putPtree<Vec>(*child2,"limit.upper",Vec::Constant(1,D2R(150.0f)));
    put<T>(*child2,"transY",1);
    putPtree<Vec3T>(*child2,"jointDirZ",Vec3T(0,0,1));
    put<T>(*child2,"geom1D.lenY",armLen);
    put<T>(*child2,"geom1D.rad",armRad);
    putPtree<Vec3T>(*child2,"meshDirY",Vec3T::UnitY());
  }
  return &pt;
}
tinyxml2::XMLElement* ArticulatedLoader::createBox(tinyxml2::XMLElement& pt,T x,T y,T z,T rad) {
  tinyxml2::XMLElement* root=addChild(pt,"joint");
  put<std::string>(*root,"type","FIX_JOINT");
  {
    //the table
    tinyxml2::XMLElement* child=addChild(*root,"joint");
    put<std::string>(*child,"type","FIX_JOINT");
    put<T>(*child,"geom3D.lenX",x);
    put<T>(*child,"geom3D.lenY",y);
    put<T>(*child,"geom3D.lenZ",z);
    put<T>(*child,"geom3D.rad",rad);
    putPtree<Vec3T>(*child,"meshDirX",Vec3T::UnitX());
    putPtree<Vec3T>(*child,"meshDirY",Vec3T::UnitY());
  }
  return &pt;
}
//build-in bodies as ArticulatedBody
ArticulatedBody ArticulatedLoader::createBird(int root,bool head) {
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  createBird(*(pt.RootElement()),root,0.25f,0.5f,0.1f,0.6f, 0.3f,0.4f,0.45f, 0.12f,0.1f,0.2f,false,head);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createBipedal(int root,bool withHand) {
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  tinyxml2::XMLElement* bodyNode=createBipedal(*(pt.RootElement()),root,0.25f,0.5f, 0.4f,0.4f,0.3f, 0.12f,0.1f,0.15f, false);
  if(withHand)
    addBipedalHand(*bodyNode,0.25f,0.5f, 0.4f,0.4f, 0.1f,0.1f);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createChain(int root,T rad,int nrLink) {
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  createChain(*(pt.RootElement()),root,nrLink,0.5f,rad,D2R(120),0,0,0,3,0.0f,0.2f,0.2f);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createSpider(int root,T lmt,T lmt2) {
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  T footLen=0.2f*sqrt(2.0f)+0.16f;
  ArticulatedLoader::createSpider(*(pt.RootElement()),root,0.25f,footLen,0.08f,D2R(10),D2R(lmt),D2R(lmt2),true);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createArm(T armLen,T armRad) {
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  createArm(*(pt.RootElement()),armLen,armRad);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createDummy() {
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  createDummy(*(pt.RootElement()));
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
}
