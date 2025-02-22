/*#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Environment/MeshExact.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/PointCloudExact.h>
#include <Environment/SphericalBBoxExact.h>
#include <Environment/CompositeShapeExact.h>
#include <Simulator/PBDSimulator.h>
#include <Simulator/XPBDSimulator.h>
#include <Simulator/ConvHullPBDSimulator.h>
#include <Simulator/MeshBasedPBDSimulator.h>
#include <Simulator/PBDMatrixSolver.h>
#include <Simulator/JointLimit.h>
#include <Simulator/SoftJoint.h>
#include <Utils/RotationUtils.h>
#include <Utils/CrossSpatialUtils.h>

using namespace PHYSICSMOTION;
namespace py=pybind11;

typedef double T;
DECL_MAT_VEC_MAP_TYPES_T
void initializeJoint(py::module& m) {
  //Joint::JOINT_TYPE
  py::enum_<Joint::JOINT_TYPE>(m, "JOINT_TYPE")
  .value("TRANS_3D", Joint::TRANS_3D)
  .value("TRANS_2D", Joint::TRANS_2D)
  .value("TRANS_1D", Joint::TRANS_1D)
  .value("ROT_3D_XYZ", Joint::ROT_3D_XYZ)
  .value("ROT_3D_EXP", Joint::ROT_3D_EXP)
  .value("BALL_JOINT", Joint::BALL_JOINT)
  .value("HINGE_JOINT", Joint::HINGE_JOINT)
  .value("FIX_JOINT", Joint::FIX_JOINT)
  .value("NR_JOINT_TYPE", Joint::NR_JOINT_TYPE);
  //Joint
  py::class_<Joint,
  std::shared_ptr<Joint>>(m,"Joint")
  .def(py::init<>())
  .def(py::init<const Joint&>())
  .def("setType",&Joint::setType)
  .def("typeToString",&Joint::typeToString)
  .def("getBB",&Joint::getBB)
  .def("getAxes",[](std::shared_ptr<Joint> joint)->std::tuple<Joint::Mat3XT,bool,bool,bool> {
    bool markX,markY,markZ;
    Joint::Mat3XT ret=joint->getAxes(markX,markY,markZ);
    return std::make_tuple(ret,markX,markY,markZ);
  })
  .def("getAxes",[](std::shared_ptr<Joint> joint,Joint::Mat3XT t)->std::tuple<Joint::Mat3XT,bool,bool,bool> {
    bool markX,markY,markZ;
    Joint::Mat3XT ret=joint->getAxes(markX,markY,markZ,&t);
    return std::make_tuple(ret,markX,markY,markZ);
  })
  .def("assemble",&Joint::assemble)
  .def("debugTransformMesh",&Joint::debugTransformMesh)
  .def("transformMass",&Joint::transformMass)
  .def("transformMesh",static_cast<void(Joint::*)(const Joint::Mat3X4T&)>(&Joint::transformMesh))
  .def("transformMesh",static_cast<void(Joint::*)(const Joint::Mat3T&,const Joint::Vec3T&)>(&Joint::transformMesh))
  .def("transformMesh",static_cast<void(Joint::*)(const Joint::Vec3T&)>(&Joint::transformMesh))
  .def("CBegEnd",&Joint::CBegEnd)
  .def("RBegEnd",&Joint::RBegEnd)
  .def("nrDOF",&Joint::nrDOF)
  .def("nrDDT",&Joint::nrDDT)
  .def("getMassC",static_cast<Joint::Mat6T(Joint::*)(const Joint::Mat3T&)const>(&Joint::getMassC))
  .def("getMassC",static_cast<Joint::Mat6T(Joint::*)()const>(&Joint::getMassC))
  .def("getMass",static_cast<Joint::Mat6T(Joint::*)(const Joint::Mat3T&)const>(&Joint::getMass))
  .def("getMass",static_cast<Joint::Mat6T(Joint::*)()const>(&Joint::getMass))
  .def("getC",&Joint::getC)
  .def("isRotational",&Joint::isRotational)
  .def("isRoot",&Joint::isRoot)
  .def("getGeomPtr",&Joint::getGeomPtr)
  .def("initL1",[](std::shared_ptr<Joint> joint,std::shared_ptr<ArticulatedBody> body) {
    return joint->initL1(*body);
  })
  .def("getL1",[](std::shared_ptr<Joint> joint,int vertexId)->Joint::Vec{return joint->getL1(vertexId);})
  //members
  .def_readwrite("children",&Joint::_children)
  .def_readwrite("parent",&Joint::_parent)
  .def_readwrite("depth",&Joint::_depth)
  .def_readwrite("typeJoint",&Joint::_typeJoint)
  .def_readwrite("mimic",&Joint::_mimic)
  .def_readwrite("offDOF",&Joint::_offDOF)
  .def_readwrite("offDDT",&Joint::_offDDT)
  .def_readwrite("class",&Joint::_class)
  .def_readwrite("limits",&Joint::_limits)
  .def_readwrite("control",&Joint::_control)
  .def_readwrite("damping",&Joint::_damping)
  .def_readwrite("trans",&Joint::_trans)
  .def_readwrite("mult",&Joint::_mult)
  .def_readwrite("offset",&Joint::_offset)
  .def_readwrite("M",&Joint::_M)
  .def_readwrite("MC",&Joint::_MC)
  .def_readwrite("MCCT",&Joint::_MCCT)
  .def_readwrite("name",&Joint::_name)
  .def_readwrite("mesh",&Joint::_mesh)
  .def_readwrite("L1",&Joint::_L1);
}
void initializeArticulatedBody(py::module& m) {
  //ArticulatedBody
  py::class_<ArticulatedBody,
  std::shared_ptr<ArticulatedBody>>(m,"ArticulatedBody")
  .def(py::init<>())
  .def(py::init<const tinyxml2::XMLElement&>())
  .def("randomize",&ArticulatedBody::randomize)
  .def("children",&ArticulatedBody::children)
  .def("commonRoot",&ArticulatedBody::commonRoot)
  .def("hasJoint",&ArticulatedBody::hasJoint)
  .def("resizeJoints",&ArticulatedBody::resizeJoints)
  .def("isLeaf",&ArticulatedBody::isLeaf)
  .def("control",&ArticulatedBody::control)
  .def("damping",&ArticulatedBody::damping)
  .def("coefLimit",&ArticulatedBody::coefLimit)
  .def("lowerLimit",static_cast<ArticulatedBody::Vec(ArticulatedBody::*)()const>(&ArticulatedBody::lowerLimit))
  .def("upperLimit",static_cast<ArticulatedBody::Vec(ArticulatedBody::*)()const>(&ArticulatedBody::upperLimit))
  .def("lowerLimit",static_cast<ArticulatedBody::Vec(ArticulatedBody::*)(ArticulatedBody::T)const>(&ArticulatedBody::lowerLimit))
  .def("upperLimit",static_cast<ArticulatedBody::Vec(ArticulatedBody::*)(ArticulatedBody::T)const>(&ArticulatedBody::upperLimit))
  .def("clampLimit",&ArticulatedBody::clampLimit)
  .def("randomPose",&ArticulatedBody::randomPose)
  .def("getT",&ArticulatedBody::getT)
  .def("getBB",static_cast<BBoxExact(ArticulatedBody::*)(const ArticulatedBody::Vec&)const>(&ArticulatedBody::getBB))
  .def("getBB",static_cast<BBoxExact(ArticulatedBody::*)(const ArticulatedBody::Mat3XT&)const>(&ArticulatedBody::getBB))
  .def("writeVTK",static_cast<void(ArticulatedBody::*)(const std::string&,const ArticulatedBody::Mat3XT&)const>(&ArticulatedBody::writeVTK))
  .def("mimic",static_cast<ArticulatedBody::Vec(ArticulatedBody::*)(ArticulatedBody::Vec)const>(&ArticulatedBody::mimic))
  .def("mimic",[](std::shared_ptr<ArticulatedBody> body)->std::tuple<ArticulatedBody::MatT,ArticulatedBody::Vec,ArticulatedBody::Vec,ArticulatedBody::Vec> {
    ArticulatedBody::MatT A;
    ArticulatedBody::Vec b,l,u;
    body->mimic(A,b,l,u);
    return std::make_tuple(A,b,l,u);
  })
  .def("movable",&ArticulatedBody::movable)
  .def("jointFromDOF",[](std::shared_ptr<ArticulatedBody> body,int id) {
    for(int j=0; j<body->nrJ(); j++)
      if(id>=body->joint(j)._offDOF && id-body->joint(j)._offDOF<body->joint(j).nrDOF())
        return body->jointSmartPtr(j);
    ASSERT_MSGV(false,"Cannot find joint from DOF %d",id)
    return body->jointSmartPtr(0);
  })
  .def("getJoint",[](std::shared_ptr<ArticulatedBody> body,int id) {
    return body->jointSmartPtr(id);
  })
  .def("setJoint",[](std::shared_ptr<ArticulatedBody> body,int id,std::shared_ptr<Joint> joint) {
    body->joint(id)=*joint;
  })
  .def("jointId",&ArticulatedBody::jointId)
  .def("rootJointId",&ArticulatedBody::rootJointId)
  .def("depth",&ArticulatedBody::depth)
  .def("nrDOF",&ArticulatedBody::nrDOF)
  .def("nrDDT",&ArticulatedBody::nrDDT)
  .def("nrJ",&ArticulatedBody::nrJ)
  .def("fillChildren",&ArticulatedBody::fillChildren)
  .def("setRootTrans",&ArticulatedBody::setRootTrans)
  .def("addBase",&ArticulatedBody::addBase)
  .def("simplify",&ArticulatedBody::simplify)
  .def("scaleMass",&ArticulatedBody::scaleMass)
  .def("totalMass",&ArticulatedBody::totalMass);
}
void initializeArticulatedUtils(py::module& m) {
  //ArticulatedUtils
  py::class_<ArticulatedUtils,
  std::shared_ptr<ArticulatedUtils>>(m,"ArticulatedUtils")
  .def(py::init<>([](std::shared_ptr<ArticulatedBody> body) {
    return std::shared_ptr<ArticulatedUtils>(new ArticulatedUtils(*body));
  }))
  .def("addBase",&ArticulatedUtils::addBase)
  .def("combine",[](std::shared_ptr<ArticulatedUtils> utils,const std::vector<std::shared_ptr<ArticulatedBody>>& bodies) {
    std::vector<ArticulatedBody> bodiesInternal;
    for(auto body:bodies)
      bodiesInternal.push_back(*body);
    utils->combine(bodiesInternal);
  })
  .def("mergeChildren",&ArticulatedUtils::mergeChildren)
  .def("fix",[](std::shared_ptr<ArticulatedUtils> utils,std::function<bool(int,std::shared_ptr<Joint>)> canFix,const Vec& DOF) {
    utils->fix([&](int jid,const Joint& joint) {
      return canFix(jid,std::shared_ptr<Joint>(new Joint(joint)));
    },DOF);
  })
  .def("eliminate",[](std::shared_ptr<ArticulatedUtils> utils,std::function<bool(int,std::shared_ptr<Joint>)> canEliminate,const Vec& DOF) {
    utils->eliminate([&](int jid,const Joint& joint) {
      return canEliminate(jid,std::shared_ptr<Joint>(new Joint(joint)));
    },DOF);
  })
  .def("simplify",[](std::shared_ptr<ArticulatedUtils> utils,std::function<bool(int,std::shared_ptr<Joint>)> canSimplify,const Vec& DOF,int nrDebug) {
    utils->simplify([&](int jid,const Joint& joint) {
      return canSimplify(jid,std::shared_ptr<Joint>(new Joint(joint)));
    },DOF,nrDebug);
  })
  .def("simplify",static_cast<ArticulatedUtils::Vec(ArticulatedUtils::*)(const ArticulatedUtils::Vec&,int)>(&ArticulatedUtils::simplify))
  .def("replaceJoint",&ArticulatedUtils::replaceJoint)
  .def("convexDecompose",static_cast<void(ArticulatedUtils::*)(ArticulatedUtils::T)>(&ArticulatedUtils::convexDecompose))
  .def("convexDecompose",static_cast<void(ArticulatedUtils::*)(std::vector<int>,int,ArticulatedUtils::T)>(&ArticulatedUtils::convexDecompose))
  .def("addBody",&ArticulatedUtils::addBody)
  .def("tessellate",&ArticulatedUtils::tessellate)
  .def("BBApproxiate",&ArticulatedUtils::BBApproxiate)
  .def("makeConvex",&ArticulatedUtils::makeConvex)
  .def("transformTorso",static_cast<ArticulatedUtils::Mat3X4T(ArticulatedUtils::*)(const ArticulatedUtils::Mat3X4T&)>(&ArticulatedUtils::transformTorso))
  .def("transformTorso",static_cast<ArticulatedUtils::Mat3X4T(ArticulatedUtils::*)(const ArticulatedUtils::Mat3T&)>(&ArticulatedUtils::transformTorso))
  .def("transformTorso",static_cast<ArticulatedUtils::Mat3X4T(ArticulatedUtils::*)(const ArticulatedUtils::Vec3T&)>(&ArticulatedUtils::transformTorso))
  .def("transformTorsoToCOM",&ArticulatedUtils::transformTorsoToCOM)
  .def("scaleBody",&ArticulatedUtils::scaleBody)
  .def("scaleMass",&ArticulatedUtils::scaleMass)
  .def("totalMass",&ArticulatedUtils::totalMass);
}
void initializeArticulatedLoader(py::module& m) {
  //ArticulatedLoader
  py::class_<ArticulatedLoader,
  std::shared_ptr<ArticulatedLoader>>(m,"ArticulatedLoader")
  .def_static("readURDF",[](const std::string& file,bool convex,bool visualMesh) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::readURDF(file,convex,visualMesh);
    return body;
  })
  //build-in bodies as ArticulatedBody
  //Bird
  .def_static("createBird",[](int root,T bodySz=0.25f,T bodyLen=0.5f,T neckSz=0.1f,T neckLen=0.6f,T footLen1=0.3f,T footLen2=0.4f,T footLen3=0.45f,T footRad1=0.12f,T footRad2=0.1f,T footRad3=0.2f,bool fixFoot=false,bool head=false) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createBird(root,bodySz,bodyLen,neckSz,neckLen,footLen1,footLen2,footLen3,footRad1,footRad2,footRad3,fixFoot,head);
    return body;
  },py::arg("root"),py::arg("bodySz")=0.25f,py::arg("bodyLen")=0.5f,py::arg("neckSz")=1.0f,py::arg("neckLen")=0.6f,py::arg("footLen1")=0.3f,py::arg("footLen2")=0.4f,py::arg("footLen3")=0.45f,py::arg("footRad1")=0.12f,py::arg("footRad2")=0.1f,py::arg("footRad3")=0.2f,py::arg("fixFoot")=false,py::arg("head")=false)
  //Bipedal
  .def_static("createBipedal",[](int root,T bodySz=0.25f,T bodyLen=0.5f,T footLen1=0.4f,T footLen2=0.4f,T footLen3=0.3f,T footRad1=0.12f,T footRad2=0.1f,T footRad3=0.15f,bool fixFoot=true,bool withHand=false,T handLen1=0.4f,T handLen2=0.4f,T handRad1=0.1f,T handRad2=0.1f) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createBipedal(root,bodySz,bodyLen,footLen1,footLen2,footLen3,footRad1,footRad2,footRad3,fixFoot,withHand,handLen1,handLen2,handRad1,handRad2);
    return body;
  },py::arg("root"),py::arg("bodySz")=0.25f,py::arg("bodyLen")=0.5f,py::arg("footLen1")=0.4f,py::arg("footLen2")=0.4f,py::arg("footLen3")=0.3f,py::arg("footRad1")=0.12f,py::arg("footRad2")=0.1f,py::arg("footRad3")=0.15f,py::arg("fixFoot")=true,py::arg("withHand")=false,py::arg("handLen1")=0.4f,py::arg("handLen2")=0.4f,py::arg("handRad1")=0.1f,py::arg("handRad2")=0.1f)
  //Chain
  .def_static("createChain",[](int root,int nr,T l=0.5f,T rad=0.1f,T rot=D2R(360),T rot0=0,T t=0,T t0=0,int geomDim=3,T yOff=0,T ratioX=1,T ratioZ=1) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createChain(root,nr,l,rad,rot,rot0,t,t0,geomDim,yOff,ratioX,ratioZ);
    return body;
  },py::arg("root"),py::arg("nr"),py::arg("l")=0.5f,py::arg("rad")=0.1f,py::arg("rot")=D2R(360),py::arg("rot0")=0,py::arg("t")=0,py::arg("t0")=0,py::arg("geomDim")=3,py::arg("yOff")=0,py::arg("ratioX")=1,py::arg("ratioZ")=1)
  //Grasper
  .def_static("createGrasper",[](T armRad) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createGrasper(armRad);
    return body;
  },py::arg("armRad"))
  //Spider
  .def_static("createSpider",[](int root,T bodySz=0.2f,T footLen=0.2f*sqrt(2.0f)+0.16f,T footRad=0.08f,T angLegLift=D2R(10),T rangeLeg1=D2R(45),T rangeLeg2=D2R(60),bool ball=true) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createSpider(root,bodySz,footLen,footRad,angLegLift,rangeLeg1,rangeLeg2,ball);
    return body;
  },py::arg("root"),py::arg("bodySz")=0.2f,py::arg("footLen")=0.2f*sqrt(2.0f)+0.16f,py::arg("footRad")=0.08f,py::arg("angLegLift")=D2R(10),py::arg("rangeLeg1")=D2R(45),py::arg("rangeLeg2")=D2R(60),py::arg("ball")=true)
  //Arm
  .def_static("createArm",[](T armLen,T armRad) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createArm(armLen,armRad);
    return body;
  })
  //Random Initialize
  .def_static("createInitJoints",[](int convexhulls,T sz) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createInitJoints(convexhulls,sz);
    return body;
  },py::arg("convexhulls"),py::arg("sz"))
  //Ball
  .def_static("createBall",[](int root,T Rad) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createBall(root,Rad);
    return body;
  },py::arg("root"),py::arg("Rad"))
  //Box
  .def_static("createBox",[](T x,T y,T z,T Rad) {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createBox(x,y,z,Rad);
    return body;
  },py::arg("x"),py::arg("y"),py::arg("z"),py::arg("Rad"))
  //Dummy
  .def_static("createDummy",[]() {
    std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
    *body=ArticulatedLoader::createDummy();
    return body;
  });
}
void initializePBDGradientInfo(py::module& m) {
  typedef PBDArticulatedGradientInfo<T> Info;
  py::class_<Info,std::shared_ptr<Info>>(m,"PBDArticulatedGradientInfo")
  .def(py::init<>([](std::shared_ptr<ArticulatedBody> body,const Info::Vec& x) {
    return std::shared_ptr<Info>(new Info(*body,x));
  }))
  .def_property_readonly("T",[](std::shared_ptr<Info> pos)->Info::Mat3XT{
    return Info::Mat3XT(pos->_TM);
  });
}
void initializeShapes(py::module& m) {
  //ShapeExact
  py::class_<ShapeExact,
  std::shared_ptr<ShapeExact>>(m,"ShapeExact")
  .def("getBB",&ShapeExact::getBB)
  .def("empty",&ShapeExact::empty)
  .def("getMesh",[](std::shared_ptr<ShapeExact> shape)->std::tuple<std::vector<ShapeExact::Vec3T>,std::vector<Eigen::Matrix<int,3,1>>> {
    std::vector<ShapeExact::Vec3T> vss;
    std::vector<Eigen::Matrix<int,3,1>> iss;
    shape->getMesh(vss,iss);
    return std::make_tuple(vss,iss);
  })
  .def("scale",&ShapeExact::scale)
  .def("closest",[](std::shared_ptr<ShapeExact> shape)->std::tuple<ShapeExact::Vec3T,ShapeExact::Vec3T,ShapeExact::Vec3T,ShapeExact::Mat3T,Eigen::Matrix<int,2,1>> {
    ShapeExact::Vec3T pt,n,normal;
    ShapeExact::Mat3T hessian;
    ShapeExact::T rad;
    Eigen::Matrix<int,2,1> feat;
    shape->closestInner(pt,n,normal,hessian,rad,feat);
    return std::make_tuple(pt,n,normal,hessian,feat);
  })
  .def("support",[](std::shared_ptr<ShapeExact> shape,ShapeExact::Vec3T D)->std::tuple<ShapeExact::Vec3T,int> {
    int id;
    ShapeExact::Vec3T ret=shape->support(D,id);
    return std::make_tuple(ret,id);
  })
  .def("project",&ShapeExact::project)
  .def("facets",&ShapeExact::facets)
  .def("edges",&ShapeExact::edges)
  .def("checkPositiveFacets",&ShapeExact::checkPositiveFacets);
  //MeshExact
  py::class_<MeshExact,std::shared_ptr<MeshExact>,ShapeExact>(m,"MeshExact")
  .def(py::init<>())
  .def(py::init<const std::string&,bool>())
  .def(py::init<std::vector<Eigen::Matrix<double,3,1>>&,
       std::vector<Eigen::Matrix<int,3,1>>&,bool>())
  .def(py::init<std::vector<Eigen::Matrix<double,3,1>>&,
       std::vector<Eigen::Matrix<double,2,1>>&,
       std::vector<Eigen::Matrix<int,3,1>>&,bool>())
  .def("bss",&MeshExact::bss)
  .def("vss",&MeshExact::vss)
  .def("tcss",&MeshExact::tcss)
  .def("iss",&MeshExact::iss)
  .def("translate",&MeshExact::translate)
  .def("transform",&MeshExact::transform)
  .def("moveMesh",&MeshExact::moveMesh)
  .def("setMesh",&MeshExact::setMesh);
  //ConvexHullExact
  py::class_<ConvexHullExact,std::shared_ptr<ConvexHullExact>,MeshExact>(m,"ConvexHullExact")
  .def(py::init<const std::string&>())
  .def(py::init<>([](std::shared_ptr<MeshExact> mesh) {
    return std::shared_ptr<ConvexHullExact>(new ConvexHullExact(*mesh));
  }))
  .def(py::init<const std::vector<Eigen::Matrix<double,3,1>>&>())
  .def("parityCheck",&ConvexHullExact::parityCheck)
  .def("ess",&ConvexHullExact::ess)
  .def("plane",&ConvexHullExact::plane)
  .def("nrPlane",&ConvexHullExact::nrPlane);
  //PointCloudExact
  py::class_<PointCloudExact,std::shared_ptr<PointCloudExact>,ShapeExact>(m,"PointCloudExact")
  .def(py::init<>())
  .def(py::init<>([](std::shared_ptr<MeshExact> mesh,PointCloudExact::T d0) {
    return std::shared_ptr<PointCloudExact>(new PointCloudExact(*mesh,d0));
  }))
  .def("vss",&PointCloudExact::vss)
  .def("translate",&PointCloudExact::translate)
  .def("transform",&PointCloudExact::transform);
  //BBoxExact
  py::class_<BBoxExact,std::shared_ptr<BBoxExact>,ShapeExact>(m,"BBoxExact")
  .def(py::init<T,T,T>())
  .def(py::init<const Vec3T&,const Vec3T&>())
  .def(py::init<const Vec3T&,const Vec3T&,const Vec3T&>())
  .def("setUnion",[](std::shared_ptr<BBoxExact> box,std::shared_ptr<BBoxExact> other) {
    box->setUnion(*other);
  })
  .def("setUnion",static_cast<void(BBoxExact::*)(const BBoxExact::Vec3T&)>(&BBoxExact::setUnion))
  .def("extendUnion",&BBoxExact::extendUnion)
  .def_property("minCorner",[](std::shared_ptr<BBoxExact> box) {
    return box->minCorner();
  },[](std::shared_ptr<BBoxExact> box,const BBoxExact::Vec3T& in) {
    box->minCorner()=in;
  })
  .def_property("maxCorner",[](std::shared_ptr<BBoxExact> box) {
    return box->maxCorner();
  },[](std::shared_ptr<BBoxExact> box,const BBoxExact::Vec3T& in) {
    box->maxCorner()=in;
  })
  .def("intersect",[](std::shared_ptr<BBoxExact> box,std::shared_ptr<BBoxExact> other) {
    return box->intersect(*other);
  })
  .def("contain",&BBoxExact::contain)
  .def("enlargedEps",&BBoxExact::enlargedEps)
  .def("enlarged",&BBoxExact::enlarged)
  .def("distToSqr",[](std::shared_ptr<BBoxExact> box,std::shared_ptr<BBoxExact> other) {
    return box->distToSqr(*other);
  })
  .def("distToSqr",static_cast<BBoxExact::T(BBoxExact::*)(const BBoxExact::Vec3T&)const>(&BBoxExact::distToSqr));
  //SphericalBBoxExact
  py::class_<SphericalBBoxExact,std::shared_ptr<SphericalBBoxExact>,BBoxExact>(m,"SphericalBBoxExact")
  .def(py::init<>())
  .def(py::init<T>())
  .def(py::init<T,T>())
  .def(py::init<T,T,T>())
  .def(py::init<T,T,T,T>())
  .def("contain",&SphericalBBoxExact::contain)
  .def("isRoundCube",&SphericalBBoxExact::isRoundCube)
  .def("isRoundBoard",&SphericalBBoxExact::isRoundBoard)
  .def("isCapsule",&SphericalBBoxExact::isCapsule)
  .def("isSphere",&SphericalBBoxExact::isSphere)
  .def("radius",&SphericalBBoxExact::radius);
  //CompositeShapeExact::Material
  py::class_<CompositeShapeExact::Material>(m,"Material")
  .def(py::init<>())
  .def_readonly("ambient",&CompositeShapeExact::Material::_ambient)
  .def_readonly("diffuse",&CompositeShapeExact::Material::_diffuse)
  .def_readonly("specular",&CompositeShapeExact::Material::_specular)
  .def_readonly("shininess",&CompositeShapeExact::Material::_shininess)
  .def_readonly("useWireframe",&CompositeShapeExact::Material::_useWireframe)
  .def_readonly("texFile",&CompositeShapeExact::Material::_texFile);
  //CompositeShapeExact
  py::class_<CompositeShapeExact,std::shared_ptr<CompositeShapeExact>,ShapeExact>(m,"CompositeShapeExact")
  .def(py::init<>())
  .def(py::init<const std::string&,bool>())
  .def(py::init<const std::vector<std::shared_ptr<ShapeExact>>&,
       const std::vector<Mat3X4T>&>())
  .def(py::init<const std::vector<std::shared_ptr<ShapeExact>>&>())
  .def("transform",&CompositeShapeExact::transform)
  .def("getGeoms",&CompositeShapeExact::getGeoms)
  .def("getTrans",&CompositeShapeExact::getTrans)
  .def("getMaterials",static_cast<const std::vector<CompositeShapeExact::Material>&(CompositeShapeExact::*)()const>(&CompositeShapeExact::getMaterials));
}
void initializeContact(py::module& m) {
  //ContactGenerator::STATUS
  py::enum_<ContactGenerator::STATUS>(m, "CONTACT_STATUS")
  .value("STATIC_STATIC", ContactGenerator::STATIC_STATIC)
  .value("STATIC_DYNAMIC", ContactGenerator::STATIC_DYNAMIC)
  .value("DYNAMIC_DYNAMIC", ContactGenerator::DYNAMIC_DYNAMIC);
  //ContactGenerator::ContactPoint
  py::class_<ContactGenerator::ContactPoint>(m,"ContactPoint")
  .def(py::init<>())
  .def("depth",&ContactGenerator::ContactPoint::depth)
  .def("swap",&ContactGenerator::ContactPoint::swap)
  .def("transform",&ContactGenerator::ContactPoint::transform)
  .def_readonly("ptA",&ContactGenerator::ContactPoint::_ptA)
  .def_readonly("ptB",&ContactGenerator::ContactPoint::_ptB)
  .def_readonly("nA2B",&ContactGenerator::ContactPoint::_nA2B)
  .def_readonly("tangentBound",&ContactGenerator::ContactPoint::_tangentBound)
  .def_readonly("ptALast",&ContactGenerator::ContactPoint::_ptALast)
  .def_readonly("ptBLast",&ContactGenerator::ContactPoint::_ptBLast)
  .def_readonly("ptAL",&ContactGenerator::ContactPoint::_ptAL)
  .def_readonly("ptBL",&ContactGenerator::ContactPoint::_ptBL)
  .def_readonly("fA",&ContactGenerator::ContactPoint::_fA)
  .def_readonly("fB",&ContactGenerator::ContactPoint::_fB)
  .def_readonly("tA2B",&ContactGenerator::ContactPoint::_tA2B);
  //ContactGenerator::ContactManifold
  py::class_<ContactGenerator::ContactManifold>(m,"ContactManifold")
  .def(py::init<>())
  .def("swap",&ContactGenerator::ContactManifold::swap)
  .def_readonly("points",&ContactGenerator::ContactManifold::_points)
  .def_readonly("sA",&ContactGenerator::ContactManifold::_sA)
  .def_readonly("sB",&ContactGenerator::ContactManifold::_sB)
  .def_readonly("sidA",&ContactGenerator::ContactManifold::_sidA)
  .def_readonly("sidB",&ContactGenerator::ContactManifold::_sidB)
  .def_readonly("jidA",&ContactGenerator::ContactManifold::_jidA)
  .def_readonly("jidB",&ContactGenerator::ContactManifold::_jidB)
  .def_readonly("tA",&ContactGenerator::ContactManifold::_tA)
  .def_readonly("tB",&ContactGenerator::ContactManifold::_tB)
  .def_readonly("DNDX",&ContactGenerator::ContactManifold::_DNDX)
  .def_readonly("x",&ContactGenerator::ContactManifold::_x);
  //ContactGenerator
  py::class_<ContactGenerator,
  std::shared_ptr<ContactGenerator>>(m,"ContactGenerator")
  .def(py::init<std::shared_ptr<ArticulatedBody>,std::vector<std::shared_ptr<ShapeExact>>>())
  .def("generateManifolds",[](std::shared_ptr<ContactGenerator> cg,T x0,bool useCCD,Mat3XT t,int status) {
    std::vector<ContactGenerator::ContactManifold> manifolds;
    cg->generateManifolds(x0,useCCD,manifolds,t,status);
    return manifolds;
  })
  .def("getExclude",&ContactGenerator::getExclude)
  .def("updateBVH",&ContactGenerator::updateBVH)
  .def("getBB",&ContactGenerator::getBB)
  .def("epsDist",&ContactGenerator::epsDist)
  .def("epsDir",&ContactGenerator::epsDir);
}
void initializeSimulator(py::module& m) {
  //SoftJoint
  typedef PBDArticulatedGradientInfo<T> Info;
  py::class_<SoftJoint>(m,"SoftJoint")
  .def(py::init<>())
  .def("posA",[](const SoftJoint& joint,std::shared_ptr<Info> info) {
    return joint.posA(*info);
  })
  .def("posB",[](const SoftJoint& joint,std::shared_ptr<Info> info) {
    return joint.posB(*info);
  })
  .def("rotA",[](const SoftJoint& joint,std::shared_ptr<Info> info,int d) {
    return joint.rotA(*info,d);
  })
  .def("rotB",[](const SoftJoint& joint,std::shared_ptr<Info> info,int d) {
    return joint.rotB(*info,d);
  })
  .def("rotAL",[](const SoftJoint& joint,int d) {
    return joint.rotAL(d);
  })
  .def("rotBL",[](const SoftJoint& joint,int d) {
    return joint.rotBL(d);
  })
  .def_readwrite("linearLimit",&SoftJoint::_linearLimit)
  .def_readwrite("angularLimit",&SoftJoint::_angularLimit)
  .def_readwrite("transA",&SoftJoint::_transA)
  .def_readwrite("transB",&SoftJoint::_transB)
  .def_readwrite("jidA",&SoftJoint::_jidA)
  .def_readwrite("jidB",&SoftJoint::_jidB);
  //Simulator::PhysicsParameter
  py::class_<Simulator::PhysicsParameter>(m,"PhysicsParameter")
  .def(py::init<>())
  .def_readwrite("isKinematic",&Simulator::PhysicsParameter::_isKinematic)
  .def_readwrite("isDesign",&Simulator::PhysicsParameter::_isDesign)
  .def_readwrite("kp",&Simulator::PhysicsParameter::_kp)
  .def_readwrite("kd",&Simulator::PhysicsParameter::_kd)
  .def_readwrite("friction",&Simulator::PhysicsParameter::_friction)
  .def_readwrite("kc",&Simulator::PhysicsParameter::_kc)
  .def("setKinFunc",[](Simulator::PhysicsParameter& param,Simulator::DOFFunction func) {
    param._kin=func;
  })
  .def("setTarPFunc",[](Simulator::PhysicsParameter& param,Simulator::DOFFunction func) {
    param._tarP=func;
  })
  .def("setTarDFunc",[](Simulator::PhysicsParameter& param,Simulator::DOFFunction func) {
    param._tarD=func;
  });
  //Simulator
  py::class_<Simulator,
  std::shared_ptr<Simulator>>(m,"Simulator")
  .def("clearJoint",&Simulator::clearJoint)
  .def("addJoint",&Simulator::addJoint)
  .def("getJoints",static_cast<const std::vector<SoftJoint>&(Simulator::*)()const>(&Simulator::getJoints))
  .def("setJoints",[](std::shared_ptr<Simulator> sim,const std::vector<SoftJoint>& joints) {
    sim->getJoints()=joints;
  })
  .def("clearShape",&Simulator::clearShape)
  .def("addShape",&Simulator::addShape)
  .def("setArticulatedBody",&Simulator::setArticulatedBody)
  .def("setHeuristcGuessStiffness",&Simulator::setHeuristcGuessStiffness)
  .def("getHeuristcGuessStiffness",&Simulator::getHeuristcGuessStiffness)
  .def("getManifolds",&Simulator::getManifolds)
  .def("getShapes",&Simulator::getShapes)
  .def("getJointPhysicsParameter",[](std::shared_ptr<Simulator> sim,int id) {
    return sim->getJointPhysicsParameter(id);
  })
  .def("setJointPhysicsParameter",[](std::shared_ptr<Simulator> sim,int id,Simulator::PhysicsParameter& param) {
    sim->getJointPhysicsParameter(id)=param;
  })
  .def("getBody",&Simulator::getBody)
  .def("getContact",&Simulator::getContactSmartPtr)
  .def_property("gravity",[](std::shared_ptr<Simulator> sim) {
    return sim->getGravity();
  },[](std::shared_ptr<Simulator> sim,const Simulator::Vec3T& g) {
    sim->setGravity(g);
  })
  .def_property("design",[](std::shared_ptr<Simulator> sim) {
    return sim->getDesign();
  },[](std::shared_ptr<Simulator> sim,const Simulator::Vec3T& d) {
    sim->setDesign(d);
  })
  .def_property("time",[](std::shared_ptr<Simulator> sim) {
    return sim->getTime();
  },[](std::shared_ptr<Simulator> sim,Simulator::T t) {
    sim->setTime(t);
  })
  .def("setOutput",&Simulator::setOutput)
  .def("step",&Simulator::step)
  .def_property("pos",[](std::shared_ptr<Simulator> sim) {
    return sim->pos();
  },[](std::shared_ptr<Simulator> sim,const Simulator::Vec& pos) {
    sim->setPos(pos);
  })
  .def_property("vel",[](std::shared_ptr<Simulator> sim) {
    return sim->vel();
  },[](std::shared_ptr<Simulator> sim,const Simulator::Vec& vel) {
    sim->setVel(vel);
  })
  .def("detectCurrentContact",&Simulator::detectCurrentContact);
  //PBDSimulator
  py::class_<PBDSimulator,std::shared_ptr<PBDSimulator>,
  Simulator>(m,"PBDSimulator")
  .def(py::init<T>())
  .def("setJTJ",&PBDSimulator::setJTJ)
  .def("setCrossTerm",&PBDSimulator::setCrossTerm);
  //XPBDSimulator
  py::class_<XPBDSimulator,std::shared_ptr<XPBDSimulator>,
  PBDSimulator>(m,"XPBDSimulator")
  .def(py::init<T>());
  //ConvHullPBDSimulator
  py::class_<ConvHullPBDSimulator,std::shared_ptr<ConvHullPBDSimulator>,
  Simulator>(m,"ConvHullPBDSimulator")
  .def(py::init<const ConvHullPBDSimulator&>())
  .def(py::init<T>())
  .def("setCustomEnergy",[&](std::shared_ptr<ConvHullPBDSimulator> sim,CustomPBDEnergy<T>::PBDEnergyForwardPython forward,CustomPBDEnergy<T>::PBDEnergyBackwardPython backward) {
    std::shared_ptr<CustomPBDEnergy<T>> energy(new CustomPBDEnergy<T>());
    energy->setForward(forward);
    energy->setBackward(backward);
    sim->setCustomEnergy(energy);
  })
  .def("getCentre",&ConvHullPBDSimulator::getCentre)
  .def("getdPTarget",&ConvHullPBDSimulator::getdPTarget)
  .def("getdDTarget",&ConvHullPBDSimulator::getdDTarget)
  .def("getdL",&ConvHullPBDSimulator::getdL)
  .def("getdLL",&ConvHullPBDSimulator::getdLL)
  .def("getdPos",&ConvHullPBDSimulator::getdPos)
  .def("getdD",&ConvHullPBDSimulator::getdD)
  .def("getConvexPoints",&ConvHullPBDSimulator::getConvexPoints)
  .def("getGlobalPoints",&ConvHullPBDSimulator::getGlobalPoints)
  .def("getJointPosGrad",[](std::shared_ptr<ConvHullPBDSimulator> sim,int k) {
    ConvHullPBDSimulator::GradInfo grad;
    sim->getJointPosGrad(grad,k);
    return grad;
  })
  .def("reset",&ConvHullPBDSimulator::reset)
  .def("resetWithPos",&ConvHullPBDSimulator::resetWithPos)
  .def("stepwithbackward",&ConvHullPBDSimulator::stepwithbackward)
  .def("backward",[](std::shared_ptr<ConvHullPBDSimulator> sim,bool debug) {
    ConvHullPBDSimulator::GradInfo grad;
    sim->backward(grad,debug);
    return grad;
  })
  .def("backward",static_cast<void(ConvHullPBDSimulator::*)()>(&ConvHullPBDSimulator::backward))
  .def_property("gTol",&ConvHullPBDSimulator::gTol,&ConvHullPBDSimulator::setGTol)
  .def_property("x0",&ConvHullPBDSimulator::x0,&ConvHullPBDSimulator::setX0)
  .def_property("hardLimit",&ConvHullPBDSimulator::hardLimit,&ConvHullPBDSimulator::setHardLimit)
  .def("setCoefBarrier",&ConvHullPBDSimulator::setCoefBarrier)
  .def("setIsDesign",&ConvHullPBDSimulator::setIsDesign)
  .def("setPD",&ConvHullPBDSimulator::setPD)
  .def("setclass",&ConvHullPBDSimulator::setclass)
  .def("getGlobalPoints",&ConvHullPBDSimulator::getGlobalPoints)
  .def("updateConvexPoints",&ConvHullPBDSimulator::updateConvexPoints);
  //MeshBasedPBDSimulator
  py::class_<MeshBasedPBDSimulator,std::shared_ptr<MeshBasedPBDSimulator>,
  ConvHullPBDSimulator>(m,"MeshBasedPBDSimulator")
  .def(py::init<const MeshBasedPBDSimulator&>())
  .def(py::init<T>());
  //CollisionGradInfo
  py::class_<CollisionGradInfo<T>>(m,"CollisionGradInfo").def(py::init<>())
  .def(py::init([](std::shared_ptr<ArticulatedBody> body,const Vec& theta) {
    return CollisionGradInfo<T>(*body,theta);
  }))
  .def("x",[](const CollisionGradInfo<T>& info)->Vec {
    return info._info._xM;
  })
  .def("reset",&CollisionGradInfo<T>::reset)
  .def("getTransformation",[&](const CollisionGradInfo<T>& info,std::shared_ptr<ArticulatedBody> body,int jid)->std::pair<Mat3T,Vec3T> {
    Mat3T R;
    Vec3T t;
    info.getTransformation(*body,jid,R,t);
    return std::make_pair(R,t);
  })
  .def("getJacobian",[&](const CollisionGradInfo<T>& info,std::shared_ptr<ArticulatedBody> body,int jid)->std::pair<Mat3XT,Mat3XT> {
    Mat3XT JV;
    Mat3XT JW;
    info.getJacobian(*body,jid,JV,JW);
    return std::make_pair(JV,JW);
  })
  .def("getJacobianDeriv",[&](const CollisionGradInfo<T>& info,std::shared_ptr<ArticulatedBody> body,int jid)->std::pair<std::vector<Mat3XT>,std::vector<Mat3XT>> {
    std::vector<Mat3XT> JV;
    std::vector<Mat3XT> JW;
    info.getJacobianDeriv(*body,jid,JV,JW);
    return std::make_pair(JV,JW);
  })
  .def_readonly("nrVex",&CollisionGradInfo<T>::_nrVex)
  .def_readonly("centre",&CollisionGradInfo<T>::_centre)
  .def_readonly("HTheta",&CollisionGradInfo<T>::_HTheta)
  .def_readonly("HThetaL",&CollisionGradInfo<T>::_HThetaL)
  .def_readonly("HThetaLL",&CollisionGradInfo<T>::_HThetaLL)
  .def_readonly("HThetaD",&CollisionGradInfo<T>::_HThetaD)
  .def_readonly("HThetaPTarget",&CollisionGradInfo<T>::_HThetaPTarget)
  .def_readonly("HThetaDTarget",&CollisionGradInfo<T>::_HThetaDTarget)
  .def_readonly("HPos",&CollisionGradInfo<T>::_HPos)
  .def_readonly("HThetaDesign",&CollisionGradInfo<T>::_HThetaDesign)
  .def_readonly("HThetaN",&CollisionGradInfo<T>::_HThetaN)
  .def_readonly("HThetaU",&CollisionGradInfo<T>::_HThetaU)
  .def_readonly("DTG",&CollisionGradInfo<T>::_DTG);
}
PYBIND11_MODULE(pyPBAD, m) {
  initializeJoint(m);
  initializeArticulatedBody(m);
  initializeArticulatedUtils(m);
  initializeArticulatedLoader(m);
  initializePBDGradientInfo(m);
  initializeShapes(m);
  initializeContact(m);
  initializeSimulator(m);
}*/
