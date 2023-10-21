#include <Environment/TriangleExact.h>
#include <Environment/MeshExact.h>
#include <Environment/BBoxExact.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/SphericalBBoxExact.h>
#include <Environment/CompositeShapeExact.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Environment/Environment.h>
#include <Utils/DebugGradient.h>
#include <Utils/RotationUtils.h>
#include <TinyVisualizer/MakeTexture.h>
#include <random>

using namespace PHYSICSMOTION;

//#define DEBUG_TRIANGLE
//#define DEBUG_MESH
//#define DEBUG_CONVEX
//#define DEBUG_BOX
//#define DEBUG_SBOX
#define DEBUG_COMPOSITE
#define DEBUG_ENV
template <typename T>
void debugIO(std::shared_ptr<T>& val) {
  val->SerializableBase::writeStr("dat");
  val.reset(new T());
  val->SerializableBase::readStr("dat");
}
void debugTrianglePointDist(int argc,char** argv) {
  using namespace DRAWER;
  Drawer drawer(argc,argv);
#define DIST 10
  drawer.setKeyFunc([&](GLFWwindow* wnd,int key,int scan,int action,int mods,bool captured) {
    typedef TriangleExact::T T;
    typedef TriangleExact::Vec3T Vec3T;
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS) {
      drawer.clear();
      std::shared_ptr<TriangleExact> et(new TriangleExact(Eigen::Matrix<double,3,1>::Random().template cast<T>(),
                                        Eigen::Matrix<double,3,1>::Random().template cast<T>(),
                                        Eigen::Matrix<double,3,1>::Random().template cast<T>()));
      std::shared_ptr<Shape> tri=visualizeTriangleExact(et);
      tri->setColor(GL_TRIANGLES,.7,.7,.7);
      drawer.addShape(tri);

      T sqrDist;
      Eigen::Matrix<int,2,1> feat;
      std::shared_ptr<MeshShape> lines(new MeshShape);
      Vec3T l0=Eigen::Matrix<double,3,1>::Random().cast<T>()*DIST;
      Vec3T l1=Eigen::Matrix<double,3,1>::Random().cast<T>()*DIST,b;
      for(int i=0; i<DIST*10; i++) {
        T alpha=T(i)/T(DIST*10);
        Vec3T pt=l0*(1-alpha)+l1*alpha,cp;
        et->calcPointDist(pt,sqrDist,cp,b,feat);
        lines->addVertex(pt.cast<GLfloat>());
        lines->addVertex(cp.cast<GLfloat>());
        lines->addIndex(Eigen::Matrix<GLuint,2,1>(i*2+0,i*2+1));
      }
      lines->setMode(GL_LINES);
      lines->setColor(GL_LINES,.7,0,0);
      drawer.addShape(lines);
    }
  });
#undef DIST
  drawer.addCamera3D(90);
  drawer.mainLoop();
}
template <typename T2>
void debugMeshPointDist(int argc,char** argv,std::shared_ptr<ShapeExact> m,bool useGetMesh=false,int nrIter=100) {
  using namespace DRAWER;
  Drawer drawer(argc,argv);
#define DIST 10
  drawer.setKeyFunc([&](GLFWwindow*,int key,int,int action,int,bool captured) {
    typedef Eigen::Matrix<T2,3,1> Vec3T;
    typedef Eigen::Matrix<T2,3,3> Mat3T;
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS) {
      drawer.clear();
      std::shared_ptr<Shape> mShape;
      if(useGetMesh)
        mShape=visualizeShapeExactGetMesh(m,true);
      else mShape=visualizeShapeExact(m,true);
      mShape->setColor(GL_LINES,.7,.7,.7);
      drawer.addShape(mShape);

      DEFINE_NUMERIC_DELTA_T(T2)
      int nIn=0,nOut=0;
      Mat3T hessian,hessian2;
      Vec3T n,normal,n2,normal2;
      Eigen::Matrix<int,2,1> feat,feat2;
      std::shared_ptr<MeshShape> linesIn(new MeshShape);
      std::shared_ptr<MeshShape> linesOut(new MeshShape);
      for(int i=0; i<nrIter; i++) {
        Vec3T pt=Vec3T::Random()*0.5+Vec3T::Constant(0.5);
        BBoxExact bb=m->getBB().enlargedEps(1);
        pt.array()=(bb.minCorner().array()*(1-pt.template cast<GEOMETRY_SCALAR>().array())+bb.maxCorner().array()*pt.template cast<GEOMETRY_SCALAR>().array()).template cast<T2>();
        T2 dist=m->closest<T2>(pt,n,normal,hessian,feat);
        std::string type=dist<0?"I":"O";
        type+=feat[0]<0?"T":feat[1]<0?"V":"E";
        for(int d=0; d<3; d++) {
          std::string name="N"+type+std::to_string(d);
          T2 dist2=m->closest<T2>(pt+Vec3T::Unit(d)*DELTA,n2,normal2,hessian2,feat2);
          DEBUG_GRADIENT(name.c_str(),normal[d],normal[d]-(dist2-dist)/DELTA)
        }
        if(!hessian.isZero())
          for(int d=0; d<3; d++) {
            std::string name="H"+type+std::to_string(d);
            m->closest<T2>(pt+Vec3T::Unit(d)*DELTA,n2,normal2,hessian2,feat2);
            DEBUG_GRADIENT(name.c_str(),(hessian*Vec3T::Unit(d)).norm(),(hessian*Vec3T::Unit(d)-(normal2-normal)/DELTA).norm())
          }
        if(dist<0) {
          linesIn->addVertex(pt.template cast<GLfloat>());
          linesIn->addVertex((pt+n).template cast<GLfloat>());
          linesIn->addIndex(Eigen::Matrix<GLuint,2,1>(nIn*2+0,nIn*2+1));
          nIn++;
        } else {
          linesOut->addVertex(pt.template cast<GLfloat>());
          linesOut->addVertex((pt+n).template cast<GLfloat>());
          linesOut->addIndex(Eigen::Matrix<GLuint,2,1>(nOut*2+0,nOut*2+1));
          nOut++;
        }
      }
      linesIn->setMode(GL_LINES);
      linesIn->setColor(GL_LINES,.7,0,0);
      if(nIn>0)
        drawer.addShape(linesIn);
      linesOut->setMode(GL_LINES);
      linesOut->setColor(GL_LINES,0,0,.7);
      if(nOut>0)
        drawer.addShape(linesOut);
    }
  });
#undef DIST
  drawer.addCamera3D(90);
  drawer.mainLoop();
}
template <typename T>
void debugEnvironment(int argc,char** argv,std::shared_ptr<Environment<T>> e,int nrIter=100) {
  using namespace DRAWER;
  Drawer drawer(argc,argv);
#define DIST 10
  std::shared_ptr<Texture> tex=drawGrid();
  drawer.setKeyFunc([&](GLFWwindow*,int key,int,int action,int,bool captured) {
    typedef typename Environment<T>::Vec3T Vec3T;
    typedef typename Environment<T>::Mat3T Mat3T;
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS) {
      drawer.clear();
      std::shared_ptr<Shape> mShape=visualizeEnvironment(e,Eigen::Matrix<GLfloat,2,1>(.25f,.25f));
      mShape->setTexture(tex);
      drawer.addShape(mShape);

      Mat3T h;
      Vec3T pt,g,g2;
      BBoxExact bb=e->getBB();
      DEFINE_NUMERIC_DELTA_T(T)
      for(int i=0; i<nrIter; i++) {
        for(int r=0; r<3; r++) {
          std::random_device rd;
          std::mt19937 gen(rd());
          std::uniform_real_distribution<> dis(0,1);
          Eigen::Matrix<T,3,1> alpha(dis(gen),dis(gen),dis(gen));
          pt[r]=(T)bb.minCorner()[r]*(1-alpha[r])+(T)bb.maxCorner()[r]*alpha[r];
        }
        Vec3T delta=Vec3T::Random();
        T p=e->phi(pt,&g),p2=e->phi(pt+delta*DELTA);
        if(g.dot(delta)!=0) {
          DEBUG_GRADIENT("phiGrad",g.dot(delta),g.dot(delta)-(p2-p)/DELTA)
        }
        g=e->phiGrad(pt,&h),g2=e->phiGrad(pt+delta*DELTA);
        if(sqrt((h*delta).squaredNorm())!=0) {
          DEBUG_GRADIENT("phiHess",sqrt((h*delta).squaredNorm()),sqrt((h*delta-(g2-g)/DELTA).squaredNorm()))
        }
      }
    }
  });
#undef DIST
  drawer.addCamera3D(90);
  drawer.mainLoop();
}
int main(int argc,char** argv) {
  typedef FLOAT T2;
#ifdef DEBUG_TRIANGLE
  debugTrianglePointDist(argc,argv);
#endif
#ifdef DEBUG_MESH
  MeshExact m("bunny.obj");
  std::shared_ptr<MeshExact> mCopy=std::dynamic_pointer_cast<MeshExact>(m.copy());
  debugIO(mCopy);
  debugMeshPointDist<T2>(argc,argv,mCopy,true);
#endif
#ifdef DEBUG_CONVEX
  MeshExact m("bunny.obj");
  ConvexHullExact mc(m);
  std::shared_ptr<ConvexHullExact> mcCopy=std::dynamic_pointer_cast<ConvexHullExact>(mc.copy());
  debugIO(mcCopy);
  debugMeshPointDist<T2>(argc,argv,mcCopy,true);
#endif
#ifdef DEBUG_BOX
  Eigen::Matrix<double,3,1> minC(1,2,3);
  Eigen::Matrix<double,3,1> maxC(4,9,7);
  BBoxExact bb(minC.template cast<GEOMETRY_SCALAR>(),
               maxC.template cast<GEOMETRY_SCALAR>());
  std::shared_ptr<BBoxExact> bbCopy=std::dynamic_pointer_cast<BBoxExact>(bb.copy());
  debugIO(bbCopy);
  debugMeshPointDist<T2>(argc,argv,bbCopy,true);
#endif
#ifdef DEBUG_SBOX
  SphericalBBoxExact sbb(0.2f,0.3f,0.4f,0.5f);
  std::shared_ptr<SphericalBBoxExact> sbbCopy=std::dynamic_pointer_cast<SphericalBBoxExact>(sbb.copy());
  debugIO(sbbCopy);
  debugMeshPointDist<T2>(argc,argv,sbbCopy,true);
#endif
#ifdef DEBUG_COMPOSITE
  CompositeShapeExact::Mat3X4T t;
  std::vector<std::shared_ptr<ShapeExact>> geoms;
  std::vector<CompositeShapeExact::Mat3X4T,Eigen::aligned_allocator<CompositeShapeExact::Mat3X4T>> trans;
  geoms.push_back(std::shared_ptr<ShapeExact>(new SphericalBBoxExact(0.2f,0.3f,0.4f,0.5f)));
  t.block<3,3>(0,0)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<GEOMETRY_SCALAR>();
  t.block<3,1>(0,3).setRandom();
  trans.push_back(t);
  geoms.push_back(std::shared_ptr<ShapeExact>(new SphericalBBoxExact(0.1f,0.1f,0.1f,0.5f)));
  t.block<3,3>(0,0)=eulerX1Y3Z2<double,Eigen::Matrix<double,3,1>>(Eigen::Matrix<double,3,1>::Random()*M_PI,NULL,NULL).template cast<GEOMETRY_SCALAR>();
  t.block<3,1>(0,3).setRandom();
  trans.push_back(t);
  CompositeShapeExact c(geoms,trans);
  std::shared_ptr<CompositeShapeExact> cCopy=std::dynamic_pointer_cast<CompositeShapeExact>(c.copy());
  debugIO(cCopy);
  debugMeshPointDist<T2>(argc,argv,cCopy,true);
#endif
#ifdef DEBUG_ENV
  EnvironmentHeight<T2> e;
  e.createStair(1,2,3,4,0.2,5);
  std::shared_ptr<EnvironmentHeight<T2>> eCopy=std::dynamic_pointer_cast<EnvironmentHeight<T2>>(e.copy());
  debugIO(eCopy);
  debugEnvironment<T2>(argc,argv,eCopy);
#endif
  return 0;
}
