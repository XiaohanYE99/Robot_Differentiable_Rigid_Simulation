#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Simulator/SimulatorVisualizer.h>
#include <Simulator/ConvHullPBDSimulator.h>
#include <Simulator/MeshBasedPBDSimulator.h>
#include <chrono>
using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  DECL_MAT_VEC_MAP_TYPES_T
  std::vector<Eigen::Matrix<double, 3, 1>> shape;
  shape.push_back(Eigen::Matrix<double, 3, 1>(.3,.3,.3));
  shape.push_back(Eigen::Matrix<double, 3, 1>(.1,.1,.8));
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(ArticulatedLoader::createPushTask(shape,true)));
  //ArticulatedLoader::createChain(*(pt.RootElement()),Joint::TRANS_2D|Joint::ROT_3D_EXP,3,0.8f,.1f,D2R(90),D2R(0),D2R(0),0,3,0,0,0);
  std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-30,-30,-200),BBoxExact::Vec3T(20,20,-.08)));
  //std::shared_ptr<ShapeExact> floor(new BBoxExact(BBoxExact::Vec3T(-15,0.62,-5),BBoxExact::Vec3T(5,5,5)));
  //simulator
  //ConvHullPBDSimulator sim(0.01f);//0.01
  MeshBasedPBDSimulator sim(0.02f);
  //PBDSimulator sim(0.01f);
  sim.setfri(.4);
  sim.setOutput(true);
  sim.setArticulatedBody(body);
  sim.addShape(floor);
  sim.setGravity(ConvHullPBDSimulator::Vec3T(0,0,-9.81));

    Vec x;
    x.setZero(12);
    /*x[0]=1.49672087e+00;
    x[1]=5.03350648e-01;
    x[2]=-1.47205678e-01;
    x[3]=1.58870264+0;
    x[4]=9.14689821e-03;
    x[5]=1.90744384e-02;
    x[6]=1.74295909e+0;
    x[7]=1.67147605e-01;
    x[8]=1.69666133e-01;
    x[9]=1.57076642e+00;
    x[10]=4.02036761e-05;
    x[11]=3.94697528e-01;*/
    x[0]=1.5;
    x[1]=0.5;
    x[2]=0.08;
    x[3]=1.58870264+0;
    x[4]=9.14689821e-03;
    x[5]=1.90744384e-02;
    x[6]=2.2;
    x[7]=.6;
    x[8]=.25;
    x[9]=1.57079633;
  sim.resetWithPos(x);
  Vec v;
  v.setZero(12);
  v[6]=-6;
  sim.setVel(v);
  //PD controller
  /*for(int k=0; k<body->nrJ(); k++) {
    //control
    Joint& J=body->joint(k);
    J._control.setOnes(J.nrDOF());
    std::shared_ptr<MeshExact> local=std::dynamic_pointer_cast<MeshExact>(J._mesh);
    if(local) std::cout<<k<<" "<<J.nrDOF()<<" "<<local->iss().size()<<std::endl;
  }*/
  //run app
  visualizeSimulator(argc,argv,sim);
  /*auto start = std::chrono::high_resolution_clock::now();
  std::vector<Vec> traj;
  traj.resize(1);
  for(int i=0;i<1;i++) {
    std::cout<<i<<std::endl;
    traj.at(i)=sim.pos();
    //sim.step();
    
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  std::cout << "Elapsed time: " << elapsed_time.count() << " seconds" << std::endl;
  int frame=0;
  render.visualize(sim,&traj,&frame,1);*/
  return 0;
}
