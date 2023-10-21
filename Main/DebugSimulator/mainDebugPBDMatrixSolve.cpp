#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Articulated/NEArticulatedGradientInfo.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Simulator/PBDMatrixSolver.h>
#include <Utils/CrossSpatialUtils.h>

using namespace PHYSICSMOTION;

int main(int argc,char** argv) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  DEFINE_NUMERIC_DELTA_T(T)
  typedef NEArticulatedGradientInfo<T> GNE;
  typedef PBDArticulatedGradientInfo<T> GPBD;
  //create body
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),0,10,0.5f,0.1f,D2R(10),0,0,0,3,0,0,0);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  ArticulatedUtils utils(*body);
  utils.assemble(*(pt.RootElement()));
  {
    //inverting mass matrix
    Vec q=Vec::Random(body->nrDOF());
    GNE infoNE(*body,q);
    GPBD infoPBD(*body,q);
    Mat3XT MRR,MRt,MtR,Mtt,G;
    Mat3XT MRR2,MRt2,MtR2,Mtt2;
    MRR.setZero(3,body->nrJ()*3);
    MRt.setZero(3,body->nrJ()*3);
    MtR.setZero(3,body->nrJ()*3);
    Mtt.setZero(3,body->nrJ()*3);
    G.setZero(3,body->nrJ()*4);
    for(int k=0; k<body->nrJ(); k++) {
      const Joint& J=body->joint(k);
      MRR.template block<3,3>(0,k*3)=-invDoubleCrossMatTrace<T>(ROTI(infoPBD._TM,k)*J._MCCT.template cast<T>()*ROTI(infoPBD._TM,k).transpose());
      MRt.template block<3,3>(0,k*3)=cross<T>(ROTI(infoPBD._TM,k)*J._MC.template cast<T>());
      MtR.template block<3,3>(0,k*3)=-cross<T>(ROTI(infoPBD._TM,k)*J._MC.template cast<T>());
      Mtt.template block<3,3>(0,k*3)=Mat3T::Identity()*J._M;
      //compare 6x6 mass matrix
      Mat6T R=Mat6T::Zero(),MPBD;
      R.template block<3,3>(0,0)=R.template block<3,3>(3,3)=ROTI(infoPBD._TM,k);
      Mat6T MNE=infoNE._IM.template block<6,6>(0,k*6);
      MPBD.template block<3,3>(0,0)=MRR.template block<3,3>(0,k*3);
      MPBD.template block<3,3>(0,3)=MRt.template block<3,3>(0,k*3);
      MPBD.template block<3,3>(3,0)=MtR.template block<3,3>(0,k*3);
      MPBD.template block<3,3>(3,3)=Mtt.template block<3,3>(0,k*3);
      MPBD=R.transpose()*MPBD*R;
      DEBUG_GRADIENT("Joint"+std::to_string(k),MPBD.norm(),(MPBD-MNE).norm())
    }
    //compare the entire mass matrix
    MatT H;
    H.setZero(body->nrDOF(),body->nrDOF());
    infoPBD.toolABZ(*body,GPBD::mapM(MRR2=MRR),GPBD::mapM(MRt2=MRt),GPBD::mapM(MtR2=MtR),GPBD::mapM(Mtt2=Mtt),GPBD::mapM(G),GPBD::mapM(H));
    infoNE.calcH(*body);
    DEBUG_GRADIENT("Mass-Matrix",H.norm(),(H-infoNE._HM).norm())
    //invert using Eigen
    MatT b=Vec::Random(body->nrDOF());
    Vec diag=Vec::Random(body->nrDOF())+Vec::Ones(body->nrDOF());
    PBDMatrixSolverEigen solEig(body);
    solEig.compute(MatT(infoNE._HM+MatT(diag.asDiagonal())));
    Vec invHbEig=solEig.solve(b);
    //invert using CRBA
    PBDMatrixSolverCRBA solCRBA(body);
    solCRBA.compute(MatT(infoNE._HM+MatT(diag.asDiagonal())));
    DEBUG_GRADIENT("CRBA",invHbEig.norm(),(invHbEig-solCRBA.solve(b)).norm())
    //invert using ABA
    PBDMatrixSolverABA solABA(body);
    solABA.compute(q,Mat6XT(infoNE._IM),diag);
    DEBUG_GRADIENT("ABA",invHbEig.norm(),(invHbEig-solABA.solve(b)).norm())
    solABA.compute(q,MRR,MRt,MtR,Mtt,diag);
    DEBUG_GRADIENT("ABA2",invHbEig.norm(),(invHbEig-solABA.solve(b)).norm())
  }
  {
    //inverting arbitrary matrix
    Vec q=Vec::Random(body->nrDOF());
    GNE infoNE(*body,q);
    GPBD infoPBD(*body,q);
    Mat3XT MRR,MRt,MtR,Mtt,G;
    Mat3XT MRR2,MRt2,MtR2,Mtt2;
    MRR.setZero(3,body->nrJ()*3);
    MRt.setZero(3,body->nrJ()*3);
    MtR.setZero(3,body->nrJ()*3);
    Mtt.setZero(3,body->nrJ()*3);
    G.setZero(3,body->nrJ()*4);
    Mat6XT I=Mat6XT::Zero(6,body->nrJ()*6);
    for(int k=0; k<body->nrJ(); k++) {
      Mat6T R=Mat6T::Zero(),MPBD=Mat6T::Random();
      MPBD=(MPBD*MPBD.transpose()).eval();
      R.template block<3,3>(0,0)=R.template block<3,3>(3,3)=ROTI(infoPBD._TM,k);
      MRR.template block<3,3>(0,k*3)=MPBD.template block<3,3>(0,0);
      MRt.template block<3,3>(0,k*3)=MPBD.template block<3,3>(0,3);
      MtR.template block<3,3>(0,k*3)=MPBD.template block<3,3>(3,0);
      Mtt.template block<3,3>(0,k*3)=MPBD.template block<3,3>(3,3);
      I.template block<6,6>(0,k*6)=R.transpose()*MPBD*R;
    }
    //compare the entire matrix
    MatT H;
    H.setZero(body->nrDOF(),body->nrDOF());
    infoPBD.toolABZ(*body,GPBD::mapM(MRR2=MRR),GPBD::mapM(MRt2=MRt),GPBD::mapM(MtR2=MtR),GPBD::mapM(Mtt2=Mtt),GPBD::mapM(G),GPBD::mapM(H));
    infoNE.calcHInner(*body,infoNE._HM,GNE::mapCM(I));
    DEBUG_GRADIENT("General-Matrix",H.norm(),(H-infoNE._HM).norm())
    //invert using Eigen
    MatT b=Vec::Random(body->nrDOF());
    Vec diag=Vec::Random(body->nrDOF())+Vec::Ones(body->nrDOF());
    PBDMatrixSolverEigen solEig(body);
    solEig.compute(MatT(infoNE._HM+MatT(diag.asDiagonal())));
    Vec invHbEig=solEig.solve(b);
    //invert using CRBA
    PBDMatrixSolverCRBA solCRBA(body);
    solCRBA.compute(MatT(infoNE._HM+MatT(diag.asDiagonal())));
    DEBUG_GRADIENT("CRBA",invHbEig.norm(),(invHbEig-solCRBA.solve(b)).norm())
    //invert using ABA
    PBDMatrixSolverABA solABA(body);
    solABA.compute(q,I,diag);
    DEBUG_GRADIENT("ABA",invHbEig.norm(),(invHbEig-solABA.solve(b)).norm())
    solABA.compute(q,MRR,MRt,MtR,Mtt,diag);
    DEBUG_GRADIENT("ABA2",invHbEig.norm(),(invHbEig-solABA.solve(b)).norm())
  }
  return 0;
}
