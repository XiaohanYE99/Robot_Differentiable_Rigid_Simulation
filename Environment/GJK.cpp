#include "GJK.h"
#include <Utils/Utils.h>
#include <Utils/Pragma.h>
#include <Utils/VTKWriter.h>
#include "DistanceFunction.h"

namespace PHYSICSMOTION {
template <int DIM>
void updateFeat(GJK::GJKPoint v[4],int& nrP,char* feat) {
  if(feat[0]==-1)
    return;
  else if(DIM>=3 && feat[2]>=0) {
    sort3(feat[0],feat[1],feat[2]);
    v[0]=v[(int)feat[0]];
    v[1]=v[(int)feat[1]];
    v[2]=v[(int)feat[2]];
    nrP=3;
    return;
  } else if(DIM>=2 && feat[1]>=0) {
    sort2(feat[0],feat[1]);
    v[0]=v[(int)feat[0]];
    v[1]=v[(int)feat[1]];
    nrP=2;
    return;
  } else if(DIM>=1 && feat[0]>=0) {
    v[0]=v[(int)feat[0]];
    nrP=1;
    return;
  }
  ASSERT_MSG(false,"Strange error in GJK, unknown feat status!")
}
void GJK::GJKPoint::calculate(const Mat3X4T& transA,const Mat3X4T& transB) {
  _ptAB =ROT(transA)*_ptAL+CTR(transA);
  _ptAB-=ROT(transB)*_ptBL+CTR(transB);
}
GJK::Vec3T GJK::computeD(const GJKPoint v[4],int nrP,T* bary,
                         const Mat3X4T& transA,
                         const Mat3X4T& transB,
                         Vec3T& pAL,Vec3T& pBL) {
  pAL=pBL=Vec3T::Zero();
  for(int d=0; d<nrP; d++) {
    pAL+=v[d]._ptAL*bary[d];
    pBL+=v[d]._ptBL*bary[d];
  }
  return (ROT(transA)*pAL+CTR(transA))-(ROT(transB)*pBL+CTR(transB));
}
GJK::T GJK::runGJK(std::shared_ptr<ShapeExact> A,
                   std::shared_ptr<ShapeExact> B,
                   const Mat3X4T& transA,
                   const Mat3X4T& transB,
                   Vec3T& pAL,Vec3T& pBL,bool* intersect) {
  int nrP;
  Vec3T cp,D;
  GJKPoint v[4];
  T dist,minDist;
  Eigen::Matrix<T,4,1> bary;
  //initialize
  nrP=1;
  v[0]._ptAL=A->support(Vec3T::Zero(),v[0]._idA);
  v[0]._ptBL=B->support(Vec3T::Zero(),v[0]._idB);
  v[0].calculate(transA,transB);
  dist=minDist=std::numeric_limits<double>::max();
  D=v[0]._ptAB;
  //main loop
  if(intersect)
    *intersect=false;
  while(true) {
    //insert new point
    v[nrP]._ptAL=A->support(-ROT(transA).transpose()*D,v[nrP]._idA);
    v[nrP]._ptBL=B->support( ROT(transB).transpose()*D,v[nrP]._idB);
    v[nrP++].calculate(transA,transB);
    //calculate distance
    if(nrP==4) {
      Eigen::Matrix<char,3,1> feat;
      Vec3T vAB[4]= {v[0]._ptAB,v[1]._ptAB,v[2]._ptAB,v[3]._ptAB};
      dist=distToSqrTetrahedron<T>(Vec3T::Zero(),vAB,Eigen::Map<Eigen::Matrix<T,4,1>>(bary.data()),cp,&feat);
      D=computeD(v,nrP,bary.data(),transA,transB,pAL,pBL);
      updateFeat<3>(v,nrP,feat.data());
      if(dist>=minDist)
        break;
      else {
        minDist=dist;
        if(feat[0]==-1) {
          if(intersect)
            *intersect=true;
          break;
        }
      }
    } else if(nrP==3) {
      Eigen::Matrix<char,2,1> feat;
      Vec3T vAB[3]= {v[0]._ptAB,v[1]._ptAB,v[2]._ptAB};
      dist=distToSqrTriangle<T>(Vec3T::Zero(),vAB,Eigen::Map<Eigen::Matrix<T,3,1>>(bary.data()),cp,&feat);
      D=computeD(v,nrP,bary.data(),transA,transB,pAL,pBL);
      updateFeat<2>(v,nrP,feat.data());
      if(dist>=minDist)
        break;
      else minDist=dist;
    } else if(nrP==2) {
      char feat;
      Vec3T vAB[2]= {v[0]._ptAB,v[1]._ptAB};
      dist=distToSqrLineSegment<T>(Vec3T::Zero(),vAB,Eigen::Map<Eigen::Matrix<T,2,1>>(bary.data()),cp,&feat);
      D=computeD(v,nrP,bary.data(),transA,transB,pAL,pBL);
      updateFeat<1>(v,nrP,&feat);
      if(dist>=minDist)
        break;
      else minDist=dist;
    } else {
      ASSERT_MSG(false,"Strange error in GJK, nrP==1!")
    }
  }
  return minDist;
}
void GJK::writeVTK(const std::string& path,
                   std::shared_ptr<MeshExact> A,
                   std::shared_ptr<MeshExact> B,
                   const Mat3X4T& transA,
                   const Mat3X4T& transB) {
  Vec3T pAL,pBL;
  bool intersect;
  runGJK(A,B,transA,transB,pAL,pBL,&intersect);
  VTKWriter<double> os("GJK",path,true);
  //write mesh
  A->writeVTK(os,transA);
  B->writeVTK(os,transB);
  //write distance
  runGJK(A,B,transA,transB,pAL,pBL,&intersect);
  //write GJK
  os.setRelativeIndex();
  std::vector<Eigen::Matrix<double,3,1>> vss;
  vss.push_back((ROT(transA)*pAL+CTR(transA)).template cast<double>());
  vss.push_back((ROT(transB)*pBL+CTR(transB)).template cast<double>());
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(0,0,1),
                 VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(1,0,1),
                 VTKWriter<double>::LINE,true);
  os.appendCells(VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(0,0,0),
                 VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(2,0,0),
                 VTKWriter<double>::POINT,true);
}
void GJK::writeVTK(const std::string& path,
                   std::shared_ptr<MeshExact> A,
                   std::shared_ptr<PointCloudExact> B,
                   const Mat3X4T& transA,
                   const Mat3X4T& transB) {
  Eigen::Matrix<double,3,1> pAL,pBL;
  runPD<double>(A,B,
                transA.template cast<double>(),
                transB.template cast<double>(),
                pAL,pBL);
  VTKWriter<double> os("PD",path,true);
  //write mesh
  A->writeVTK(os,transA);
  B->writeVTK(os,transB);
  //write GJK
  os.setRelativeIndex();
  std::vector<Eigen::Matrix<double,3,1>> vss;
  vss.push_back((ROT(transA)*pAL.template cast<T>()+CTR(transA)).template cast<double>());
  vss.push_back((ROT(transB)*pBL.template cast<T>()+CTR(transB)).template cast<double>());
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(0,0,1),
                 VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(1,0,1),
                 VTKWriter<double>::LINE,true);
  os.appendCells(VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(0,0,0),
                 VTKWriter<double>::IteratorIndex<Eigen::Matrix<int,2,1>>(2,0,0),
                 VTKWriter<double>::POINT,true);
}
}
