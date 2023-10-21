#include "DistanceFunction.h"
#include <Utils/VTKWriter.h>

namespace PHYSICSMOTION {
template <typename T>
void writePointTetrahedronVTK(const std::string& path,
                              const Eigen::Matrix<T,3,1>& pt,
                              const Eigen::Matrix<T,3,1> v[4],
                              const Eigen::Matrix<T,3,1>& cp) {
  std::vector<Eigen::Matrix<double,3,1>> vss;
  std::vector<Eigen::Matrix<int,2,1>> iss;
  VTKWriter<double> os("PointToTetrahedron",path,true);
  //write lines
  vss.push_back(v[0].template cast<double>());
  vss.push_back(v[1].template cast<double>());
  vss.push_back(v[2].template cast<double>());
  vss.push_back(v[3].template cast<double>());
  vss.push_back(pt.template cast<double>());
  vss.push_back(cp.template cast<double>());
  iss.push_back(Eigen::Matrix<int,2,1>(0,1));
  iss.push_back(Eigen::Matrix<int,2,1>(0,2));
  iss.push_back(Eigen::Matrix<int,2,1>(0,3));
  iss.push_back(Eigen::Matrix<int,2,1>(1,2));
  iss.push_back(Eigen::Matrix<int,2,1>(1,3));
  iss.push_back(Eigen::Matrix<int,2,1>(2,3));
  iss.push_back(Eigen::Matrix<int,2,1>(4,5));
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::LINE);
  //write points
  iss.push_back(Eigen::Matrix<int,2,1>(0,-1));
  iss.push_back(Eigen::Matrix<int,2,1>(1,-1));
  iss.push_back(Eigen::Matrix<int,2,1>(2,-1));
  iss.push_back(Eigen::Matrix<int,2,1>(3,-1));
  iss.push_back(Eigen::Matrix<int,2,1>(4,-1));
  iss.push_back(Eigen::Matrix<int,2,1>(5,-1));
  os.appendCells(iss.begin(),iss.end(),VTKWriter<double>::POINT);
}
bool sameFeat(Eigen::Matrix<char,3,1> a,Eigen::Matrix<char,3,1> b) {
  sort3(a[0],a[1],a[2]);
  sort3(b[0],b[1],b[2]);
  return a==b;
}
template <typename T>
void debugDistToSqrTetrahedron() {
  DECL_MAT_VEC_MAP_TYPES_T
  using namespace std;
  Eigen::Matrix<char,3,1> featList[15]= {
    Eigen::Matrix<char,3,1>(-1,-1,-1),
    //vert
    Eigen::Matrix<char,3,1>(0,-1,-1),
    Eigen::Matrix<char,3,1>(1,-1,-1),
    Eigen::Matrix<char,3,1>(2,-1,-1),
    Eigen::Matrix<char,3,1>(3,-1,-1),
    //edge
    Eigen::Matrix<char,3,1>(0,1,-1),
    Eigen::Matrix<char,3,1>(0,2,-1),
    Eigen::Matrix<char,3,1>(0,3,-1),
    Eigen::Matrix<char,3,1>(1,2,-1),
    Eigen::Matrix<char,3,1>(1,3,-1),
    Eigen::Matrix<char,3,1>(2,3,-1),
    //face
    Eigen::Matrix<char,3,1>(0,1,2),
    Eigen::Matrix<char,3,1>(0,1,3),
    Eigen::Matrix<char,3,1>(0,2,3),
    Eigen::Matrix<char,3,1>(1,2,3),
  };
  for(int i=0; i<15; i++) {
    while(true) {
      Vec4T bary;
      Vec3T pt=Vec3T::Random(),cp;
      Vec3T v[4]= {
        Vec3T::Random(),
        Vec3T::Random(),
        Vec3T::Random(),
        Vec3T::Random()
      };
      Eigen::Matrix<char,3,1> feat;
      distToSqrTetrahedron(pt,v,Eigen::Map<Vec4T>(bary.data()),cp,&feat);
      if(sameFeat(feat,featList[i])) {
        writePointTetrahedronVTK("PointToTetCase"+std::to_string(i)+".vtk",pt,v,cp);
        break;
      }
    }
  }
}
//instance
template void debugDistToSqrTetrahedron<double>();
}
