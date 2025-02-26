#ifndef GJK_POLYTOPE_H
#define GJK_POLYTOPE_H

#include <Articulated/ArticulatedBody.h>

namespace PHYSICSMOTION {
template <typename T>
struct CollisionGradInfo;
template <typename T>
class GJKPolytope : public SerializableBase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef double Point[3];
  using SerializableBase::read;
  using SerializableBase::write;
  GJKPolytope();
  GJKPolytope(int JID,const ArticulatedBody& body,const CollisionGradInfo<T>& info);
  GJKPolytope(int JID,std::shared_ptr<MeshExact> mesh,const CollisionGradInfo<T>& info);
  GJKPolytope(std::shared_ptr<MeshExact> mesh);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  int jid() const;
  void setVertexId(const Eigen::Matrix<int,2,1>& off);
  Eigen::Matrix<int,2,1> getVertexId() const ;
  const Mat3XT& globalVss() const;
  const Mat3X4T& trans() const;
  std::shared_ptr<MeshExact> mesh() const;
  const std::vector<Node<int,BBoxExact>>& getBVH() const;
  void resetGlobalVss(const CollisionGradInfo<T>* info,bool buildBVH=true);
  void writeVTK(const std::string& path) const;
  void writeVTK(VTKWriter<double>& os) const;
  static void writeConfigVTK(const std::string& path,
                             const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,
                             GJKPolytope<T>::Point wp1,GJKPolytope<T>::Point wp2);
  static void writeConfig(const std::string& path,
                          const GJKPolytope<T>& p1,const GJKPolytope<T>& p2);
  static void readConfig(const std::string& path);
  static Vec2T project(const GJKPolytope<T>& p,const Vec3T& n);
  static T distance(const GJKPolytope<T>& A,const GJKPolytope<T>& B,Point pA,Point pB);
 private:
  int _jid,_jidP;
  int _offBeg,_offEnd;
  Mat3X4T _trans;
  Mat3XT _globalVss;
  std::shared_ptr<MeshExact> _mesh;
  std::vector<Node<int,BBoxExact>> _bvh;
};
}
#endif
