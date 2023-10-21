#ifndef GJK_POLYTOPE_H
#define GJK_POLYTOPE_H

#include "EntityId.h"

namespace PHYSICSMOTION {
template <typename T>
class GJKPolytope : public SerializableBase {
 public:
  DECL_MAT_VEC_MAP_TYPES_T
  typedef double Point[3];
  using SerializableBase::read;
  using SerializableBase::write;
  GJKPolytope();
  GJKPolytope(int JID,const ArticulatedBody& body,const CollisionGradInfo<T>& info);
  GJKPolytope(std::shared_ptr<MeshExact> mesh);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  int jid() const;
  std::shared_ptr<MeshExact> mesh() const;
  void writeVTK(const std::string& path,const ArticulatedBody* body) const;
  void writeVTK(VTKWriter<T>& os,const ArticulatedBody* body) const;
  static void writeConfigVTK(const std::string& path,
                             const GJKPolytope<T>& p1,const GJKPolytope<T>& p2,
                             GJKPolytope<T>::Point wp1,GJKPolytope<T>::Point wp2,
                             const ArticulatedBody* body);
  static void writeConfig(const std::string& path,
                          const GJKPolytope<T>& p1,const GJKPolytope<T>& p2);
  static void readConfig(const std::string& path);
  static Vec2T project(const GJKPolytope<T>& p,const Vec3T& n);
  static T distance(const GJKPolytope<T>& A,const GJKPolytope<T>& B,Point pA,Point pB);
 private:
  int _jid,_n;
  Mat3X4T _trans;
  std::shared_ptr<MeshExact> _mesh;
  std::vector<double> _coord;
};
}
#endif
