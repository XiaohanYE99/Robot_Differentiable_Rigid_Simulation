#ifndef ARTICULATED_VISUALIZER_H
#define ARTICULATED_VISUALIZER_H

#include "ArticulatedBody.h"
#include <TinyVisualizer/Drawer.h>

namespace PHYSICSMOTION {
extern std::shared_ptr<DRAWER::Shape> visualizeArticulated
(std::shared_ptr<ArticulatedBody> b,
 const Eigen::Matrix<GLfloat,3,1>& colorBody=Eigen::Matrix<GLfloat,3,1>(.3f,.3f,.3f),
 const Eigen::Matrix<GLfloat,3,1>& colorBB=Eigen::Matrix<GLfloat,3,1>(.3f,.2f,.1f),
 const Eigen::Matrix<GLfloat,3,1>& colorBBC=Eigen::Matrix<GLfloat,3,1>(0.f,0.f,0.f),
 bool wire=false,bool convex=false);
extern std::shared_ptr<DRAWER::Shape> visualizeArticulated
(std::shared_ptr<ArticulatedBody> b,
 const std::vector<Eigen::Matrix<GLfloat,3,1>>& colorBody,
 const Eigen::Matrix<GLfloat,3,1>& colorBB=Eigen::Matrix<GLfloat,3,1>(.3f,.2f,.1f),
 const Eigen::Matrix<GLfloat,3,1>& colorBBC=Eigen::Matrix<GLfloat,3,1>(0.f,0.f,0.f),
 bool wire=false,bool convex=false);
extern void highlightJoint
(std::shared_ptr<DRAWER::Shape> root,std::shared_ptr<ArticulatedBody> b,
 const std::vector<char>* selected,
 const std::vector<char>* excluded,
 const std::vector<char>* markFoot,
 const Eigen::Matrix<GLfloat,3,1>& colorBody=Eigen::Matrix<GLfloat,3,1>(.3f,.3f,.3f));
extern int getRayIntersectedJoint(std::shared_ptr<DRAWER::Shape> root,std::shared_ptr<ArticulatedBody> b,const Eigen::Matrix<GLfloat,6,1>& ray);
extern void updateArticulatedBody(std::shared_ptr<DRAWER::Shape> shape,std::shared_ptr<ArticulatedBody> b,const ArticulatedBody::Mat3XT& M);
extern void updateArticulatedBody(std::shared_ptr<DRAWER::Shape> shape,std::shared_ptr<ArticulatedBody> b,const ArticulatedBody::Vec& x);
extern void compareVisualMeshURDF(int argc,char** argv,const std::string& file,bool wire=false,
                                  const Eigen::Matrix<GLfloat,3,1>& color=Eigen::Matrix<GLfloat,3,1>(.7f,.7f,.7f));
}

#endif
