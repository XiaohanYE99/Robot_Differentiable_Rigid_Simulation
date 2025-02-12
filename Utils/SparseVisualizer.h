#ifndef SPARSE_VISUALIZER_H
#define SPARSE_VISUALIZER_H

#include <Eigen/Sparse>
#include <TinyVisualizer/Drawer.h>
#include <TinyVisualizer/MakeMesh.h>
#include <TinyVisualizer/Box2DShape.h>

namespace PHYSICSMOTION {
template <typename T,int O,typename I>
void visualizeSparse(int argc,char** argv,const Eigen::SparseMatrix<T,O,I>& a) {
  using namespace DRAWER;
  Drawer drawer(argc,argv);
  //frame
  std::shared_ptr<MeshShape> frame=makeSquare(false,Eigen::Matrix<GLfloat,2,1>(a.cols()/2.,a.rows()/2.));
  frame->setColorDiffuse(GL_LINE_LOOP,.7,.7,.7);
  std::shared_ptr<Box2DShape> movedFrame(new Box2DShape);
  movedFrame->setLocalTransform(a.cols()/2.,a.rows()/2.,1);
  movedFrame->addShape(frame);
  drawer.addShape(movedFrame);
  //entries
  std::shared_ptr<MeshShape> entry=makeSquare(true,Eigen::Matrix<GLfloat,2,1>(.5f,.5f));
  entry->setColorDiffuse(GL_QUADS,.7,.7,.7);
  std::shared_ptr<MeshShape> entryFrame=makeSquare(false,Eigen::Matrix<GLfloat,2,1>(.5f,.5f));
  entryFrame->setColorDiffuse(GL_LINE_LOOP,.5,.5,.5);
  for(int k=0; k<a.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(a,k); it; ++it) {
      std::shared_ptr<Box2DShape> movedEntry(new Box2DShape);
      movedEntry->setLocalTransform(it.col()+.5f,a.rows()-it.row()-.5f,1);
      movedEntry->addShape(entryFrame);
      movedEntry->addShape(entry);
      drawer.addShape(movedEntry);
    }
  //kick off app
  drawer.addCamera2D(a.rows()+a.cols());
  drawer.mainLoop();
}
}

#endif
