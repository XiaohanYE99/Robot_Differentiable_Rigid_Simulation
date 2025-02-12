#include <Simulator/Simulator.h>

#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/CameraExportPlugin.h>
#include <TinyVisualizer/CaptureGIFPlugin.h>
#include <TinyVisualizer/ImGuiPlugin.h>
#include <TinyVisualizer/ShadowAndLight.h>
#include <TinyVisualizer/Camera3D.h>
#include <TinyVisualizer/Povray.h>
#include <imgui/imgui.h>

namespace PHYSICSMOTION {
typedef FLOAT T;
DECL_MAT_VEC_MAP_TYPES_T
extern void visualizeSimulator(int argc,char** argv,Simulator& sim,int nSub=1,float maxStiffness=1e4f);

//void visualize(Simulator& sim,const std::vector<Vec>* theta,int* frame);
class SimulatorVisualizer {
 public:
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  std::shared_ptr<DRAWER::Drawer> _drawer;
  std::shared_ptr<ArticulatedBody> _body;
  Simulator& _sim;
  Eigen::Matrix<float,3,1> _cBody;
  Eigen::Matrix<float,3,1> _cBB;
  Eigen::Matrix<float,3,1> _cBBC;
  SimulatorVisualizer(Simulator& sim);
  virtual ~SimulatorVisualizer();
  int _width,_height;
  float _lightsize;
  Eigen::Matrix<float,3,1>_LightDiffuse;
  std::vector<Eigen::Matrix<float,3,1>>_ArticulateDiffuse;
  void setLightDiffuse(Eigen::Matrix<float,3,1>LightDiffuse);
  void setArticulateDiffuse(std::vector<Eigen::Matrix<float,3,1>>ArticulateDiffuse);
  void setLightSize(float lightsize);
  void visualize(const std::vector<Vec>* theta,int* frame);
};
}
