#include "SimulatorVisualizer.h"
#include <Articulated/ArticulatedVisualizer.h>
#include <Environment/EnvironmentVisualizer.h>
#include <Simulator/Simulator.h>
#include <Simulator/PBDSimulator.h>

#include <TinyVisualizer/FirstPersonCameraManipulator.h>
#include <TinyVisualizer/CameraExportPlugin.h>
#include <TinyVisualizer/CaptureGIFPlugin.h>
#include <TinyVisualizer/ImGuiPlugin.h>
#include <TinyVisualizer/ShadowAndLight.h>
#include <TinyVisualizer/Camera3D.h>
#include <imgui/imgui.h>

namespace PHYSICSMOTION {
void visualizeSimulator(int argc,char** argv,Simulator& sim,int nSub,float maxStiffness) {
  typedef FLOAT T;
  DECL_MAT_VEC_MAP_TYPES_T
  using namespace DRAWER;

  //visualizer
  Drawer drawer(argc,argv);
  GLfloat forceCoef=1.f,IAlpha;
  Eigen::Matrix<GLfloat,6,1> ray;
  bool doSim=false,wire=false,lastWire=true,recording=false;
  std::shared_ptr<CameraExportPlugin> exporter(new CameraExportPlugin(GLFW_KEY_2,GLFW_KEY_3,"camera.dat"));
  std::shared_ptr<CaptureGIFPlugin> recorder(new CaptureGIFPlugin(GLFW_KEY_1,"record.gif",drawer.FPS()));
  drawer.addPlugin(exporter);
  drawer.addPlugin(recorder);

  std::shared_ptr<ArticulatedBody> body=sim.getBody();
  Eigen::Matrix<GLfloat,3,1> cBody(.3f,.3f,.3f);
  Eigen::Matrix<GLfloat,3,1> cBB(.3f,.2f,.1f);
  Eigen::Matrix<GLfloat,3,1> cBBC(.0f,.0f,.0f);
  std::shared_ptr<MeshShape> drag(new MeshShape);
  std::shared_ptr<CompositeShape> contact(new CompositeShape);
  std::shared_ptr<CompositeShape> bodyShapeC(new CompositeShape);
  std::shared_ptr<CompositeShape> envShapeC(new CompositeShape);
  auto bodyShape=visualizeArticulated(sim.getBody(),cBody,cBB,cBBC,false);
  auto envShape=visualizeEnvironment(sim.getShapes(),false);
  auto bodyWireShape=visualizeArticulated(sim.getBody(),cBody,cBB,cBBC,true);
  auto envWireShape=visualizeEnvironment(sim.getShapes(),true);
  bodyWireShape->setLineWidth(2);
  envWireShape->setLineWidth(2);

  drawer.addShape(bodyShapeC);
  drawer.addShape(envShapeC);
  drawer.addShape(contact);
  drawer.addCamera3D(90);
  drawer.getCamera3D()->setManipulator(std::shared_ptr<CameraManipulator>(new FirstPersonCameraManipulator(drawer.getCamera3D())));
  drawer.addLightSystem();
  drawer.getLight()->lightSz(10);
  auto bb=sim.getContact().getBB();
  GLfloat range=std::max<GLfloat>((GLfloat)bb.maxCorner().cwiseAbs().maxCoeff(),(GLfloat)bb.minCorner().cwiseAbs().maxCoeff());
  for(int x=-1; x<=1; x+=2)
    for(int y=-1; y<=1; y+=2)
      for(int z=-1; z<=1; z+=2)
        drawer.getLight()->addLight(Eigen::Matrix<GLfloat,3,1>(x,y,z)*range,
                                    Eigen::Matrix<GLfloat,3,1>(.5,.5,.5),
                                    Eigen::Matrix<GLfloat,3,1>(.5,.5,.5),
                                    Eigen::Matrix<GLfloat,3,1>(.5,.5,.5));
  drawer.setFrameFunc([&](std::shared_ptr<SceneNode>&) {
    if(wire!=lastWire) {
      lastWire=wire;
      if(wire) {
        bodyShapeC->addShape(bodyWireShape);
        envShapeC->addShape(envWireShape);
        bodyShapeC->removeChild(bodyShape);
        envShapeC->removeChild(envShape);
      } else {
        bodyShapeC->addShape(bodyShape);
        envShapeC->addShape(envShape);
        bodyShapeC->removeChild(bodyWireShape);
        envShapeC->removeChild(envWireShape);
      }
    }
    if(wire)
      updateArticulatedBody(bodyWireShape,sim.getBody(),ArticulatedBody::Vec(sim.pos().template cast<ArticulatedBody::T>()));
    else updateArticulatedBody(bodyShape,sim.getBody(),ArticulatedBody::Vec(sim.pos().template cast<ArticulatedBody::T>()));
    visualizeContact(contact,sim.getManifolds(),5,forceCoef);
    if(!doSim)
      return;
    for(int i=0; i<nSub; i++)
      sim.step();
    if(!sim.getDrags().empty()) {
      auto& d=sim.getDrags()[0];
      PBDArticulatedGradientInfo<T> pos(*body,sim.pos());
      //drag
      drag->clear();
      drag->addVertex(d._pt.template cast<GLfloat>());
      drag->addVertex((ROTI(pos._TM,d._jid)*d._ptL+CTRI(pos._TM,d._jid)).template cast<GLfloat>());
      drag->addIndexSingle(0);
      drag->addIndexSingle(1);
      drag->setMode(GL_LINES);
      drag->setColorDiffuse(GL_LINES,.5f,1.f,.5f);
      drag->setLineWidth(2);
    }
  });
  drawer.setKeyFunc([&](GLFWwindow*,int key,int,int action,int,bool captured) {
    if(captured)
      return;
    else if(key==GLFW_KEY_R && action==GLFW_PRESS)
      doSim=!doSim;
    else if(key==GLFW_KEY_E && action==GLFW_PRESS)
      wire=!wire;
  });
  drawer.setMouseFunc([&](GLFWwindow* wnd,int button,int action,int,bool captured) {
    if(captured || wire)
      return;
    else if(button==GLFW_MOUSE_BUTTON_2 && action==GLFW_PRESS) {
      sim.clearDrag();
      IAlpha=100;
      double x=0,y=0;
      glfwGetCursorPos(wnd,&x,&y);
      ray=drawer.getCameraRay(x,y);
      for(int i=0,k=0; i<body->nrJ(); i++)
        if(body->joint(i)._mesh) {
          auto joint=std::dynamic_pointer_cast<CompositeShape>(bodyShape)->getChild(k++);
          if(joint->rayIntersect(ray,IAlpha)) {
            PBDArticulatedGradientInfo<T> pos(*body,sim.pos());
            Simulator::DragEnergy d;
            d._jid=i;
            d._pt=(ray.template segment<3>(0)+ray.template segment<3>(3)*IAlpha).template cast<T>();
            d._ptL=ROTI(pos._TM,i).transpose()*(d._pt-CTRI(pos._TM,i));
            //drag
            drag->clear();
            drag->addVertex(d._pt.template cast<GLfloat>());
            drag->addVertex(d._pt.template cast<GLfloat>());
            drag->addIndexSingle(0);
            drag->addIndexSingle(1);
            drag->setMode(GL_LINES);
            drag->setColorDiffuse(GL_LINES,.5f,1.f,.5f);
            drag->setLineWidth(2);
            drawer.addShape(drag);
            sim.addDrag(d);
            break;
          }
        }
    } else if(button==GLFW_MOUSE_BUTTON_2 && action==GLFW_RELEASE) {
      drawer.removeShape(drag);
      sim.clearDrag();
    }
  });
  drawer.setMotionFunc([&](GLFWwindow* wnd,double,double,bool captured) {
    if(captured)
      return;
    else if(!sim.getDrags().empty()) {
      auto& d=sim.getDrags()[0];
      double x=0,y=0;
      glfwGetCursorPos(wnd,&x,&y);
      ray=drawer.getCameraRay(x,y);
      PBDArticulatedGradientInfo<T> pos(*body,sim.pos());
      d._pt=(ray.template segment<3>(0)+ray.template segment<3>(3)*IAlpha).template cast<T>();
      //drag
      drag->clear();
      drag->addVertex(d._pt.template cast<GLfloat>());
      drag->addVertex((ROTI(pos._TM,d._jid)*d._ptL+CTRI(pos._TM,d._jid)).template cast<GLfloat>());
      drag->addIndexSingle(0);
      drag->addIndexSingle(1);
      drag->setMode(GL_LINES);
      drag->setColorDiffuse(GL_LINES,.5f,1.f,.5f);
      drag->setLineWidth(2);
    }
  });
  drawer.addPlugin(std::shared_ptr<Plugin>(new ImGuiPlugin([&]() {
    drawer.getCamera3D()->getManipulator()->imGuiCallback();
    ImGui::Begin("Simulator");
    ImGui::Checkbox("wireframe",&wire);
    ImGui::Checkbox("simulating",&doSim);
    ImGui::SameLine();
    if(ImGui::Button("reset position")) {
      sim.setPos(Vec::Zero(body->nrDOF()));
      sim.setVel(Vec::Zero(body->nrDOF()));
      sim.setTime(0);
    }
    ImGui::SameLine();
    if(ImGui::Button("detect collision"))
      sim.detectCurrentContact();
    ImGui::SameLine();
    recording=recorder->recording();
    if(ImGui::Checkbox("record",&recording)) {
      if(recording)
        recorder->startRecording();
      else recorder->stopRecording();
    }
    if(ImGui::Button("save camera"))
      exporter->saveCamera();
    ImGui::SameLine();
    if(ImGui::Button("load camera"))
      exporter->loadCamera();
    ImGui::SliderFloat("force amplitude",&forceCoef,0.f,10.f);
    ImGui::Separator();
    //friction coefficient
    GLfloat fri=(GLfloat)sim.getJointPhysicsParameter(0)._friction;
    if(ImGui::SliderFloat("friction coefficient",&fri,0.f,1.f))
      for(int k=0; k<body->nrJ(); k++)
        sim.getJointPhysicsParameter(k)._friction=fri;
    //contact stiffness
    GLfloat kc=(GLfloat)sim.getJointPhysicsParameter(0)._kc;
    if(ImGui::SliderFloat("contact stiffness",&kc,0.f,maxStiffness))
      for(int k=0; k<body->nrJ(); k++)
        sim.getJointPhysicsParameter(k)._kc=kc;
    //gravity
    GLfloat g[3]= {0,0,-9.81f};
    if(ImGui::SliderFloat3("gravity",g,-10.f,10.f))
      sim.setGravity(Vec3T(g[0],g[1],g[2]));
    ImGui::SameLine();
    if(ImGui::Button("reset")) {
      g[0]=0.f;
      g[1]=0.f;
      g[2]=-9.81f;
      sim.setGravity(Vec3T(g[0],g[1],g[2]));
    }
    ImGui::End();
  })));
  drawer.mainLoop();
}
}
