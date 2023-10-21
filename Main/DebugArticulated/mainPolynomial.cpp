#include <Utils/SOSPolynomial.cpp>

using namespace PHYSICSMOTION;

typedef PolyXAM PolyXA;
void debugPolynomial() {
  PolyXA a("2*a1^2 +3*a2^2 +2*a1*a2 +2*a3^2 +3*a4^2 +2*a3*a4",true);
  std::cout << "a=" << a.toString() << std::endl;
  PolyXA x("2*x1^2 +3*x2^2 +2*x1*x2 +2*x3^2 +3*x4^2 +2*x3*x4",true);
  std::cout << "x=" << x.toString() << std::endl;
  std::cout << "a*x=" << (a*x).toString() << std::endl;
  //basic op:
  PolyXA ax=a*x;
  std::cout << "ax+x5*a5=" << (ax+PolyXA("x5*a5",true)).toString() << std::endl;
  std::cout << "ax-x5*a5=" << (ax-PolyXA("x5*a5",true)).toString() << std::endl;
  std::cout << "ax*x5*a5=" << (ax*PolyXA("x5*a5",true)).toString() << std::endl;
  //remove zero
  ax+=PolyXA("0.000001*x7*a7",true);
  std::cout << "ax+0.000001*x7*a7=" << ax.toString() << std::endl;
  std::cout << "removeZero=" << ax.removeZero(1e-3f).toString() << std::endl;
  std::cout << "removeVariableId(0)=" << ax.removeVariableId(0).toString() << std::endl;
  std::cout << "varRemap=" << ax.varRemap({5,6,7,8,9,10,11,12}).toString() << std::endl;
  std::cout << "rename=" << ax.rename({12,11,10,9,8,7,6,5}).toString() << std::endl;
  std::cout << "affineTransXId=" << ax.affineTransXId<'x'>(50,10).toString() << std::endl;
  //transform
  std::cout << "linearConstraint=" << ax.linearConstraint(1,PolyXA("5*x2+6*x3",true)).toString() << std::endl;
  {
    std::unordered_map<int,typename ScalarOfT<PolyXA>::Type> lib;
    lib[1]=0.2;
    std::cout << "linearTransform=" << ax.linearTransform(lib).toString() << std::endl;
  }
  {
    std::unordered_map<int,PolyXA> lib;
    lib[0]=PolyXA("a0",true);
    lib[1]=PolyXA("5*x2+6*x3",true);
    lib[2]=PolyXA("5*x2+6*x3",true);
    lib[3]=PolyXA("a3",true);
    lib[4]=PolyXA("a4",true);
    lib[5]=PolyXA("a5",true);
    lib[6]=PolyXA("a6",true);
    lib[7]=PolyXA("a7",true);
    std::cout << "linearTransform=" << ax.linearTransform(lib).toString() << std::endl;
  }
  //gradient/hessian
  PolyXA::VECP g=ax.gradientV();
  for(int i=0; i<g.size(); i++)
    std::cout << "gradient[" << i << "]=" << g[i].toString() << std::endl;
  PolyXA::MATP h=ax.hessianM();
  for(int i=0; i<h.rows(); i++)
    for(int j=0; j<h.cols(); j++)
      std::cout << "hessian[" << i << "," << j << "]=" << h(i,j).toString() << std::endl;
}
void debugPolynomialSolve() {
  PolyXA a("a",true);
  PolyXA x0("x0",true);
  PolyXA x1("x1",true);
  PolyXA x2("x2",true);
  PolyXA x3("x3",true);
  PolyXA x01=x0*(PolyXA(1)-a)+x1*a;
  PolyXA x12=x1*(PolyXA(1)-a)+x2*a;
  PolyXA x23=x2*(PolyXA(1)-a)+x3*a;
  PolyXA x012=x01*(PolyXA(1)-a)+x12*a;
  PolyXA x123=x12*(PolyXA(1)-a)+x23*a;
  PolyXA bezier=x012*(PolyXA(1)-a)+x123*a;
  PolyXA bezierRef("1.5*a^3+1.2*a^2+1.3*a+1.7",true);
  std::cout << "BezierCurve=" << bezier.toString() << std::endl;
  std::cout << "TargetCurve=" << bezierRef.toString() << std::endl;
  PolyXA::COLD sol=PolyXA::solve({bezier}, {bezierRef});
  std::cout << "Coefficient=" << sol.transpose() << std::endl;
  std::cout << "Plugged-in=" << bezier.eval<'x'>(sol).toString() << std::endl;
}
int main(int argc,char** argv) {
  debugPolynomial();
  debugPolynomialSolve();
  return 0;
}
