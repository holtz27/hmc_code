#ifndef SVMN_H_INCLUDED
#define SVMN_H_INCLUDED
mat svmnsim(vec beta, double mu, double phi, double sigma2, double y0);
vec mvrgaussian(int n);

class SVMNHMC{
public:

double lposttheta1(vec theta1, vec vH);
double Gtheta1(vec theta1, mat &G, mat &invG, mat &Gparmu, mat &Gparomega, mat &Gpargamma);
vec gradphitheta1(vec theta1, vec vH,mat invG,mat Gmu,mat Gomega,mat Ggamma);
vec nuHtheta1(mat Gmu, mat Gomega, mat Ggamma, mat invG, vec p);
vec updatetheta1(vec theta1prev, vec vH, int nfix, int nlpfrog, double epsilon);
vec gradHmomtheta1(mat invG, vec p);

double Gtheta2(vec theta2,vec vH,mat &G, mat &invG, mat &Gparbeta0,mat &Gpardelta,mat &Gparbeta2);
double lposttheta2(vec theta2, vec vH);
vec gradphitheta2(vec theta2, vec vH,mat invG,mat Gbeta0,mat Gdelta,mat Gbeta2);
vec gradHmomtheta2(mat invG,vec p);
vec nuHtheta2(mat Gbeta0,mat Gdelta,mat Gbeta2,mat invG, vec p);
vec updatetheta2(vec theta2prev,vec vH,int nfix,int nlpfrog,double epsilon);

double lpostvH(vec vH,vec theta1,vec theta2);
vec gradlpvH(vec vH,vec theta1,vec theta2);
vec updatevH(vec vHprev,vec theta1,vec theta2,int nlpgrog,double epsilon);

vec vY;
};

#endif // SVMN_H_INCLUDED
