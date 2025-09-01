#include <iostream>
#include <cmath>
using namespace std;

double N(double x) {
    if (x < 0) return (1-N(-x));
    double p = 0.2316419;
    double b1 = 0.319381530;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double t = 1/(1+p*x);
    return 1-1/sqrt(2* M_PI)*exp(-x*x/2)*(b1*t+b2*pow(t,2)+b3*pow(t,3)+b4*pow(t,4)+b5*pow(t,5));

}

double callEuroBS(double t, double x, double K, double T, double r, double sigma) {
    double theta = T-t;
    double d1 = (log(x/K)+(r+sigma*sigma/2)*theta)/(sigma*sqrt(theta));
    double d2 = d1-sigma*sqrt(theta);
    return x*N(d1)-K*exp(-r*theta)*N(d2);
}

double putEuroBS(double t, double x, double K, double T, double r, double sigma) {
    double theta = T-t;
    double d1 = (log(x/K)+(r+sigma*sigma/2)*theta)/(sigma*sqrt(theta));
    double d2 = d1-sigma*sqrt(theta);
    return K*exp(-r*theta)*N(-d2)-x*N(-d1);
}
double deltaCallEuroBS(double t, double x, double K, double T, double r, double sigma) {
    double theta = T-t;
    double d1 = (log(x/K)+(r+sigma*sigma/2)*theta)/(sigma*sqrt(theta));
    return N(d1);
}

double deltaPutEuroBS(double t, double x, double K, double T, double r, double sigma) {
    double theta = T-t;
    double d1 = (log(x/K)+(r+sigma*sigma/2)*theta)/(sigma*sqrt(theta));
    return N(d1)-1;
}

double f(double x, double alpha, double K, double T, double r, double sigma) {
    return abs(alpha)*(K-putEuroBS(0,x,K,T,r,sigma))/(deltaPutEuroBS(0,x,K,T,r,sigma)+1+abs(alpha));
}

double xstar(double epsilon, double alpha, double K, double T, double r, double sigma) {
    int n = ceil(log2(K/epsilon));
    double a = 0;
    double b = K;
    for (int i = 1; i <= n; i++) {
        double mid = (a+b)/2;
        double fx = f(mid, alpha, K, T, r, sigma)-mid;
        if (fx>0) {
            a = mid;
        }else {
            b = mid;
        }
    }
    return (a+b)/2;
}

double putAmerBSMacMilan(double t, double x, double K, double T, double r, double sigma, double epsilon) {
    double sigma2 = sigma*sigma;
    double alpha = 1/2-r/sigma2-sqrt(pow(r/sigma2+1/2,2)+2/(T*sigma2));
    double x_star = xstar(epsilon,alpha,K,T,r,sigma);
    double vx;
    if (x >= x_star) {
        double lambda = (max(K-x_star,0.0)-putEuroBS(0,x_star,K,T,r,sigma))/pow(x_star,alpha);
        vx = lambda*pow(x,alpha);
    }else {
        vx = max(K-x,0.0)-putEuroBS(0,x,K,T,r,sigma);
    }
    return vx+putEuroBS(0,x,K,T,r,sigma);
}

int main() {
    double x = 100;
    double K = 110;
    double T = 2;
    double r = 0.03;
    double sigma = 0.2;
    cout << putEuroBS(0,x,K,T,r,sigma) << endl;
    cout << putAmerBSMacMilan(0,x,K,T,r,sigma,1e-3) << endl;
    return 0;
}