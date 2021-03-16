// document for AUP exercise tasks as preparation for the exam WS20/21
#include<fstream>
#include<iostream>
#include<math.h>
#include<vector>
#include<string>

using namespace std;

// function declaration for each task
void t1(); void t2(); void t3(); void t4(); void t5(); void t6(); void t7(); void t8();

const double PI = 3.1415;

// subroutine declarations for task 1
struct Circle{
    double r;
    double xm;
    double ym;
};
Circle circle(double,double,double,double,double,double);

// subroutine declarations for task 2

// subroutine declarations for task 3
double integration(vector<double>,vector<double>);
// subroutine declarations for task 4
vector<double> gauss(vector<vector<double>>, vector<double>);
// subroutine declarations for task 5
double cw(double);
// subroutine declarations for task 6
vector<double> grenzschicht(vector<double>, double);
double H(double);
// subroutine declarations for task 7

// subroutine declarations for task 8


int main(){
    int task = 8;  
    switch(task){
        case 1: t1(); break;
        case 2: t2(); break;
        case 3: t3(); break;
        case 4: t4(); break;
        case 5: t5(); break;
        case 6: t6(); break;
        case 7: t7(); break;
        case 8: t8(); break;
        default: cout << "Not a valid task number!" << endl; break;
    }
    return 0;
}

// implementation of task functions
void t1(){
    // 3 points build one circle. Each subsequent circle takes the last point of the last circle as first point 
    // to form a curve. S.t. the formula for successful number of points to circles is  2 * k + 1 with k > 0 
    // and k being the number of circles/ circle center points. A circle is defined by its midpoint coordinates
    // and its radius. For the given set of points in kreis.dat the midpoints and radius should be determined.
    vector<double> x,y;
    double v1,v2;
    ifstream file ("data/kreis.dat");
    if(file.is_open()){
        while(!file.eof()){
            file >> v1 >> v2;
            x.push_back(v1);
            y.push_back(v2);
        }
        file.close();
    }
    int n = x.size();
    Circle ct;
    ofstream out ("data/RMP.dat");
    if(out.is_open()){
        for(int i=0;i<x.size()-2;i=i+2){
            ct = circle(x[i],y[i],x[i+1],y[i+1],x[i+2],y[i+2]);
            out << ct.r << "\t" << ct.xm << "\t" << ct.ym << "\t" << endl;
        }
        out.close();
    }
    cout << "Execution of Task 1 successful!" << endl;
    cout << "Results of Task 1 in RMP.dat!" << endl;
}

Circle circle(double x1,double y1,double x2,double y2,double x3,double y3){
    // variable declaration
    struct Circle result;
    double xm,ym,r,x12,x13,xq12,xq13,y12,y13,yq12,yq13,a,b;
    // calculation of variable combinations for formulas
    x12 = x1 - x2; x13 = x1 - x3; y12 = y1 - y2; y13 = y1 - y3;
    xq12 = pow(x1,2) - pow(x2,2); xq13 = pow(x1,2) - pow(x3,2); yq12 = pow(y1,2) - pow(y2,2); yq13 = pow(y1,2) - pow(y3,2);
    // calculation of end terms for xm
    a = x12/y12 - x13/y13;
    b = (xq12+yq12)/(2*y12) - (xq13+yq13)/(2*y13);
    // calculation of midpoint using formulas
    xm = b/a; ym = (-xq12 + 2*x12*xm-yq12)/(-2*y12);
    // setting structure properties
    result.xm = xm; result.ym = ym; result.r = sqrt(pow(x1-xm,2)+pow(y1-ym,2));
    return result;
}

void t2(){
    double a = 1;
    double b = 1.5;
    double x = (a+b)/2;
    double fx = sin(x*x/(2*PI));
    double fd;
    int n=0;
    while(abs(fx)>0.000001){
        fx = sin(x*x/2*PI);
        fd = x * PI * cos(x*x/2*PI);
        x = x - fx / fd;
        if(n>10000){
            cout << "No satisfying null point was found before max-iter!" << endl;
            break;
        }
    }
    cout << "Execution of Task 2 successful!" << endl;
    cout << "The null point of sin(x²/2*PI) is: x=" << x << "!" << endl; 
}

void t3(){
    // fourier formula is a linear combination of sin and cos functions with variable weight
    // and set multiple of a base frequency. The coefficients can be determined in a discrete 
    // system: b_k = T/2 \integ_0^T{f(t)*sin(k*w*t)dt} and a_k = T/2 \integ_0^T{f(t)*cos(k*w*t)dt}
    // for the coefficient k with k >= 0. Therefore we need an integration scheme for discrete data
    // points. We can choose trapez integration or integration by simpson rule with boundary 
    // cases being implemented with trapez rule.  
    double T = 5;
    double o = (2 * 3.1415)/T;
    int k = 100; // number of coefficients
    vector<double> ak(k), bk(k);
    ifstream file ("data/signal.dat");
    vector<double> y,t;
    double v1,v2;
    if(file.is_open()){
        while(!file.eof()){
            file >> v1 >> v2;
            t.push_back(v1);
            y.push_back(v2);
        }
        file.close();
    }
    // coefficient calculation
    bk[0] = 0;
    vector<double> ycos(y.size(),0), ysin(y.size(),0);
    for(int i=0;i<k;i++){
        // value iteration for integral values
        for(int h=0;h<y.size();++h){
            ycos[h] = y[h] * cos(i*o*t[h]);
            ysin[h] = y[h] * sin(i*o*t[h]);
        }
        // loop for bks using integration
        ak[i] = (2 / T) * integration(t,ycos);
        // loop for aks
        if(i!=0){
            bk[i] = (2 / T) * integration(t,ysin);
        }
    }
    // change coefficient for frequency > 3000 Hz corresponding to k > 3
    /*for(int i=3;i<k;i++){
        ak[i] /= 2;
        bk[i] /= 2;
    }*/
    // computation of dampened signal values
    double value;
    vector<double> ya(y.size());
    for(int i=0;i<y.size();i++){
        value = 0;
        for(int j=1;j<k;j++){
            value += bk[j] * sin(j*o*t[i]);
            value += ak[j] * cos(j*o*t[i]);
        }
        ya[i] = (ak[0]/2) + value;
    }

    // file output
    ofstream out ("data/damped_signal.dat");
    if(out.is_open()){
        for(int i=0;i<y.size();i++){
            out << t[i] << "\t" << y[i] << "\t" << ya[i] << endl;
        }
        out.close();
    }
    cout << "Execution of Task 3 successful!" << endl;
}

double integration(vector<double> x, vector<double> y){
    // simpson integration with trapez boundary
    double A = 0;
    int size=x.size();
    // trapez integration for outsiders
    while(size%3!=0){
        A += 0.5 * (x[size-1]-x[size-2]) * (y[size-1] + y[size-2]); 
        size--;
    }
    for(int i=1;i<x.size();i+=2){
        A+=(x[i+1]-x[i-1])/6*(y[i-1]+4*y[i]+y[i+1]);
    }
    return A;
}

void t4(){
    int n = 11; // polynomial degree
    // reading data from file
    vector<double> x,y;
    double v1,v2;
    ifstream file ("data/wing.dat");
    if(file.is_open()){
        while(!file.eof()){
            file >> v1 >> v2;
            x.push_back(v1);
            y.push_back(v2);
        }
        file.close();
    }
    int m = x.size();
    // construction of system of equations
    vector<double> a(n,0), rhs(n,0);
    vector<vector<double>> X(n,vector<double>(n,0));
    for(int i=0; i<m;i++){
        X[0][0] += sqrt(x[i]*sqrt(x[i]));
        rhs[0] += sqrt(x[i])*y[i];
    }
    for (int row = 1; row < n; ++row) {
        // construct first row,col and rhs
        for (int i = 0; i < m; ++i) {
            X[row][0] += pow(x[i],row) * sqrt(x[i]);
            X[0][row] += pow(x[i],row) * sqrt(x[i]); // symmetry
            rhs[row] += pow(x[i],row) * y[i];
        }
        for (int col = 1; col < n; ++col) {
            for (int i = 0; i < m; ++i) {
                X[row][col] += pow(x[i],row)*pow(x[i],col);
            }
        }
        
    }
    a = gauss(X,rhs);
    // output file
    int d = 100;
    double dx = (x.back()-x.front()) / (d-1);
    double xo,yo;

    ofstream out ("data/poly_wing.dat");
    if(out.is_open()){
        for(int k=0;k<d;k++){
            xo = k*dx+x.front();
            yo = a[0] * sqrt(xo);
            for(int i=1;i<n;i++){
                yo += a[i] * pow(xo,i);
            }
            out << xo << "\t" << yo << endl;
        }
        out.close();
    }
    cout << "Polynomial approximation of wing.dat saved in poly_wing.dat!" << endl;
    cout << "Execution of Task 4 successful!" << endl;
}

vector<double> gauss(vector<vector<double>> M, vector<double> R) {
    int n = R.size();
    vector<double> X(n, 0.0);
    int i, j, k;
    //forward substitution
    for (i = 0; i < n - 1; i++) {
        for (k = i + 1; k < n; k++) {
            for (j = i; j < n; j++) {
                X[j] = -M[i][j] * M[k][i] / M[i][i];
            }
            R[k] = R[k] - R[i] * M[k][i] / M[i][i];
            for (j = i; j <= n - 1; j++) {
                M[k][j] = M[k][j] + X[j];
            }
        }
    }
  // backwards substitution
    X[n - 1] = R[n - 1] / M[(n - 1)][n - 1];
    for (i = n - 2; i >= 0; i--) {
        X[i] = R[i];
        for (j = i + 1; j < n; j++) {
            X[i] = X[i] - M[i][j] * X[j];
        }
        X[i] = X[i] / M[i][i];
        }
    return X;
}

void t5(){
    double cw1 = 0.2, cw2 = 5;
    double C = 0.1, g = 9.81; // C ~ 1.22 * 1 / (2*70) 
    vector<double> X,Y,T;
    // initial conditions
    X.push_back(5000);
    Y.push_back(0);
    T.push_back(0);
    // k values
    double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y;
    double kx,ky;
    double dt = 0.1;
    double xi, yi;
    // iteration loop
    while(X.back() >= 0){
        T.push_back(T.back()+dt);    
        xi = X.back();
        yi = Y.back();
        // k
        k1y = g - C * cw(xi) * pow(yi,2);
        k1x = yi;
        k2y = g - C * cw(xi+dt/2*k1x) * pow(yi+dt/2*k1y,2);
        k2x = yi + dt/2*k1y;
        k3y = g - C * cw(xi+dt/2*k2x) * pow(yi+dt/2*k2y,2);
        k3x = yi + dt/2*k2y;
        k4y = g - C * cw(xi+dt*k3x) * pow(yi+dt*k3y,2);
        k4x = yi + dt*k3y;
        ky = (k1y+2*k2y+2*k3y+k4y)/6;
        kx = (k1x+2*k2x+2*k3x+k4x)/6;
        // append vectors
        X.push_back(X.back()-kx*dt);
        Y.push_back(Y.back()+ky*dt);
    }
    // output
    ofstream out ("data/jump.dat");
    if(out.is_open()){
        for(int i=0;i<X.size();i++){
            out << T[i] << "\t" << X[i] << "\t" << Y[i] << endl;
        }
        out.close();
    }
    cout << "Results saved in jump.dat" << endl;
    cout << "Execution of Task 5 successful!" << endl;
}

double cw(double x){
    if(x>2000){
        return 0.2;
    }else{
        return 5;
    }
}

void t6(){
    double hv = 3;
    int maxIter = 1000, count = 0;
    // intial guesses 
    double u1 = 1, u2 = 10;
    double uk;
    double h1 = H(u1) - hv; double h2 = H(u2) - hv;
    // calculation loop
    while(abs(h2) > 0.001 || count < maxIter){
        count++;
        uk = u2 - h2/(h2-h1)*(u2-u1); // newton step with approximated gradient
        cout << (h2-h1) << " " << (u2-u1) << " " << h2 << endl;
        // update parameters
        u1 = u2;
        h1 = h2;
        u2 = uk;
        h2 = H(uk) - hv;
        cout << (h2-h1) << " " << (u2-u1) << " " << h2 << endl;
        break;
    }   
    cout << "Optimales U: " << uk << "!" << endl;
    cout << "Execution of Task 6 successful!" << endl;
}
// integral function already implemented as integration

vector<double> grenzschicht(vector<double> Y, double U) {
    int n = Y.size();
    vector<double> u(n, 1); // real function not provided
    for (int i = 0; i < n; ++i) {
        u[i] = U*(1-(1-pow(Y[i]/0.1,2))); 
    }
    return u;
}

double H(double U){
    int n = 101; // discretisation for integral
    vector<double> u(n),y(n);
    vector<double> intDelta(n), intTheta(n);
    double theta, delta;
    double diff = 0.1;
    double dy = (diff+0) / (n-1);
    // initialize y
    for(int i=0;i<n;i++){
        y[i]= (i*dy);
    }
    // set u(y)
    u = grenzschicht(y,U);
    // calculate theta, delta
    for(int i=0;i<n;i++){
        intDelta[i] = 1-u[i]/U;
        intTheta[i] = (1-(u[i]/U))*(u[i]/U);
    }
    delta = integration(y, intDelta);
    theta = integration(y, intTheta);

    return delta/theta;
}

void t7(){
    double rho, lambda, cq, v, l, area, peri, alpha;
    double dx, dt, ca, cb;
    rho = 1000.0;
    lambda = 0.6;
    l = 20.0;
    cq = 4200.0;
    v = 0.1;
    area = 3.14 * pow(10, -4);
    peri = 0.0628;
    alpha = 2.5;

    ca = lambda / (rho * cq);
    cb = (alpha * peri) / (rho * area * cq);

    int diskret_x = 101;
    int diskret_t = 20000;

    dx = l / (diskret_x - 1);
    dt = 0.01;
  
    vector<double> T(diskret_x, 288.0);

    // boundary conditions
    T.front() = 318.0;
    vector<double> TB = T;
  
    ofstream out;
    for (int j = 0; j <= diskret_t; j++) {
        for (int i = 1; i < diskret_x - 1; ++i) {
            T[i] = ca * dt * ((TB[i-1] + 2. * TB[i] + TB[i + 1]) / (dx * dx)) - cb * (TB[i] - 288.0) * dt - v * (TB[i] - TB[i-1]) * dt / dx + TB[i];
            TB = T;
        }
        if (j % 2000 == 0) {
        out.open(string("data/output_") + std::to_string(j));
        for (int i = 0; i < diskret_x; ++i) {
            out << i * dx << "," << T[i] << std::endl;
        }
        out.close();
    }
  }
    cout << "Execution of Task 8 successful!" << endl;
}

void t8(){
    vector<double> y; int n = 101; // discretise
    double dy = 1./(n-1.);
    vector<double> u(n,0);
    for(int i=0;i<n;i++){
        y.push_back(i*dy);    
    }
    
    // explicit
    for(int i = 0;i<100000;i++){
        for(int j=1;j<n-1;j++){ // exclude boundaries!
            u[j] = (u[j+1] + u[j-1] + pow(dy,2))/2;
        }
    }
    // output
    ofstream out ("data/canal_explicit.dat");
    if(out.is_open()){
        for(int i=0;i<n;i++){
            out << y[i] << "," << u[i] << endl;
        }
        out.close();
    }

    // implicit
    vector<vector<double>> A(n,vector<double>(n,0));
    vector<double> rhs(n,0); double r = pow(dy,2);
    A[0][0] = -2; A[0][1] = 1;
    A[n-1][n-1] = -2; A[n-1][n-2] = 1;
    rhs[0] = 0; rhs[n-1] = 0;
    for(int i=1;i<n-1;i++){
        rhs[i] = -r;
        A[i][i] = -2;
        A[i][i+1] = 1;
        A[i][i-1] = 1;
    }
    u = gauss(A,rhs);
    out.open("data/canal_implicit.dat");
    if(out.is_open()){
        for(int i=0;i<n;i++){
            out << y[i] << "," << u[i] << endl;
        }
        out.close();
    }
    cout << "Execution of Task 7 successful!" << endl;  
}
