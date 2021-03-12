#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h> 
#include <string>

using namespace std;

// output for circle function
struct str_circle{ 
    double R;
    double xm;
    double ym;
};
struct koeff{
    vector<double> ak;
    vector<double> bk;
};

// function declarations
str_circle circle(double,double,double,double,double,double);
void circle(bool);

double integ(vector<double>,vector<double>,string);
void integ();
vector<double> func(vector<double>);

koeff fourier(vector<double>,vector<double>);
void fourier();
void makefile();

void solver(bool);
vector<double> gauss(vector<double>,vector<double>);

void ls(bool);
vector<double> ls(vector<double>, vector<double>,int);

void spline();
vector<double> spline(vector<double>,vector<double>, vector<double>);

// main function
int main(){
    cout << "This is the program containing the lecture tasks!" << endl;
    string program = "spline";
    if(program == "circle"){
        circle(false);
    } else if (program=="integration"){
        integ();
    } else if (program=="fourier"){
        fourier();
    } else if (program=="solver"){
        // gauss only
        solver(true);
    } else if(program=="least-square"){
        ls(true);
    } else if(program=="spline"){
        // not working -> there is a mistake in the calculation of coefficients giving nan values on first iteration
        spline();
    }
    //system("pause"); // activate for pause at the program end when using external consol!
    return 0;
}


// midpoint and radius calculation of three points
str_circle circle(double x1,double x2,double x3,double y1 ,double y2,double y3){
    // variable declaration
    str_circle result;
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
    result.xm = xm; result.ym = ym; result.R = (pow(x1-xm,2)+pow(y1-ym,2));
    return result;
}

void circle(bool input){
    // variable declaration
    str_circle c;
    double x1,x2,x3,y1,y2,y3;
    // input from console
    if(input){
        cin >> x1 >> x2 >> x3 >> y1 >> y2 >> y3; 
    } else{x1= 0;x2= 5;x3= 10;y1= 5;y2= 0;y3= 5;}
    // point switch for mathematical definition 
    if(x1-x2==0){
        double tmp;
        tmp = x2; x2 = x3; x3 = tmp; tmp = y2; y2 = y3; y3 = tmp;
    }
    if(y1-y2==0){
        double tmp;
        tmp = x1; x1 = x3; x3 = tmp; tmp = y1; y1 = y3; y3 = tmp;
    } 
    if(y1-y3==0){
        double tmp;
        tmp = x1; x1 = x2; x2 = tmp; tmp = y1; y1 = y2; y2 = tmp;
    } 
    // check valid input and compute function
    if((x1-x2==0)||(y1-y2==0)||(y1-y3==0)){
        cout << "Invalid input" << endl;
    }else{
        c = circle(x1,x2,x3,y1,y2,y3);
        // console output construct
        cout << "Radius: " << c.R << "\t" << "(Xm,Ym): " << c.xm << "," << c.ym << endl;
    }
}

// integration of custom function via different integration rules
vector<double> func(vector<double> x){
    vector<double> y(x.size());
    for(int i=0;i<x.size();i++){
        // customizable function
        y[i] = pow(x[i], 3) - 2* pow(x[i], 2) -x[i] + 2; // x^3 - 2x^2- x + 2
    }
    return y;
}

double integ(vector<double> x,vector<double> y,string method){
    double A = 0;
    if(method=="trapez"){
        for(int i=0;i<x.size()-1;i++){
            A += 0.5 * (x[i+1]-x[i]) * (y[i+1] + y[i]); 
        }
    }else if(method=="left-bound"){
        for(int i=0;i<x.size()-1;i++){
            A += y[i] * (x[i+1]-x[i]);
        } 
    }else if(method=="simpson"){
        int size=x.size();
        // trapez integration for outsiders
        while(size%3!=0){
            A += 0.5 * (x[size-1]-x[size-2]) * (y[size-1] + y[size-2]); 
            size--;
        }
        for(int i=1;i<x.size();i+=2){
            A+=(x[i+1]-x[i-1])/6*(y[i-1]+4*y[i]+y[i+1]);
        }
    }else{
        cout << "Not a valid integration scheme!" << endl;
    }
    return A;
}

void integ(){
    // create x and y values based on dx and integral boundaries
    double dx = 0.5, l = -3, u = 7; // parameters
    vector<double> x, y;
    for(double i=l;i<=u;i+=dx){
        x.push_back(i);
    }
    y = func(x);
    // selection of rule for integration
    string method1 = "left-bound", method2 = "trapez", method3="simpson";
    // return values
    double a1, a2, a3;
    a1 = integ(x,y,method1);
    a2 = integ(x,y,method2);
    a3 = integ(x,y,method3);
    cout << "Left-Bound: " << a1 << endl << "Trapez: " << a2 << endl << "Simpson: "<< a3 << endl;
}

// fourier analysis of data file
koeff fourier(vector<double> x, vector<double> y){
    // define constants for fourier 
    double T = 5, om = 2*3.14159 / T;
    int k = 100; // number of koefficients
    // define variables for calculation
    koeff c;
    vector<double> ak, bk, xc ,yc ,ycos, ysin;
    double v1,v2 = 0,i1,i2;
    // cut data matrixes
    int size = x.size();
    for(int i=0;i<=size;i++){
        if(i>=size){
            cout << "Oscillation period exceeds provided data! Provide a wider range of data." << endl;
            return c;
        }else{
            if(x[i]>=T){
                for(int j=0;j<=i;j++){
                    xc.push_back(x[j]);
                    yc.push_back(y[j]);
                }
                break;
            }
        }
    }
    // calculate fourier coefficients
    for(int i=0;i<k;i++){
        ycos.clear(); ysin.clear();
        for(int j=0;j<yc.size();j++){
            ycos.push_back(y[j]*cos(i*om*xc[j]));
            ysin.push_back(y[j]*sin(i*om*xc[j]));
        }
        v1 = 2/T * integ(xc,ycos,"simpson");
        ak.push_back(v1);
        if(i!=0){
            v2 = 2/T * integ(xc,ysin,"simpson");
        }
        bk.push_back(v2);
    }
    // return koefficients
    c.ak = ak; c.bk = bk;
    return c;
}

void fourier(){
    // create input file
    // declare variables for storing data
    vector<double> x, y;
    double v1,v2;  
    // reading in from file  
    ifstream file ("fourier.dat");
    if(file.is_open()){
        while(!file.eof()){
            file >> v1 >> v2;
            x.push_back(v1);
            y.push_back(v2);
        }
        file.close();
    }else{
        cout << "File could not be opened!" << endl;
    }
    // fourier analysis
    koeff k = fourier(x,y);
    if (k.ak.size()>5){
        ofstream file ("koefficients.dat");
        if(file.is_open()){
            for(int i=0;i<k.ak.size();i++){
                file << k.ak[i] << "\t" << k.bk[i] << endl;
            }
            file.close();
        }
        cout << "Koefficients saved in \"koefficients.dat\"!" << endl;
    }else{
        cout << "The coefficients are" << endl << "ak" << "\t" << "bk" << endl;
        for(int i=0;i<k.ak.size();i++){
            cout << k.ak[i] << "\t" << k.bk[i] << endl;
        }
    }
}

// gauss solver for system of equations
void solver(bool given){
    // declare systems of equations
    vector<double> A, b, x;
    if(given){
        A = {1,-1,-2,-2,1,-6,1,0,-2};
        b = {0,0,3};
    }
    if(floor(A.size()/b.size())!=b.size()){
        cout << "System of equations is not constructed properly!" << endl;
    }else{
        x = gauss(A,b);
    }
    cout << "Solution:" << endl;
    for(int i=0;i<x.size();i++){
        cout << "X" << i << ": \t" << x[i] << endl;
    }
}

vector<double> gauss(vector<double> A,vector<double> b){
    vector<double> x(b.size());
    double f = 0;
    int size = b.size();
    // construct upper triangular matrix by variable elimination
    for(int i=0;i<size-1;i++){ // variabel elminiation iteration
        for(int j=i+1;j<size;j++){ // row iteration
            f = A[j*size+i]/A[i*size+i]; // factor of elimination
            A[j*size+i] = 0; // set element to eliminate 0
            for(int k=i+1;k<size;k++){
                A[j*size+k] = A[j*size+k] - A[i*size+k] * f;
            }
            b[j] = b[j] - b[i] * f;
        }
    }
    // solve upper triangular matrix
    x[size-1] = (b[size-1]/A[(size-1)*size-1]); // first variable
    for(int i=size-2;i>=0;i--){ // iteration over variables backwards
        x[i] = b[i];
        for(int j= i+1; j<size;j++){
            x[i] = x[i] - A[ i * size + j] * x[j];
        }
        x[i] = x[i]/A[i * size + i];
    }
    return x;
}

// polynomial least square method
void makefile(){
    // create x and y values based on dx and integral boundaries
    double dx = 0.1, l = -1, u = 5; // parameters, dont change lower bound!!!
    vector<double> x, y;
    for(double i=l;i<=u;i+=dx){
        x.push_back(i);
    }
    y = func(x);
    // write data matrixes to a file
    ofstream file ("function.dat");
    if(file.is_open()){
        for(int i=0;i<x.size();i++){
            file << x[i] << "\t" << y[i] << endl;
        }
        file.close();
    }
}

void ls(bool make){
    int n = 5; //order of polynomial
    // creating datafile from function
    if(make){
        makefile();
    }
    //reading input data
    ifstream file ("xy_data.dat");
    vector<double> x, y;
    double v1, v2;
    if(file.is_open()){
        while(!file.eof()){
            file >> v1 >> v2;
            x.push_back(v1);
            y.push_back(v2);
        }
        file.close();
    }else{
        cout << "File could not be read!" << endl;
    }
    // ls scheme
    vector<double> c = ls(x,y,n);
    // output of coefficients   
    cout << "The coefficients for the polynomial of order " << n << " are:" << endl;
    for(int i=0;i<c.size();i++){
        cout << "c" << i << ": \t" << c[i] << endl;
    }
    cout << endl;
}

vector<double> ls(vector<double> x, vector<double> y, int n){
    vector<double> A(pow(n+1,2)), b(n+1);
    vector<double> c(n+1);
    for(int i=0;i<n+1;i++){ // iteration over coefficients
        b[i] = 0; // set initial rhs to zero
        for(int j=0;j<n+1;j++){ // iteration over matrix columns
            A[i*(n+1)+j] = 0;
            for(int h=0;h<x.size();h++){
                A[i*(n+1)+j] += pow(x[h],i) * pow(x[h],j);
            }
        }
        for(int h=0;h<x.size();h++){
            b[i] += y[h] * pow(x[h],i);
        }
    }
    c = gauss(A,b);
    return c;
}

// spline approximation for function points
vector<double> tom_solve(vector<double> a, vector<double> b, vector<double> c, vector<double> d){
    vector<double> x(a.size());
    int iu = a.size(), il = 0;
    for(int i=il+1; i<iu; i++){
        d[i] = d[i] - b[i]/d[i-1] * a[i-1];
        c[i] = c[i] - b[i]/d[i-1] * c[i-1];
    }
    x[iu] = c[iu]/d[iu];
    for(int i=1;i<=(iu-il);i++)
    {
        x[iu - i] = (c[iu - i] - a[iu - i] * x[iu - i+1])/d[iu - i];
    }
    return x;
}

void spline(){
    // reading data from file
    vector<double> x,y;
    double v1,v2;
    ifstream file ("xy_data.dat");
    if(file.is_open()){
        while(!file.eof()){
            file >> v1 >> v2;
            x.push_back(v1);
            y.push_back(v2);
        }
        file.close();
    }else{
        cout << "File could not be properly opened!" << endl;
    }
    // make discrete spline approximation with n points
    int n = 1000;
    vector<double> xt(n), yt(n);
    double dx = (x[x.size()-1] - x[0])/n;
    for(int i=0;i<xt.size();i++){
        xt[i] = x[0] + i * dx;
    }
    yt = spline(x,y,xt);
	// output to file
    ofstream out ("splinepoints.dat");
    if(out.is_open()){
        for(int i=0;i < xt.size();i++){
            out << xt[i] << "\t" << yt[i] << endl;
        }
        out.close();
        cout << "The " << n << "-splinepoints saved in splinepoints.dat!" << endl;
    }
}

vector<double> spline(vector<double> x,vector<double> y, vector<double> xs){
    vector<double> ys;
    int n = x.size();
    // declare sub diagonals
    vector<double> a(n),b(n),c(n),d(n),rhs(n);
    // set diagonal boundaries
    a[0] = 0; b[0] = 0; d[0] = 1; rhs[0] = 0;
    a[n-1] = 0; b[n-1] = 0; d[n-1] = 1; rhs[n-1] = 0;
    for(int i=1;i<n-1;i++){ // number of splines
        d[i] = 2* (2*x[i]-x[i-1]-x[i+1]);
        b[i] = x[i] - x[i-1];
        a[i] = x[i+1] - x[i]; 
        rhs[i] = 6/(x[i+1]-x[i])*(y[i+1]-y[i]) - 6/(x[i]-x[i-1])*(y[i]-y[i-1]);
    }
    // solve system of equations
    vector<double> l = tom_solve(a,b,d,rhs);
    // calculate spline coefficients
    for(int i=0; i < n-1; i++){
	  	a[i] = 1/(6 * (x[i+1]-x[i])) * (l[i+1] - l[i]);
	  	b[i] = l[i]/2;
	  	c[i] = (y[i+1] - y[i])/(x[i+1]-x[i]) - (x[i+1]-x[i])/6 * (l[i+1] + 2*l[i]);
	  	d[i] = y[i];
	}
	// calculate spline approximated y
    for(int k=0;k<xs.size();k++){
        for(int i=0; i < n-1; i++){
	  	    if((xs[k] >= x[i]) && (xs[k]<=x[i+1])){
                cout << a[i] << " " << b[i] << " " << c[i] << " " << d[i];
	  		    ys[k] = a[i] * pow( xs[k] - x[i], 3) + b[i] * pow( xs[k] - x[i], 2 ) + c[i] * ( xs[k] - x[i] ) + d[i];
	  	    }
	    }
    }
    return ys;
}

