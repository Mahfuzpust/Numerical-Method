//bisection
#include<bits/stdc++.h>
using namespace std;
#define e 0.01

double func(double x)
{
    return x*x*x - x*x + 2;
}

 double c;

void bisection(double a, double b)
{
    if (func(a) * func(b) >= 0)
    {
        cout << "You have not assumed right a and b\n";
        return;
    }

    double c = a;
    while ((b-a) >= e)
    {
        c = (a+b)/2;

        if (func(c) == 0.0){
            cout <<"Root : "<< c <<endl;
            break;
        }
        else if (func(c)*func(a) < 0){
            cout <<"Root : "<< c <<endl;
            b = c;
        }
        else{
            cout <<"Root : "<< c <<endl;
            b = c;
        }
    }
    cout << endl;
    cout << "The value of root is : " << c;
}

int main()
{
    //double a =-200, b = 300;
    double a,b;
    cin >> a >>b;
    bisection(a, b);
    return 0;
}

#include<iostream>
#include<iomanip>
#include<math.h>

#define f(x) cos(x) - x * exp(x)

using namespace std;

int main()
{
	 float x0, x1, x, f0, f1, f, e;
	 int step = 1;

     cout<< setprecision(6)<< fixed;

	 up:
	 cout<<"Enter first guess: ";//a = 0
     cin>>x0;
     cout<<"Enter second guess: ";//b =1
     cin>>x1;
     cout<<"Enter tolerable error: ";// e=0.00001
     cin>>e;

	 f0 = f(x0);
	 f1 = f(x1);

	 if( f0 * f1 > 0.0)
	 {
		  cout<<"Incorrect Initial Guesses."<< endl;
		  goto up;
	 }
     cout<< endl<<"****************"<< endl;
	 cout<<"Bisection Method"<< endl;
	 cout<<"****************"<< endl;
	 do
	 {
		  x = (x0 + x1)/2;
		  f = f(x);

		  cout<<"Iteration-"<< step<<":\t x = "<< setw(10)<< x<<" and f(x) = "<< setw(10)<< f(x)<< endl;

		  if( f0 * f < 0)
		  {
			   x1 = x;
		  }
		  else
		  {
			   x0 = x;
		  }
		  step = step + 1;
	 }while(fabs(f)>e);

	 cout<< endl<< "Root is: "<<  x<< endl;

	 return 0;
}


//false

#include<bits/stdc++.h>
using namespace std;
#define MAX_ITER 1000000

double func(double x)
{
    return x*x*x - x*x + 2;
}

void regulaFalsi(double a, double b)
{
    if (func(a) * func(b) >= 0)
    {
        cout << "You have not assumed right a and b\n";
        return;
    }

    double c = a;

    for (int i=0; i < MAX_ITER; i++)
    {
        c = (a*func(b) - b*func(a))/ (func(b) - func(a));
        if (func(c)==0)
            break;

        else if (func(c)*func(a) < 0)
            b = c;
        else
            a = c;
    }
    cout << "The value of root is : " << c;
}

int main()
{
    double a =-200, b = 300;
    regulaFalsi(a, b);
    return 0;
}

#include<iostream>
#include<iomanip>
#include<math.h>

#define f(x) cos(x) - x * exp(x)

using namespace std;

int main()
{
	 float x0, x1, x, f0, f1, f, e;
	 int step = 1;

     cout<< setprecision(6)<< fixed;

	 up:
	 cout<<"Enter first guess: ";
     cin>>x0;
     cout<<"Enter second guess: ";
     cin>>x1;
     cout<<"Enter tolerable error: ";
     cin>>e;

	 f0 = f(x0);
	 f1 = f(x1);

	 if( f0 * f1 > 0.0)
	 {
		  cout<<"Incorrect Initial Guesses."<< endl;
		  goto up;
	 }
     cout<< endl<<"*********************"<< endl;
	 cout<<"False Position Method"<< endl;
	 cout<<"*********************"<< endl;
	 do
	 {
		  x = x0 - (x0-x1) * f0/ (f0-f1);
		  f = f(x);

		  cout<<"Iteration-"<< step<<":\t x = "<< setw(10)<< x<<" and f(x) = "<< setw(10)<< f(x)<< endl;

		  if( f0 * f < 0)
		  {
			   x1 = x;
			   f1 = f;
		  }
		  else
		  {
			   x0 = x;
			   f0 = f;
		  }
		  step = step + 1;
	 }while(fabs(f)>e);

	 cout<< endl<<"Root is: "<< x<< endl;

	 return 0;
}

//newton raphson

#include<bits/stdc++.h>
#define EPSILON 0.001
using namespace std;

double func(double x)
{
    return x*x*x - x*x + 2;
}

double derivFunc(double x)
{
    return 3*x*x - 2*x;
}

void newtonRaphson(double x)
{
    double h = func(x) / derivFunc(x);
    while (abs(h) >= EPSILON)
    {
        h = func(x)/derivFunc(x);

        // x(i+1) = x(i) - f(x) / f'(x)
        x = x - h;
    }

    cout << "The value of the root is : " << x;
}

int main()
{
    double x0 = -20;
    newtonRaphson(x0);
    return 0;
}


#include<iostream>
#include<iomanip>
#include<math.h>
#include<stdlib.h>


#define    f(x)    3*x - cos(x) -1

#define   g(x)   3 + sin(x)

using namespace std;

int main()
{
	 float x0, x1, f0, f1, g0, e;
	 int step = 1, N;

     cout<< setprecision(6)<< fixed;

	 cout<<"Enter initial guess: ";
	 cin>>x0;
	 cout<<"Enter tolerable error: ";
	 cin>>e;
	 cout<<"Enter maximum iteration: ";
	 cin>>N;

	 cout<< endl<<"*********************"<< endl;
	 cout<<"Newton Raphson Method"<< endl;
	 cout<<"*********************"<< endl;
	 do
	 {
		  g0 = g(x0);
		  f0 = f(x0);
		  if(g0 == 0.0)
		  {
			   cout<<"Mathematical Error.";
			   exit(0);
		  }


		  x1 = x0 - f0/g0;


		  cout<<"Iteration-"<< step<<":\t x = "<< setw(10)<< x1<<" and f(x) = "<< setw(10)<< f(x1)<< endl;
		  x0 = x1;

		  step = step+1;

		  if(step > N)
		  {
			   cout<<"Not Convergent.";
			   exit(0);
		  }

		  f1 = f(x1);

	 }while(fabs(f1)>e);

	 cout<< endl<<"Root is: "<< x1;
	 return 0;
}


//secant


#include <bits/stdc++.h>
using namespace std;

float f(float x)
{
    float f = pow(x, 3) + x - 1;
    return f;
}

void secant(float x1, float x2, float E)
{
    float n = 0, xm, x0, c;
    if (f(x1) * f(x2) < 0) {
        do {
            x0 = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));

            c = f(x1) * f(x0);

            x1 = x2;
            x2 = x0;

            n++;

            if (c == 0)
                break;
            xm = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
        } while (fabs(xm - x0) >= E);

        cout << "Root of the given equation=" << x0 << endl;
        cout << "No. of iterations = " << n << endl;
    } else
        cout << "Can not find a root in the given interval";
}


int main()
{
    float x1 = 0, x2 = 1, E = 0.0001;
    secant(x1, x2, E);
    return 0;
}


#include<iostream>
#include<iomanip>
#include<math.h>
#include<stdlib.h>


#define    f(x)    x*x*x - 2*x - 5

using namespace std;

int main()
{
	 float x0, x1, x2, f0, f1, f2, e;
	 int step = 1, N;

   cout<< setprecision(6)<< fixed;


	 cout<<"Enter first guess: ";
	 cin>>x0;
	 cout<<"Enter second guess: ";
	 cin>>x1;
	 cout<<"Enter tolerable error: ";
	 cin>>e;
	 cout<<"Enter maximum iteration: ";
	 cin>>N;

     cout<< endl<<"**************"<< endl;
	 cout<<"Secant Method"<< endl;
	 cout<<"**************"<< endl;
	 do
	 {
		  f0 = f(x0);
		  f1 = f(x1);
		  if(f0 == f1)
		  {
			   cout<<"Mathematical Error.";
			   exit(0);
		  }

		  x2 = x1 - (x1 - x0) * f1/(f1-f0);
		  f2 = f(x2);

		  cout<<"Iteration-"<< step<<":\t x2 = "<< setw(10)<< x2<<" and f(x2) = "<< setw(10)<< f(x2)<< endl;

		  x0 = x1;
		  f0 = f1;
		  x1 = x2;
		  f1 = f2;

		  step = step + 1;

		  if(step > N)
		  {
			   cout<<"Not Convergent.";
			   exit(0);
		  }
	 }while(fabs(f2)>e);

	 cout<< endl<<"Root is: "<< x2;

	 return 0;
}





//forward

#include <bits/stdc++.h>
using namespace std;

float u_cal(float u, int n)
{
    float temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u - i);
    return temp;
}

int fact(int n)
{
    int f = 1;
    for (int i = 2; i <= n; i++)
        f =f *i;
    return f;
}

int main()
{
    int n = 5;
    float x[] = { 0.10,0.15,0.20,0.25,0.30};

    float y[n][n];
    y[0][0] = 0.1003;
    y[1][0] = 0.1511;
    y[2][0] = 0.2027;
    y[3][0] = 0.2553;
    y[4][0] = 0.3093;


    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++)
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
    }

    for (int i = 0; i < n; i++) {
        cout << setw(4) << x[i]
             << "\t";
        for (int j = 0; j < n - i; j++)
            cout << setw(4) << y[i][j]
                 << "\t";
        cout << endl;
    }

    float value = 52;

    float sum = y[0][0];
    float u = (value - x[0]) / (x[1] - x[0]);
    for (int i = 1; i < n; i++) {
        sum = sum + (u_cal(u, i) * y[0][i]) /
                                 fact(i);
    }

    cout << "\n Value at " << value << " is "
         << sum << endl;
    return 0;
}



#include<iostream>

using namespace std;

int main()
{
 float x[20], y[20][20];
 int i,j, n;

 /* Input Section */
 cout << "Enter number of data? " << endl;
 cin >> n;

 cout << "Enter data: " << endl;
 for(i = 0; i < n ; i++)
 {
  cout << "x[" << i << "] = ";
  cin >> x[i];
  cout << "y[" << i <<"] = ";
  cin >> y[i][0];
 }

 /* Generating Forward Difference Table */
 for(i = 1; i < n; i++)
 {
  for(j = 0; j < n-i; j++)
  {
   y[j][i] = y[j+1][i-1] - y[j][i-1];
  }
 }

 /* Displaying Forward Difference Table */
 cout << endl << "FORWARD DIFFERENCE TABLE" << endl;
 for(i = 0; i < n; i++)
 {
  cout << x[i];
  for(j = 0; j < n-i ; j++)
  {
   cout << "\t" << y[i][j];
  }
  cout << endl;
 }

 return 0;
}


//backward

#include <bits/stdc++.h>
using namespace std;

float u_cal(float u, int n)
{
    float temp = u;
    for (int i = 1; i < n; i++)
        temp = temp * (u + i);
    return temp;
}

int fact(int n)
{
    int f = 1;
    for (int i = 2; i <= n; i++)
        f = f*i;
    return f;
}

int main()
{
    int n = 5;
    float x[] = { 1891, 1901, 1911,
                  1921, 1931 };

    float y[n][n];
    y[0][0] = 46;
    y[1][0] = 66;
    y[2][0] = 81;
    y[3][0] = 93;
    y[4][0] = 101;

    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--)
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++)
            cout << setw(4) << y[i][j]
                 << "\t";
        cout << endl;
    }

    float value = 1925;

    float sum = y[n - 1][0];
    float u = (value - x[n - 1]) / (x[1] - x[0]);
    for (int i = 1; i < n; i++) {
        sum = sum + (u_cal(u, i) * y[n - 1][i]) /
                                     fact(i);
    }

    cout << "\n Value at " << value << " is "
         << sum << endl;
    return 0;
}


#include<iostream>

using namespace std;

int main()
{
 float x[20], y[20][20];
 int i,j, n;

 /* Input Section */
 cout << "Enter number of data? " << endl;
 cin >> n;

 cout << "Enter data: " << endl;
 for(i = 0; i < n ; i++)
 {
  cout << "x[" << i << "] = ";
  cin >> x[i];
  cout << "y[" << i <<"] = ";
  cin >> y[i][0];
 }

 /* Generating Backward Difference Table */
 for(i = 1; i < n; i++)
 {
  for(j = n-1; j > i-1; j--)
  {
   y[j][i] = y[j][i-1] - y[j-1][i-1];
  }
 }

 /* Displaying Backward Difference Table */
 cout << endl << "BACKWARD DIFFERENCE TABLE" << endl;

 for(i = 0; i < n; i++)
 {
  cout << x[i];
  for(j = 0; j <= i ; j++)
  {
   cout << "\t" << y[i][j];
  }
  cout << endl;
 }

 return 0;
}


//lagrange

#include<bits/stdc++.h>
using namespace std;

struct Data
{
    int x, y;
};

double interpolate(Data f[], int xi, int n)
{
    double result = 0;

    for (int i=0; i<n; i++)
    {
        double term = f[i].y;
        for (int j=0;j<n;j++)
        {
            if (j!=i)
                term = term*(xi - f[j].x)/double(f[i].x - f[j].x);
        }

        result += term;
    }

    return result;
}

int main()
{
    Data f[] = {{0,2}, {1,3}, {2,12}, {5,147}};

    cout << "Value of f(3) is : " << interpolate(f, 3, 5);

    return 0;
}


#include<iostream>
#include<conio.h>

using namespace std;

int main()
{
	 float x[100], y[100], xp, yp=0, p;
	 int i,j,n;

	 /* Input Section */
	 cout<<"Enter number of data: ";
	 cin>>n;
	 cout<<"Enter data:"<< endl;
	 for(i=1;i<=n;i++)
	 {
		  cout<<"x["<< i<<"] = ";
		  cin>>x[i];
		  cout<<"y["<< i<<"] = ";
		  cin>>y[i];
	 }
	 cout<<"Enter interpolation point: ";
	 cin>>xp;

	 /* Implementing Lagrange Interpolation */
	 for(i=1;i<=n;i++)
	 {
		  p=1;
		  for(j=1;j<=n;j++)
		  {
			   if(i!=j)
			   {
			    	p = p* (xp - x[j])/(x[i] - x[j]);
			   }
		  }
		  yp = yp + p * y[i];
	 }
	 cout<< endl<<"Interpolated value at "<< xp<< " is "<< yp;

	 return 0;
}


//simsom1/3


#include <iostream>
#include <math.h>
using namespace std;

float func(float x)
{
    return log(x);
}

float simpsons_(float ll, float ul, int n)
{
    float h = (ul - ll) / n;

    float x[10], fx[10];

    for (int i = 0; i <= n; i++) {
        x[i] = ll + i * h;
        fx[i] = func(x[i]);
    }

    float res = 0;
    for (int i = 0; i <= n; i++) {
        if (i == 0 || i == n)
            res += fx[i];
        else if (i % 2 != 0)
            res += 4 * fx[i];
        else
            res += 2 * fx[i];
    }
    res = res * (h / 3);
    return res;
}

int main()
{
    float lower_limit = 4; // Lower limit
    float upper_limit = 5.2; // Upper limit
    int n = 6; // Number of interval
    cout << "The  output is : " <<  simpsons_(lower_limit, upper_limit, n);
    return 0;
}


//tropizoidal

#include<iostream>
#include<math.h>

/* Define function here */
#define f(x) 1/(1+pow(x,2))

using namespace std;
int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;

 /* Input */
 cout<<"Enter lower limit of integration: ";//0
 cin>>lower;
 cout<<"Enter upper limit of integration: ";//6
 cin>>upper;
 cout<<"Enter number of sub intervals: ";//6
 cin>>subInterval;

 /* Calculation */

 /* Finding step size */
 stepSize = (upper - lower)/subInterval;

 /* Finding Integration Value */
 integration = f(lower) + f(upper);

 for(i=1; i<= subInterval-1; i++)
 {
  k = lower + i*stepSize;
  integration = integration + 2 * (f(k));
 }

 integration = integration * stepSize/2;

 cout<< endl<<"Required value of integration is: "<< integration;

 return 0;
}


//euler

#include <iostream>
using namespace std;

float func(float x, float y)
{
    return (x + y + x * y);
}

void euler(float x0, float y, float h, float x)
{
    float temp = -0;

    while (x0 < x) {
        temp = y;
        y = y + h * func(x0, y);
        x0 = x0 + h;
    }

    cout << "Approximate solution at x = "
         << x << "  is  " << y << endl;
}

int main()
{
    float x0 = 0;
    float y0 = 1;
    float h = 0.025;

    float x = 0.1;

    euler(x0, y0, h, x);
    return 0;
}


//runge kutta


#include<iostream>

#define f(x,y) (y*y-x*x)/(y*y+x*x)

using namespace std;

#define f(x,y) (y*y-x*x)/(y*y+x*x)

using namespace std;
int main()
{
 float x0, y0, xn, h, yn, k1, k2, k3, k4, k;
 int i, n;

 cout<<"Enter Initial Condition"<< endl;
 cout<<"x0 = ";
 cin>> x0;
 cout<<"y0 = ";
 cin >> y0;
 cout<<"Enter calculation point xn = ";
 cin>>xn;
 cout<<"Enter number of steps: ";
 cin>> n;

 h = (xn-x0)/n;

 cout<<"\nx0\ty0\tyn\n";
 cout<<"------------------\n";
 for(i=0; i < n; i++)
 {
  k1 = h * (f(x0, y0));
  k2 = h * (f((x0+h/2), (y0+k1/2)));
  k3 = h * (f((x0+h/2), (y0+k2/2)));
  k4 = h * (f((x0+h), (y0+k3)));
  k = (k1+2*k2+2*k3+k4)/6;
  yn = y0 + k;
  cout<< x0<<"\t"<< y0<<"\t"<< yn<< endl;
  x0 = x0+h;
  y0 = yn;
 }

 cout<<"\nValue of y at x = "<< xn<< " is " << yn;

 return 0;
}



/*
Enter Initial Condition
x0 = 0
y0 = 1
Enter calculation point xn = 0.6
Enter number of steps: 3

x0      y0      yn
------------------
0       1       1.196
0.2     1.196   1.37527
0.4     1.37527 1.53311

Value of y at x = 0.6 is 1.533
*/
