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
