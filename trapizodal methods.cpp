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
