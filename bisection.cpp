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
