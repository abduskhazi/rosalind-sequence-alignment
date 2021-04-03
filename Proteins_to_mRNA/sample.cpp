#include<iostream>

using namespace std;

int main()
{
    auto x = 2;
    auto modx = 2;
    for(auto i = 1; i < 15; ++i)
    {
        cout << i << endl;
        x = x * i;
        modx = (modx * i ) % 100;
        cout << "x: " << x << endl;
        cout << "modxL " << modx << endl;
        cout << endl;
    }

    return 0;
}

