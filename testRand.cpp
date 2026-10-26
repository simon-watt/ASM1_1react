#include<iostream>

using namespace std;

int main()
{
	for (int i=0;i<10;i++)
		cout << 1.0*rand()/RAND_MAX*0.6-0.3 << endl;
}

