#include<iostream>

using namespace std;

int main()
{
	long long i=1;

	for (int j=0;j<20;j++)
	{
		i*=61;
		cout << j << " " << i << endl;
	}


}
