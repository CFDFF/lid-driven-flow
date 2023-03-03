#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

using namespace std;


void Show(int a[1][2])
{

	cout << a[0] [0]<< endl;
	cout << a[0][1] << endl;



}

void test1()
{

	int m = 3;
	ostringstream name;
	name << "cavity_" << m<< ".dat";
	cout << name.str().c_str() << endl;




}

int main1()
{
	test1();



	//int a[2][2][2] = 
	//{
	//	21, 12,    34, 24,
	//	87, 6,     3, 4,
	// };

	//Show(a[1]);

	////cout << sizeof(a) << endl;
	////cout << sizeof(a[0][0][0]) << endl;





	return 0;

}