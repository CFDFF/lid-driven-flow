//顶盖驱动流的格子Boltzmann模拟


#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

using namespace std;



const int Q = 9;               //D2Q9模型
const int NX = 256;         //x方向网格数量
const int NY = 256;         //y方向
const double U = 0.1;      //顶盖流速





int e[Q][2] = { {0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1} };              //离散速度的设置
double w[Q] = { 4.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 36,1.0 / 36,1.0 / 36,1.0 / 36 };          //权系数？
double rho[NX + 1][NY + 1],            //密度场
u[NX + 1][NY + 1][2], u0[NX + 1][NY + 1][2],       //  前后两个时层       速度场 
f[NX + 1][NY + 1][Q], F[NX + 1][NY + 1][Q];        //演化前后       两个分布函数
int i, j, k, ip, jp, n;
double c, Re, dx, dy, Lx, Ly, dt, rho0, P0, tau_f, niu, error;


//函数声明
void init();                                                                      //初始化流场
double feq(int k, double rho, double u[2]);
void evolution();
void output(int m);
void Error();

int main()
{
	using namespace std;
	init();
	for (n = 0;; n++)
	{
		evolution();

		if (n % 100 == 0)
		{
			Error();
			cout << "The" << n << "th computation result: " << endl << "The u,v of point(NX/2, NY/2) is " << setprecision(6) << u[NX / 2][NY / 2][0] << "," << u[NX / 2][NY / 2][1] << endl;
			cout << "The max relative error of uv is: " << setiosflags(ios::scientific) << error << endl;
			if (n >= 1000)
			{
				if (n % 1000 == 0)
				{
					output(n);
				}
				if (error < 1.0e-6)
				{
					break;
				}
			}
		}
	}
	return 0;
}


void init()
{
	//初始化参数
	dx = 1.0;
	dy = 1.0;
	Lx = dx * (double)NX;
	Ly = dy * (double)NY;
	dt = dx;
	c = dx / dt;
	rho0 = 1.0;
	Re = 1000;
	niu = U * Lx / Re;
	tau_f = 3.0 * niu + 0.5;
	cout << "tau_f=  " << tau_f << endl;


	//初始化流场
	for (i = 0; i <= NX; i++)
	{
		for (j = 0; j <= NY; j++)
		{
			u[i][j][0] = 0;
			u[i][j][1] = 0;
			rho[i][j] = rho0;
			u[i][NY][0] = U;
			for (k = 0; k < Q; k++)
			{
				f[i][j][k] = feq(k, rho[i][j], u[i][j]);
			}
		}
	}
}

//计算平衡态分布函数
double feq(int k, double rho, double u[2])
{
	double eu, uv, feq;
	eu = (e[k][0] * u[0] + e[k][1] * u[1]);
	uv = (u[0] * u[0] + u[1] * u[1]);
	feq = w[k] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
	return feq;
}

//分布函数的演化
void evolution()
{
	for (i = 1; i < NX; i++)
	{
		for (j = 1; j < NY; j++)
		{
			for (k = 0; k < Q; k++)
			{
				ip = i - e[k][0];
				jp = j - e[k][1];
				F[i][j][k] = f[ip][jp][k] + (feq(k, rho[ip][jp], u[ip][jp]) - f[ip][jp][k]) / tau_f;
			}
		}
	}

	//计算宏观物理量
	for (i = 1; i < NX; i++)
	{
		for (j = 1; j < NY; j++)
		{
			u0[i][j][0] = u[i][j][0];
			u0[i][j][1] = u[i][j][1];
			rho[i][j] = 0;
			u[i][j][0] = 0;
			u[i][j][1] = 0;
			for (k = 0; k < Q; k++)
			{
				f[i][j][k] = F[i][j][k];
				rho[i][j] += f[i][j][k];
				u[i][j][0] += e[k][0] * f[i][j][k];
				u[i][j][1] += e[k][1] * f[i][j][k];
			}
			u[i][j][0] /= rho[i][j];
			u[i][j][1] /= rho[i][j];
		}
	}




	//边界处理
	for (j = 1; j < NY; j++)
	{
		for (k = 0; k < Q; k++)
		{
			//左壁面
			rho[NX][j] = rho[NX - 1][j];
			f[NX][j][k] = feq(k, rho[NX][j], u[NX][j]) + f[NX - 1][j][k] - feq(k, rho[NX - 1][j], u[NX - 1][j]);
			//右壁面
			rho[0][j] = rho[1][j];
			f[0][j][k] = feq(k, rho[0][j], u[0][j]) + f[1][j][k] - feq(k, rho[1][j], u[1][j]);

		}
	}
	for (i = 1; i < NX; i++)
	{
		for (k = 0; k < Q; k++)
		{
			//下壁面
			rho[i][0] = rho[i][1];
			f[i][0][k] = feq(k, rho[i][0], u[i][0]) + f[i][1][k] - feq(k, rho[i][1], u[i][1]);

			//上面
			rho[i][NY] = rho[i][NY - 1];
			u[i][NY][0] = U;
			f[i][NY][k] = feq(k, rho[i][NY], u[i][NY]) + f[i][NY - 1][k] - feq(k, rho[i][NY - 1], u[i][NY - 1]);
		}
	}
}



void Error()
{
	double temple1, temple2;
	temple1 = 0;
	temple2 = 0;
	for (i = 1; i < NX; i++)
	{
		for (j = 1; j < NY; j++)
		{
			temple1 += ((u[i][j][0] - u0[i][j][0]) * (u[i][j][0] - u0[i][j][0]) + (u[i][j][1] - u0[i][j][1]) * (u[i][j][1] - u0[i][j][1]));
			temple2 += ((u[i][j][0] * u[i][j][0]) + (u[i][j][1] * u[i][j][1]));
		}	
	}
		temple1 = sqrt(temple1);
		temple2 = sqrt(temple2);
		error = temple1 / (temple2 + 1e-30);

}


void output(int m)
{
	ostringstream name;
	name << "cavity_" << m << ".dat";
	ofstream out(name.str().c_str());

	out << "Title= \"    LBM Lid Driven Flow \"\n" << "VARIABLES= \" X\", \" Y\", \" U\", \" V\"\n" << "ZONE T=\"BOX\",I= " << NX + 1 << ",J= " << NY + 1 << ",F=POINT" << endl;
	for (j = 0; j <= NY; j++)
	{
		for (i = 0; i <= NX; i++)
		{
			out << i << "  " << j << "  " << u[i][j][0] << "  " << u[i][j][1] << endl;
		}
	}
}

