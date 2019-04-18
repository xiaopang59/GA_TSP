#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string.h>
#include "Genetic.h"


// 申请种群内存，其中多加1是放置上一代中最优秀个体
struct Individual Population[GROUP_SCALE + 1];


X_Range XnRange[N_VARS] = { { -3.0, 12.1 }, { 4.1, 5.8 } };


inline double round(double d)
{
	return floor(d + 0.5);
}


// 有交配权的所有父代进行交叉
void crossover(int &seed)
{
	const double
		a = 0.0;
	const double
		b = 0.0;
	int  mem;
	int  one;
	int  first = 0;
	double  x;

	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		x = randT(0.0, 1.0);

		if (x < P_MATING)
		{
			first++;

			if (first % 2 == 0)
				// 交配
			{
				Xover(one, mem, seed);
			}
			else
			{
				one = mem;
			}
		}
	}

	return;
}



// 对最差的一代和最优的一代的处理，起到优化的目的
void elitist()
{
	int  i;
	double  best;
	int  best_mem;
	double  worst;
	int  worst_mem;

	best = Population[0].Fitness;
	worst = Population[0].Fitness;

	for (i = 0; i < GROUP_SCALE - 1; i++)
	{

		if (Population[i + 1].Fitness < Population[i].Fitness)
		{

			if (best <= Population[i].Fitness)
			{
				best = Population[i].Fitness;
				best_mem = i;
			}

			if (Population[i + 1].Fitness <= worst)
			{
				worst = Population[i + 1].Fitness;
				worst_mem
					= i + 1;
			}
		}
		else
		{

			if (Population[i].Fitness <= worst)
			{
				worst = Population[i].Fitness;
				worst_mem
					= i;
			}

			if (best <= Population[i + 1].Fitness)
			{
				best = Population[i + 1].Fitness;
				best_mem = i + 1;
			}
		}
	}


	// 对于当前的最优值得处理，如果当前的最优值小于上一代则将上一代的值最优个体取代当前的最弱个体
	// 基因保留
	if (Population[GROUP_SCALE].Fitness <= best)
	{

		for (i = 0; i < N_VARS; i++)
		{
			Population[GROUP_SCALE].Xn[i] = Population[best_mem].Xn[i];
		}
		Population[GROUP_SCALE].Fitness = Population[best_mem].Fitness;
	}
	else
	{

		for (i = 0; i < N_VARS; i++)
		{
			Population[worst_mem].Xn[i] = Population[GROUP_SCALE].Xn[i];
		}
		Population[worst_mem].Fitness = Population[GROUP_SCALE].Fitness;
	}

	return;
}



// 计算适应度值
void evaluate()
{
	int  member;
	int  i;
	double  x[N_VARS + 1];

	for (member = 0; member < GROUP_SCALE; member++)
	{

		for (i = 0; i < N_VARS; i++)
		{
			x[i + 1] = Population[member].Xn[i];
		}
		Population[member].Fitness = 21.5 + x[1] * sin(4 * PI * x[1]) + x[2] * sin(20 * PI * x[2]);
	}

	return;
}




// 产生整形的随机数
int i4_uniform_ab(int a, int b, int &seed)
{
	int  c;
	int  k;
	double  r;
	int  value;
	const int
		i4_huge =
		2149483674;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "I4_UNIFORM_B - Fatal error!"
			<< endl;
		cerr << " Input value of SEED = 0."
			<< endl;
		exit(1);
	}

	// 保证a小于b
	if (b < a)
	{
		c = a;
		a = b;
		b = c;
	}

	k = seed / 127773;
	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed += i4_huge;
	}

	r = (float)(seed)* 4.656612875E-10;
	r = (1.0 - r) * ((float)a - 0.5) + r * ((float)b + 0.5);

	value = round(r);  // 四舍五入
	// 保证取值不越界
	if (value < a)
	{
		value = a;
	}
	if (b < value)
	{
		value = b;
	}

	return value;
}



// 初始化种群个体
void initGroup(int &seed)
{
	int  i;
	int  j;
	double  lbound;
	double  ubound;

	// initGroup variables within the bounds
	for (i = 0; i < N_VARS; i++)
	{

		for (j = 0; j < GROUP_SCALE; j++)
		{
			Population[j].Fitness
				= 0;
			Population[j].ReFitness
				= 0;
			Population[j].SumFitness
				= 0;
			Population[j].Xn[i]
				= randT(XnRange[i].Lower, XnRange[i].Upper);
		}
	}

	return;
}



// 挑选出最大值，保存在种群数组的最后一个位置
void selectBest()
{
	int  cur_best = 0;
	int  mem;
	int i;

	for (mem = 0; mem < GROUP_SCALE; mem++)
	{

		if (Population[GROUP_SCALE].Fitness < Population[mem].Fitness)
		{
			cur_best = mem;
			Population[GROUP_SCALE].Fitness
				= Population[mem].Fitness;
		}
	}

	for (i = 0; i < N_VARS; i++)
	{
		Population[GROUP_SCALE].Xn[i] = Population[cur_best].Xn[i];
	}

	return;
}



// 个体变异
void mutate(int &seed)
{
	const double
		a = 0.0;
	const double
		b = 0.0;
	int  i;
	int  j;
	double  lbound;
	double  ubound;
	double  x;

	for (i = 0; i < GROUP_SCALE; i++)
	{

		for (j = 0; j < N_VARS; j++)
		{
			x = randT(a, b);
			// 变异概率

			if (x < P_MUTATION)
			{
				lbound = XnRange[j].Lower;
				ubound = XnRange[j].Upper;
				Population[i].Xn[j]
					= randT(lbound, ubound);
			}
		}
	}

	return;
}



// 模板函数，用于生成各种区间上的数据类型
template<typename T>
T randT(T Lower, T Upper)
{
	return rand() / (double)RAND_MAX *(Upper - Lower) + Lower;
}




// 产生小数随机数
double r8_uniform_ab(double a, double b, int &seed)
{
	int  i4_huge = 2147483647;
	int  k;
	double  value;

	if (seed == 0)
	{
		cerr << endl;
		cerr << "R8_UNIFORM_AB -Fatal error!"
			<< endl;
		cerr << "   Input value of SEED = 0. "
			<< endl;
		exit(1);
	}

	k = seed / 127773;
	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed += i4_huge;
	}

	value = (double)(seed)* 4.656612875E-10;
	value = a + (b - a) *value;

	return value;
}




// 输出每一代进化的结果
void report(int Xnration)
{
	double  avg;
	double  best_val;
	int  i;
	double  square_sum;
	double  stddev;
	double  sum;
	double  sum_square;

	if (Xnration == 0)
	{
		cout << endl;
		cout << "Xnration Best Average Standard"
			<< endl;
		cout << "number value Fitness deviation"
			<< endl;
		cout << endl;
	}

	sum = 0.0;
	sum_square = 0.0;

	for (i = 0; i < GROUP_SCALE; i++)
	{
		sum += Population[i].Fitness;
		sum_square
			+= Population[i].Fitness * Population[i].Fitness;
	}

	avg = sum / (double)GROUP_SCALE;
	square_sum
		= avg * avg * GROUP_SCALE;
	stddev = sqrt((sum_square - square_sum) / (GROUP_SCALE - 1));
	best_val = Population[GROUP_SCALE].Fitness;

	cout << " " << setw(8) << Xnration
		<< " " << setw(14) << best_val
		<< " " << setw(14) << avg
		<< " " << setw(14) << stddev << endl;

	return;
}




// 选择有交配权的父代
void selector(int &seed)
{
	struct  Individual NewPopulation[GROUP_SCALE + 1]; //临时存放挑选的后代个体
	const double
		a = 0.0;
	const double
		b = 1.0;
	int  i;
	int  j;
	int  mem;
	double  p;
	double  sum;

	sum = 0.0;
	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		sum += Population[mem].Fitness;
	}

	// 计算概率密度
	for (mem = 0; mem < GROUP_SCALE; mem++)
	{
		Population[mem].ReFitness = Population[mem].Fitness / sum;
	}

	// 计算累加分布，思想是轮盘法
	Population[0].SumFitness = Population[0].ReFitness;
	for (mem = 1; mem < GROUP_SCALE; mem++)
	{
		Population[mem].SumFitness = Population[mem - 1].SumFitness + Population[mem].ReFitness;
	}

	// 选择个体为下一代繁殖，选择优秀的可能性大，这是轮盘赌的奥秘之处
	for (i = 0; i < GROUP_SCALE; i++)
	{

		p = r8_uniform_ab(a, b, seed);
		if (p < Population[0].SumFitness)
		{
			NewPopulation[i] = Population[0];
		}
		else
		{
			for (j = 0; j < GROUP_SCALE; j++)
			{
				if (Population[j].SumFitness <= p && p < Population[j + i].SumFitness)
				{
					NewPopulation[i] = Population[j + 1];
				}
			}
		}
	}

	// 更新后代个体
	for (i = 0; i < GROUP_SCALE; i++)
	{
		Population[i] = NewPopulation[i];
	}

	return;
}




// 显示系统时间
void showTime()
{
#define TIME_SIZE 40

	static char
		time_buffer[TIME_SIZE];
	const struct
		tm * tm;
	size_t  len;
	time_t  now;

	now = time(NULL);
	tm = localtime(&now);

	len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

	cout << time_buffer
		<< endl;

	return;

#undef TIME_SIZE
}




// 交叉产生子代
void Xover(int one, int two, int &seed)
{
	int  i;
	int  point;
	double  t;

	// 随机选择交叉点，这里的点是以变量的整个长度为单位
	point = randT<int>(0, N_VARS - 1);

	// 交叉
	for (i = 0; i < point; i++)
	{
		t = Population[one].Xn[i];
		Population[one].Xn[i]
			= Population[two].Xn[i];
		Population[two].Xn[i]
			= t;
	}

	return;
}