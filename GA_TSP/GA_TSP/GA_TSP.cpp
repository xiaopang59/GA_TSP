#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>
#include <time.h>

using namespace std;

const int		cities	= 50;	// 城市个数
const int		MAXX	= 100;	// 迭代次数
const double	pc		= 0.8;	// 交配概率
const double	pm		= 0.05;	// 变异概率
const int		num		= 50;	// 种群的大小

int bestsolution;
// 最优染色体
double distance1[cities][cities];
// 城市之间的距离

struct _group
{
	int  city[cities]; // 城市的顺序
	double  adapt; // 适应度
	double  p; // 在种群中的幸存概率
}group[num], grouptemp[num];


//struct point
//{
// double	x, y;
//}Cpoint[cities];

// 随机产生cities个城市之间的相互距离
void init()
{
	int i, j;
	memset(distance1, 0, sizeof(distance1));
	srand((unsigned)time(NULL));
	for (i = 0; i < cities; i++)
	{
		for (j = i + 1; j < cities; j++)
		{
			distance1[i][j] = rand() % 100;
			distance1[j][i] = distance1[i][j];
		}
	}

	// 打印距离矩阵
	cout << "城市的距离矩阵如下: " << endl;
	for (i = 0; i < cities; i++)
	{
		for (j = 0; j < cities; j++)
		{
			cout << distance1[i][j] << "\t";
		}
		cout << endl;
	}
}

// 随机产生初始种群
void groupproduce()
{
	int i, j, t, k, flag;

	// 初始化
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < cities; j++)
		{
			group[i].city[j] = -1;
		}
	}

	srand((unsigned)time(NULL));
	for (i = 0; i < num; i++)
	{
		// 产生10个不相同的数字
		for (j = 0; j < cities;)
		{
			t = rand() % cities;
			flag = 1;
			for (k = 0; k < j; k++)
			{
				if (group[i].city[k] == t)
				{
					flag = 0;
					break;
				}
			}
			if (flag)
			{
				group[i].city[j] = t;
				j++;
			}
		}
	}

	// 打印种群基因
	cout << "初始的种群" << endl;
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < cities; j++)
		{
			cout << group[i].city[j] << " ";
		}
		cout << endl;
	}
}

// 评价函数，找出最优染色体
void Fitness()
{
	int  i, j;
	int  n1, n2;
	double  sumdistance;
	int  biggestsum = 0;
	double  biggestp = 0;

	for (i = 0; i < num; i++)
	{
		sumdistance = 0;
		for (j = 1; j < cities; j++)
		{
			n1 = group[i].city[j - 1];
			n2 = group[i].city[j];
			sumdistance += distance1[n1][n2];
		}
		sumdistance += distance1[group[i].city[cities - 1]][group[i].city[0]];
		group[i].adapt = sumdistance;
		// 每条染色体的路径总和
		biggestsum += sumdistance;
		// 种群的总路径
	}

	// 计算染色体的幸存能力，路径越短生存概率越大
	for (i = 0; i < num; i++)
	{
		group[i].p = 1 - (double)group[i].adapt / (double)biggestsum;
		biggestp += group[i].p;
	}

	for (i = 0; i < num; i++)
	{
		group[i].p /= biggestp;
	}// 在种群中的幸存概率，总和为1
	// 求最佳路径
	bestsolution = 0;
	for (i = 0; i < num; i++)
	{
		if (group[i].p > group[bestsolution].p)
		{
			bestsolution = i;
		}
	}
}

// 选择
void Select()
{
	int  i, j, temp;
	double  gradient_p[num]; // 梯度概率
	double  select_p[num]; // 选择染色体的随机概率
	int  select[num];

	// 初始化梯度概率
	for (i = 0; i < num; i++)
	{
		gradient_p[i]
			= 0.0;
		select_p[i]
			= 0.0;
	}
	gradient_p[0] = group[0].p;

	for (i = 1; i < num; i++)
	{
		gradient_p[i] = gradient_p[i - 1] + group[i].p;
	}

	srand((unsigned)time(NULL));

	// 随机产生染色体的存活概率
	for (i = 0; i < num; i++)
	{
		select_p[i] = rand() % 100;
		select_p[i] /= 100;
	}

	// 选择能生存的染色体
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < num; j++)
		{
			if (select_p[i] < gradient_p[j])
			{
				select[i] = j;
				// 第i个位置存放第j个染色体
				break;
			}
		}
	}

	// 拷贝种群
	for (i = 0; i < num; i++)
	{
		grouptemp[i].adapt
			= group[i].adapt;
		grouptemp[i].p
			= group[i].p;
		for (j = 0; j < cities; j++)
		{
			grouptemp[i].city[j] = group[i].city[j];
		}
	}

	// 数据更新
	for (i = 0; i < num; i++)
	{
		temp = select[i];
		group[i].adapt
			= grouptemp[temp].adapt;
		group[i].p
			= grouptemp[temp].p;
		for (j = 0; j < cities; j++)
		{
			group[i].city[j] = grouptemp[temp].city[j];
		}
	}
}

// 交配，对每个染色体产生交配概率，满足交配率的染色体进行交配
void Crosser()
{
	int  i, j, k, kk;
	int  t; // 参与交配的染色体的个数
	int  point1, point2, temp; // 交配断点
	int  pointnum;
	int  temp1, temp2;
	int  map1[cities], map2[cities];
	double  crosser_p[num]; // 染色体的交配概率
	int  crosser_flag[num]; // 染色体的可交配情况
	int  kkk, flag = 0;

	// 初始化
	for (i = 0; i < num; i++)
	{
		crosser_flag[i] = 0;
	}

	// 随机产生交配概率
	srand((unsigned)time(NULL));
	for (i = 0; i < num; i++)
	{
		crosser_p[i] = rand() % 100;
		crosser_p[i] /= 100;
	}

	// 确定可以交配的染色体
	t = 0;
	for (i = 0; i < num; i++)
	{
		if (crosser_p[i] < pc)
		{
			crosser_flag[i] = 1;
			t++;
		}
	}

	t = t / 2 * 2;
	// t必须为偶数，产生t/2个0-9交配断点
	srand((unsigned)time(NULL));
	temp1 = 0;
	// temp1号染色体和temp2号染色体交配
	for (i = 0; i < t / 2; i++)
		// 如果有5个染色体需要交配，但是实际上t/2代表只有4个染色体会真正地交配，剩下的1个才加上5个不需要交配的染色体直接进入下一代。
	{
		point1 = rand() % cities;
		// 交配点1
		point2 = rand() % cities;
		// 交配点2

		// 选出一个需要交配的染色体1
		for (j = temp1; j < num; j++)
		{
			if (crosser_flag[j] == 1)
			{
				temp1 = j;
				break;
			}
		}

		// 选出另一个需要交配的染色体2与1交配
		for (j = temp1 + 1; j < num; j++)
		{
			if (crosser_flag[j] == 1)
			{
				temp2 = j;
				break;
			}
		}

		// 进行基因交配
		if (point1 > point2)
			// 保证point1 <= point2
		{
			temp = point1;
			point1 = point2;
			point2 = point1;
		}

		// 初始化
		memset(map1, -1, sizeof(map1));
		memset(map2, -1, sizeof(map2));

		// 断点之间的基因产生映射
		for (k = point1; k <= point2; k++)
		{
			map1[group[temp1].city[k]] = group[temp2].city[k];
			map2[group[temp2].city[k]] = group[temp1].city[k];
		}

		// 断点两边的基因互换
		for (k = 0; k < point1; k++)
		{
			temp = group[temp1].city[k];
			group[temp1].city[k] = group[temp2].city[k];
			group[temp2].city[k] = temp;
		}
		for (k = point2 + 1; k < cities; k++)
		{
			temp = group[temp1].city[k];
			group[temp1].city[k] = group[temp2].city[k];
			group[temp2].city[k] = temp;
		}

		// 处理染色体1产生的冲突基因
		for (k = 0; k < point1; k++)
		{
			for (kk = point1; kk <= point2; kk++)
			{
				if (group[temp1].city[k] == group[temp1].city[kk])
				{
					group[temp1].city[k] = map1[group[temp1].city[k]];
					// 如果相等则进行映射操作
					// find
					for (kkk = point1; kkk <= point2; kkk++)
					{
						if (group[temp1].city[k] == group[temp1].city[kkk])
							// 考虑如果映射一次仍然具有相同的城市，则在一次映射操作
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
						// flag不断判断统一染色体中是否还存在相同的城市
					{
						kk = point1 - 1;
						flag = 0;
					}
					else
					{
						flag = 0;
						break;
					}
				}
			}
		}
		for (k = point2 + 1; k < cities; k++)
		{
			for (kk = point1; kk <= point2; kk++)
			{
				if (group[temp1].city[k] == group[temp1].city[kk])
				{
					group[temp1].city[k] = map1[group[temp1].city[k]];
					// 如果相等则进行映射操作
					// find
					for (kkk = point1; kkk <= point2; kkk++)
					{
						if (group[temp1].city[k] == group[temp1].city[kkk])
							// 考虑如果映射一次仍然具有相同的城市，则在一次映射操作
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
						// flag不断判断统一染色体中是否还存在相同的城市
					{
						kk = point1 - 1;
						flag = 0;
					}
					else
					{
						flag = 0;
						break;
					}
				}
			}
		}

		// 处理染色体2产生的冲突基因
		for (k = 0; k < point1; k++)
		{
			for (kk = point1; kk <= point2; kk++)
			{
				if (group[temp2].city[k] == group[temp2].city[kk])
				{
					group[temp2].city[k] = map2[group[temp2].city[k]];
					// 如果相等则进行映射操作
					// find
					for (kkk = point1; kkk <= point2; kkk++)
					{
						if (group[temp2].city[k] == group[temp2].city[kkk])
							// 考虑如果映射一次仍然具有相同的城市，则在一次映射操作
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
						// flag不断判断统一染色体中是否还存在相同的城市
					{
						kk = point1 - 1;
						flag = 0;
					}
					else
					{
						flag = 0;
						break;
					}
				}
			}
		}
		for (k = point2 + 1; k < cities; k++)
		{
			for (kk = point1; kk <= point2; kk++)
			{
				if (group[temp2].city[k] == group[temp2].city[kk])
				{
					group[temp2].city[k] = map2[group[temp2].city[k]];
					// 如果相等则进行映射操作
					// find
					for (kkk = point1; kkk <= point2; kkk++)
					{
						if (group[temp2].city[k] == group[temp2].city[kkk])
							// 考虑如果映射一次仍然具有相同的城市，则在一次映射操作
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
						// flag不断判断统一染色体中是否还存在相同的城市
					{
						kk = point1 - 1;
						flag = 0;
					}
					else
					{
						flag = 0;
						break;
					}
				}
			}
		}

		temp1 = temp2 + 1;
	}
}

void Mutation()
{
	int  i, j;
	int  t;
	int  temp1, temp2, point;
	double  mutation_p[num]; // 染色体的变异概率
	int  mutation_flag[num]; // 染色体的变异情况
	// 初始化
	for (i = 0; i < num; i++)
	{
		mutation_flag[i] = 0;
	}
	// 随机产生变异概率
	srand((unsigned)time(NULL));
	for (i = 0; i < num; i++)
	{
		mutation_p[i] = rand() % 100;
		mutation_p[i] /= 100;
	}

	// 确定可以变异的染色体
	t = 0;
	for (i = 0; i < num; i++)
	{
		if (mutation_p[i] < pm)
		{
			mutation_flag[i] = 1;
			t++;
		}
	}

	// 变异操作，即交换染色体的两个节点
	srand((unsigned)time(NULL));
	for (i = 0; i < num; i++)
	{
		if (mutation_flag[i] == 1)
		{
			temp1 = rand() % cities;
			temp2 = rand() % cities;
			point = group[i].city[temp1];
			group[i].city[temp1] = group[i].city[temp2];
			group[i].city[temp2] = point;
		}
	}
}

int main()
{
	int  i, j, t;
	init();
	groupproduce();
	// 初始种群评价
	Fitness();
	t = 0;
	while (t++ < MAXX)
	{
		Select();
		Crosser();
		Mutation();
		Fitness();
	}

	// 最终种群的评价
	cout << endl << "输出最终的种群评价: " << endl;
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < cities; j++)
		{
			cout << group[i].city[j] << " ";
		}
		cout << endl << "adapt: " << group[i].adapt << "\tp: " << group[i].p << endl;
	}

	cout << "最优解为" << bestsolution << "条染色体" << endl;
	system("pause");
	return 0;
}
