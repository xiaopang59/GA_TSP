#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>
#include <time.h>

using namespace std;

const int cities	= 10; // ���и���
const int MAXX		= 100; // ��������
const double pc		= 0.8; // �������
const double pm		= 0.05; // �������
const int num		= 50; // ��Ⱥ�Ĵ�С

int bestsolution;
// ����Ⱦɫ��
double distance1[cities][cities];
// ����֮��ľ���

struct _group
{
	int  city[cities]; // ���е�˳��
	double  adapt; // ��Ӧ��
	double  p; // ����Ⱥ�е��Ҵ����
}group[num], grouptemp[num], groupsum[3 * num];


struct point
{
	double  x, y;
}Cpoint[cities];

// �������cities������֮����໥����
void init()
{
	//int i, j;
	//memset(distance1, 0, sizeof(distance1));
	//srand((unsigned)time(NULL));
	//for (i = 0; i < cities; i++)
	//{
	//  for (j = i + 1; j < cities; j++)
	//  {
	//  distance1[i][j] = rand() % 100;
	//  distance1[j][i] = distance1[i][j];
	//  }
	//}

	ifstream f1("city.txt", ios::in);
	for (int i = 0; i < cities; i++)
	{
		//if(ll % 2 != 0)
		f1 >> Cpoint[i].x;
		//if(ll % 2 == 0)
		f1 >> Cpoint[i].y;
		//ll++;
	}
	f1.close();
	double tmp1 = 0.0;
	for (int m = 0; m < cities; m++)
	{
		for (int n = 0; n < cities; n++)
		{
			if (m == n)
			{
				distance1[m][n] = 0;
			}
			else
			{
				tmp1 = sqrt(pow(Cpoint[m].x - Cpoint[n].x, 2) + pow(Cpoint[m].y - Cpoint[n].y, 2));
				distance1[m][n] = tmp1;
				distance1[n][m] = tmp1;
			}
		}
	}

	// ��ӡ�������
	cout << "���еľ����������: " << endl;
	for (int i = 0; i < cities; i++)
	{
		for (int j = 0; j < cities; j++)
		{
			cout << distance1[i][j] << "\t";
		}
		cout << endl;
	}
}

// ���������ʼ��Ⱥ
void groupproduce()
{
	int i, j, t, k, flag;

	// ��ʼ��
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
		// ����10������ͬ������
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

	// ��ӡ��Ⱥ����
	cout << "��ʼ����Ⱥ" << endl;
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < cities; j++)
		{
			cout << group[i].city[j] << " ";
		}
		cout << endl;
	}
}

// ���ۺ������ҳ�����Ⱦɫ��
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
		// ÿ��Ⱦɫ���·���ܺ�
		biggestsum += sumdistance;
		// ��Ⱥ����·��
	}

	// ����Ⱦɫ����Ҵ�������·��Խ���������Խ��
	for (i = 0; i < num; i++)
	{
		group[i].p = 1 - (double)group[i].adapt / (double)biggestsum;
		biggestp += group[i].p;
	}

	for (i = 0; i < num; i++)
	{
		group[i].p /= biggestp;
	}// ����Ⱥ�е��Ҵ���ʣ��ܺ�Ϊ1
	// �����·��
	bestsolution = 0;
	for (i = 0; i < num; i++)
	{
		if (group[i].p > group[bestsolution].p)
		{
			bestsolution = i;
		}
	}
}

// ������Ӧ��
void  CalFitness(_group *Group)
{
	int  j;
	int  n1, n2;
	double  sumdistance = 0;
	for (j = 1; j < cities; j++)
	{
		n1 = Group->city[j - 1];
		n2 = Group->city[j];
		sumdistance += distance1[n1][n2];
	}
	sumdistance += distance1[Group->city[cities - 1]][Group->city[0]];
	Group->adapt = sumdistance;
	// ÿ��Ⱦɫ���·���ܺ�
}

// ѡ��
void Select()
{
	int  i, j, temp;
	double  gradient_p[num]; // �ݶȸ���
	double  select_p[num]; // ѡ��Ⱦɫ����������
	int  select[num];

	// ��ʼ���ݶȸ���
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

	// �������Ⱦɫ��Ĵ�����
	for (i = 0; i < num; i++)
	{
		select_p[i] = rand() % 100;
		select_p[i] /= 100;
	}

	// ѡ���������Ⱦɫ��
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < num; j++)
		{
			if (select_p[i] < gradient_p[j])
			{
				select[i] = j;
				// ��i��λ�ô�ŵ�j��Ⱦɫ��
				break;
			}
		}
	}

	// ������Ⱥ
	for (i = 0; i < num; i++)
	{
		groupsum[i].adapt
			= group[i].adapt;
		groupsum[i].p
			= group[i].p;
		for (j = 0; j < cities; j++)
		{
			groupsum[i].city[j] = group[i].city[j];
		}
	}

	// ���ݸ���
	for (i = 0; i < num; i++)
	{
		temp = select[i];
		grouptemp[i].adapt
			= group[temp].adapt;
		grouptemp[i].p
			= group[temp].p;
		for (j = 0; j < cities; j++)
		{
			grouptemp[i].city[j] = group[temp].city[j];
		}
	}
}

// ���䣬��ÿ��Ⱦɫ�����������ʣ����㽻���ʵ�Ⱦɫ����н���
int Crosser()
{
	int  i, j, k, kk;
	int  t; // ���뽻���Ⱦɫ��ĸ���
	int  point1, point2, temp; // ����ϵ�
	int  pointnum;
	int  temp1, temp2;
	int  map1[cities], map2[cities];
	double  crosser_p[num]; // Ⱦɫ��Ľ������
	int  crosser_flag[num]; // Ⱦɫ��Ŀɽ������
	int  kkk, flag = 0;

	// ��ʼ��
	for (i = 0; i < num; i++)
	{
		crosser_flag[i] = 0;
	}

	// ��������������
	srand((unsigned)time(NULL));
	for (i = 0; i < num; i++)
	{
		crosser_p[i] = rand() % 100;
		crosser_p[i] /= 100;
	}

	// ȷ�����Խ����Ⱦɫ��
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
	// t����Ϊż��������t/2��0-9����ϵ�
	srand((unsigned)time(NULL));
	temp1 = 0;
	// temp1��Ⱦɫ���temp2��Ⱦɫ�彻��
	for (i = 0; i < t / 2; i++)
		// �����5��Ⱦɫ����Ҫ���䣬����ʵ����t/2����ֻ��4��Ⱦɫ��������ؽ��䣬ʣ�µ�1���ż���5������Ҫ�����Ⱦɫ��ֱ�ӽ�����һ����
	{
		point1 = rand() % cities;
		// �����1
		point2 = rand() % cities;
		// �����2

		// ѡ��һ����Ҫ�����Ⱦɫ��1
		for (j = temp1; j < num; j++)
		{
			if (crosser_flag[j] == 1)
			{
				temp1 = j;
				break;
			}
		}

		// ѡ����һ����Ҫ�����Ⱦɫ��2��1����
		for (j = temp1 + 1; j < num; j++)
		{
			if (crosser_flag[j] == 1)
			{
				temp2 = j;
				break;
			}
		}

		// ����һ������Ⱦɫ����и���
		groupsum[num + 2 * i].adapt
			= grouptemp[temp1].adapt;
		groupsum[num + 2 * i].p
			= grouptemp[temp1].p;
		for (j = 0; j < cities; j++)
		{
			groupsum[num + 2 * i].city[j] = grouptemp[temp1].city[j];
		}

		// ���ڶ�������Ⱦɫ����и���
		groupsum[num + 2 * i + 1].adapt
			= grouptemp[temp2].adapt;
		groupsum[num + 2 * i + 1].p
			= grouptemp[temp2].p;
		for (j = 0; j < cities; j++)
		{
			groupsum[num + 2 * i + 1].city[j] = grouptemp[temp2].city[j];
		}

		// ���л�����
		if (point1 > point2)
			// ��֤point1 <= point2
		{
			temp = point1;
			point1 = point2;
			point2 = point1;
		}

		// ��ʼ��
		memset(map1, -1, sizeof(map1));
		memset(map2, -1, sizeof(map2));

		// �ϵ�֮��Ļ������ӳ��
		for (k = point1; k <= point2; k++)
		{
			map1[grouptemp[temp1].city[k]] = grouptemp[temp2].city[k];
			map2[grouptemp[temp2].city[k]] = grouptemp[temp1].city[k];
		}

		// �ϵ����ߵĻ��򻥻�
		for (k = 0; k < point1; k++)
		{
			temp = groupsum[num + 2 * i].city[k];
			groupsum[num + 2 * i].city[k] = groupsum[num + 2 * i + 1].city[k];
			groupsum[num + 2 * i + 1].city[k] = temp;
		}
		for (k = point2 + 1; k < cities; k++)
		{
			temp = groupsum[num + 2 * i].city[k];
			groupsum[num + 2 * i].city[k] = groupsum[num + 2 * i + 1].city[k];
			groupsum[num + 2 * i + 1].city[k] = temp;
		}

		// ����Ⱦɫ��1�����ĳ�ͻ����
		for (k = 0; k < point1; k++)
		{
			for (kk = point1; kk <= point2; kk++)
			{
				if (groupsum[num + 2 * i].city[k] == groupsum[num + 2 * i].city[kk])
				{
					groupsum[num + 2 * i].city[k] = map1[groupsum[num + 2 * i].city[k]];
					// �����������ӳ�����
					// find
					for (kkk = point1; kkk <= point2; kkk++)
					{
						if (groupsum[num + 2 * i].city[k] == groupsum[num + 2 * i].city[kkk])
							// �������ӳ��һ����Ȼ������ͬ�ĳ��У�����һ��ӳ�����
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
						// flag�����ж�ͳһȾɫ�����Ƿ񻹴�����ͬ�ĳ���
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
				if (groupsum[num + 2 * i].city[k] == groupsum[num + 2 * i].city[kk])
				{
					groupsum[num + 2 * i].city[k] = map1[groupsum[num + 2 * i].city[k]];
					// �����������ӳ�����
					// find
					for (kkk = point1; kkk <= point2; kkk++)
					{
						if (groupsum[num + 2 * i].city[k] == groupsum[num + 2 * i].city[kkk])
							// �������ӳ��һ����Ȼ������ͬ�ĳ��У�����һ��ӳ�����
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
						// flag�����ж�ͳһȾɫ�����Ƿ񻹴�����ͬ�ĳ���
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

		// ����Ⱦɫ��2�����ĳ�ͻ����
		for (k = 0; k < point1; k++)
		{
			for (kk = point1; kk <= point2; kk++)
			{
				if (groupsum[num + 2 * i + 1].city[k] == groupsum[num + 2 * i + 1].city[kk])
				{
					groupsum[num + 2 * i + 1].city[k] = map2[groupsum[num + 2 * i + 1].city[k]];
					// �����������ӳ�����
					// find
					for (kkk = point1; kkk <= point2; kkk++)
					{
						if (groupsum[num + 2 * i + 1].city[k] == groupsum[num + 2 * i + 1].city[kkk])
							// �������ӳ��һ����Ȼ������ͬ�ĳ��У�����һ��ӳ�����
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
						// flag�����ж�ͳһȾɫ�����Ƿ񻹴�����ͬ�ĳ���
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
				if (groupsum[num + 2 * i + 1].city[k] == groupsum[num + 2 * i + 1].city[kk])
				{
					groupsum[num + 2 * i + 1].city[k] = map2[groupsum[num + 2 * i + 1].city[k]];
					// �����������ӳ�����
					// find
					for (kkk = point1; kkk <= point2; kkk++)
					{
						if (groupsum[num + 2 * i + 1].city[k] == groupsum[num + 2 * i + 1].city[kkk])
							// �������ӳ��һ����Ȼ������ͬ�ĳ��У�����һ��ӳ�����
						{
							flag = 1;
							break;
						}
					}
					if (flag == 1)
						// flag�����ж�ͳһȾɫ�����Ƿ񻹴�����ͬ�ĳ���
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

		// ������Ӧ���Ӵ�Ⱦɫ�����Ӧ��
		CalFitness(&groupsum[num + 2 * i]);
		CalFitness(&groupsum[num + 2 * i + 1]);
	}

	return 2 * i;
}

int Mutation(int CrosserNum)
{
	int  i, j;
	int  t;
	int  temp1, temp2, point;
	double  mutation_p[num]; // Ⱦɫ��ı������
	int  mutation_flag[num]; // Ⱦɫ��ı������
	int  count = 0;
	// ��ʼ��
	for (i = 0; i < num; i++)
	{
		mutation_flag[i] = 0;
	}
	// ��������������
	srand((unsigned)time(NULL));
	for (i = 0; i < num; i++)
	{
		mutation_p[i] = rand() % 100;
		mutation_p[i] /= 100;
	}

	// ȷ�����Ա����Ⱦɫ��
	t = 0;
	for (i = 0; i < num; i++)
	{
		if (mutation_p[i] < pm)
		{
			mutation_flag[i] = 1;
			t++;
		}
	}

	// ���������������Ⱦɫ��������ڵ�
	srand((unsigned)time(NULL));
	for (i = 0; i < num; i++)
	{
		if (mutation_flag[i] == 1)
		{
			// ������Ⱦɫ����и���
			groupsum[num + CrosserNum + count].adapt
				= group[i].adapt;
			groupsum[num + CrosserNum + count].p
				= group[i].p;
			for (j = 0; j < cities; j++)
			{
				groupsum[num + CrosserNum + count].city[j] = group[i].city[j];
			}

			temp1 = rand() % cities;
			temp2 = rand() % cities;
			point = group[i].city[temp1];
			groupsum[num + CrosserNum + count].city[temp1] = groupsum[num + CrosserNum + count].city[temp2];
			groupsum[num + CrosserNum + count].city[temp2] = point;

			count++;
		}
	}

	return t;
}

// ������Ⱥ������滻
void GroupReplace()
{
	int i, j;

	for (int i = 0; i < num; i++)
	{
		group[i].adapt
			= groupsum[i].adapt;
		group[i].p
			= groupsum[i].p;
		for (j = 0; j < cities; j++)
		{
			group[i].city[j] = groupsum[i].city[j];
		}
	}
}

// ���㸸����ѡ��ĸ���
void CalProSelect()
{
	int  i, j;
	int  biggestsum = 0;
	double  biggestp = 0;

	for (i = 0; i < num; i++)
	{
		biggestsum += group[i].adapt;
		// ��Ⱥ����·��
	}

	// ����Ⱦɫ����Ҵ�������·��Խ���������Խ��
	for (i = 0; i < num; i++)
	{
		group[i].p = 1 - (double)group[i].adapt / (double)biggestsum;
		biggestp += group[i].p;
	}

	for (i = 0; i < num; i++)
	{
		group[i].p /= biggestp;
	}// ����Ⱥ�е��Ҵ���ʣ��ܺ�Ϊ1
	// �����·��
	bestsolution = 0;
	for (i = 0; i < num; i++)
	{
		if (group[i].p > group[bestsolution].p)
		{
			bestsolution = i;
		}
	}
}

bool GroupCompare(_group &Group1, _group &Group2)
{
	return Group1.adapt < Group2.adapt;
}

int main()
{
	int  i, j, t;
	int  CrosserNum, MutationNum;
	init();
	groupproduce();
	// ��ʼ��Ⱥ����
	Fitness();
	t = 0;
	while (t++ < MAXX)
	{
		Select();
		CrosserNum = Crosser();
		MutationNum = Mutation(CrosserNum);
		sort(groupsum, groupsum + num + CrosserNum + MutationNum, GroupCompare);
		// ����Ⱥ���塢����Ⱥ��ͱ���Ⱥ���������
		GroupReplace();
		// �滻����Ⱥ��
		CalProSelect();
		// ���㸸����ѡ��ĸ���
	}

	// ������Ⱥ������
	cout << endl << "������յ���Ⱥ����: " << endl;
	for (i = 0; i < num; i++)
	{
		for (j = 0; j < cities; j++)
		{
			cout << group[i].city[j] << " ";
		}
		cout << endl << "adapt: " << group[i].adapt << "\tp: " << group[i].p << endl;
	}

	cout << "���Ž�Ϊ" << bestsolution << "��Ⱦɫ��" << endl;

	system("pause");
	return 0;
}
