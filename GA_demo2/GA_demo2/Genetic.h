#pragma once

using namespace std;

#define PI 3.1415926535897923846


// �Ŵ��㷨��������Ⱥ��ģ( 0~100 )����ֳ��������������������������ʡ��������
#define GROUP_SCALE	50
#define MAX_GENS	500
#define N_VARS		2
#define P_MATING	0.8
#define P_MUTATION	0.15

struct Individual
{
	double  Xn[N_VARS]; // ��ű���ֵ
	double  Fitness; // ��Ӧֵ
	double  ReFitness; // ��Ӧֵ�����ܶ�
	double  SumFitness; // �ۼӷֲ���Ϊ����ת
};

struct X_Range
{
	double  Upper; // �������Ͻ�ȡֵ
	double  Lower; // �������½�ȡֵ
};


template<typename T>
T randT(T Lower, T Upper); // ���������������������


void crossover(int &seed);
void elitist(); // ������
void evaluate();


void initGroup(int &seed);


void selectBest();
void mutate(int &seed);


double r8_uniform_ab(double a, double b, int &seed);
int i4_uniform_ab(int a, int b, int c);


void report(int Xnration);
void selector(int &seed);
void showTime();
void Xover(int one, int two, int &seed);
