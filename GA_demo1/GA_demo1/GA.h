#pragma once

#include <vector>
#include <math.h>

using namespace std;


const double pai = 3.1415926;



class Genome
{
public:
	friend class GenAlg;
	friend class GenEngine;

	Genome() : fitness(0){}

	Genome(vector<double>vec, double f) : vecGenome(vec), fitness(f){}
	// ��Ĵ�������ʼ���Ĺ��캯��

private:
	vector<double>
		vecGenome; // dFitness���ڴ洢�Ըû������Ӧ������
	double  fitness; // ����޲�����ʼ������
};



// �Ŵ��㷨
class GenAlg
{
public:
	friend class GenEngine;

	// ���캯��
	GenAlg();

	// ��ʼ������ 
	void Reset();

	// ��ʼ������
	void init(
		const int popsize,
		const double
		MutRate,
		const double
		CrossRate,
		const int
		GenLength,
		const double
		LeftPoint,
		const double
		RightPoint);

	// ����TotalFitness, BestFitness, WorstFitness, AverageFitness�ȱ���
	void CalculateBestWorstAvTot();

	// ���̶�ѡ����
	Genome GetChromoRoulette();

	// ������캯��
	void Mutate(vector<double> &chromo);

	// �ú���������һ������
	void Epoch(vector<Genome> &vecNewPop);
	Genome GetBestFitness();
	double GetAverageFitness();

private:
	vector<Genome>
		vecPop; // ����������洢ÿһ�������Ⱦɫ��
	int  popSize; // �˿�(��Ⱥ)����
	int  chromoLength; // ÿһ��Ⱦɫ��Ļ��������Ŀ
	double  totalFitness; // ���и����Ӧ����Ӧ�����ֵ��ܺ�
	double  bestFitness; // �����и��嵱������Ӧ�ĸ������Ӧ������
	double  averageFitness; // ���и������Ӧ�����ֵ�ƽ��ֵ
	double  worstFitness; // �����и��嵱�����Ӧ�ĸ������Ӧ������
	Genome  fittestGenome; // ����Ӧ�ĸ�����m_vecPop���������������
	double  mutationRate; // ����ͻ��ĸ��ʣ�һ�����0.05��0.3֮��
	double  crossoverRate; // ���򽻲�ĸ���һ����Ϊ0.7
	int  generation; // �����ļ�����
	double  maxPerturbation; // �����첽��
	double  leftPoint;
	double  rightPoint;
};



class Curve
{
public:
	double function(const vector<double> &input)
	{
		double x = input[0];
		double output
			= x * sin(10 * pai * x) + 2.0;
		return output;
	}

private:

};



// �Ŵ���������
class GenEngine
{
public:
	GenEngine(
		const int &popsize,
		const double
		&mutationRate,
		const double
		&crossoverRate,
		const int
		&numGen,
		const int
		&generation,
		const double
		&leftPoint,
		const double
		&rightPoint) : genAlg(), curve(), m_population()
	{
		g_popsize
			= popsize;
		g_dMutationRate
			= mutationRate;
		g_dCrossoverRate
			= crossoverRate;
		g_numGen = numGen;
		g_Generation
			= generation;
		g_LeftPoint
			= leftPoint;
		g_RightPoint
			= rightPoint;
		bestFitness
			= 0;
		bestSearch
			= 0;
	}

	void OnStartGenAlg();

	// ����ÿһ�����������
	void report(const int &genNum);

private:
	GenAlg  genAlg;
	Curve  curve;
	vector<Genome>
		m_population;
	int  g_popsize;
	double  g_dMutationRate;
	double  g_dCrossoverRate;
	int  g_numGen;
	int  g_Generation;
	double  g_LeftPoint;
	double  g_RightPoint;
	double  bestFitness;
	double  bestSearch;
	double  averageFitness;
};


