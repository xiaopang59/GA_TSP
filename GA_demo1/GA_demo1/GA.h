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
	// 类的带参数初始化的构造函数

private:
	vector<double>
		vecGenome; // dFitness用于存储对该基因的适应性评估
	double  fitness; // 类的无参数初始化参数
};



// 遗传算法
class GenAlg
{
public:
	friend class GenEngine;

	// 构造函数
	GenAlg();

	// 初始化变量 
	void Reset();

	// 初始化函数
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

	// 计算TotalFitness, BestFitness, WorstFitness, AverageFitness等变量
	void CalculateBestWorstAvTot();

	// 轮盘赌选择函数
	Genome GetChromoRoulette();

	// 基因变异函数
	void Mutate(vector<double> &chromo);

	// 该函数产生新一代基因
	void Epoch(vector<Genome> &vecNewPop);
	Genome GetBestFitness();
	double GetAverageFitness();

private:
	vector<Genome>
		vecPop; // 这个容器将存储每一个个体的染色体
	int  popSize; // 人口(种群)数量
	int  chromoLength; // 每一条染色体的基因的总数目
	double  totalFitness; // 所有个体对应的适应性评分的总和
	double  bestFitness; // 在所有个体当中最适应的个体的适应性评分
	double  averageFitness; // 所有个体的适应性评分的平均值
	double  worstFitness; // 在所有个体当中最不适应的个体的适应性评分
	Genome  fittestGenome; // 最适应的个体在m_vecPop容器里面的索引号
	double  mutationRate; // 基因突变的概率，一般介于0.05和0.3之间
	double  crossoverRate; // 基因交叉的概率一般设为0.7
	int  generation; // 代数的计数器
	double  maxPerturbation; // 最大变异步长
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



// 遗传运算引擎
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

	// 报告每一代的运行情况
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


