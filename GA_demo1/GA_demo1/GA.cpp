#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "GA.h"

using namespace std;



double random()
{
	double randNum = rand() * 1.0 / RAND_MAX;
	return randNum;
}



GenAlg::GenAlg()
{

}



void GenAlg::init(const int
	popsize,
	const double
	MutRate,
	const double
	CrossRate,
	const int
	GenLength,
	const double
	LeftPoint,
	const double
	RightPoint)
{
	popSize = popsize;
	mutationRate
		= MutRate;
	crossoverRate
		= CrossRate;
	chromoLength
		= GenLength;
	totalFitness
		= 0;
	generation
		= 0;
	bestFitness
		= 0.0;
	worstFitness
		= 999999;
	averageFitness
		= 0;
	maxPerturbation
		= 0.004;
	leftPoint
		= LeftPoint;
	rightPoint
		= RightPoint;

	// 清空种群容器，以初始化
	vecPop.clear();

	for (int i = 0; i < popSize; i++)
	{
		// 类的构造函数已经把适应性评分初始化为0
		vecPop.push_back(Genome());

		// 把所有的基因编码初始化为函数区间内的随机数。
		for (int j = 0; j <chromoLength; j++)
		{
			vecPop[i].vecGenome.push_back(random() * (rightPoint - leftPoint) + leftPoint);
		}
	}
}



void GenAlg::Reset()
{
	totalFitness
		= 0;
	averageFitness
		= 0;
}



void GenAlg::CalculateBestWorstAvTot()
{
	for (int i = 0; i < popSize; i++)
	{
		// 累计适应性分数
		totalFitness += vecPop[i].fitness;

		if (vecPop[i].fitness >= bestFitness)
		{
			bestFitness
				= vecPop[i].fitness;
			fittestGenome
				= vecPop[i];
		}

		if (vecPop[i].fitness <= worstFitness)
		{
			worstFitness = vecPop[i].fitness;
		}
	}

	averageFitness = totalFitness / popSize;
}



Genome GenAlg::GetChromoRoulette()
{
	// 产生一个0到人口总适应性评分总和之间的随机数
	// 中 m_dTotalFitness 记录了整个种群的适应性分数总和)
	double Slice = (random()) * totalFitness;

	// 这个基因将承载转盘所选出来的那个个体.
	Genome TheChosenOne;

	// 累计适应性分数的和.
	double FitnessSoFar = 0;

	// 遍历总人口里面的每一条染色体
	for (int i = 0; i < popSize; i++)
	{
		// 累计适应性分数
		FitnessSoFar += vecPop[i].fitness;

		// 如果累计分数大于随机数，就选择此时的基因.
		if (FitnessSoFar >= Slice)
		{
			TheChosenOne = vecPop[i];
			break;
		}
	}

	// 返回轮盘选出来的个体基因
	return TheChosenOne;
}



void GenAlg::Mutate(vector<double> &chromo)
{
	// 遵循预定的突变概率，对基因进行突变
	for (int i = 0; i < chromo.size(); i++)
	{
		// 如果发生突变的话
		if (random() < mutationRate)
		{
			// 使该权值增加或减少一个很小的随机数值
			chromo[i] += ((random() - 0.5) * maxPerturbation);

			// 限定范围
			if (chromo[i] < leftPoint)
			{
				chromo[i] = leftPoint;
			}
			else if (chromo[i] > rightPoint)
			{
				chromo[i] = rightPoint;
			}

			// 以上代码非基因变异的一般性代码只是用来保证基因编码的可行性。
		}
	}
}



// 此函数产生新的一代，见证着整个进化的全过程.
// 以父代种群的基因组容器作为参数传进去，该函数将往该容器里放入新一代的基因组(当然是经过了优胜劣汰的)
void GenAlg::Epoch(vector<Genome> &vecNewPop)
{
	// 用类的成员变量来存储父代的基因组(在此之前 m_vecPop 存储的是不带估值的所有基因组)
	vecPop = vecNewPop;

	// 初始化相关变量
	Reset();

	// 为相关变量赋值
	CalculateBestWorstAvTot();

	// 清空装载新种群的容器
	vecNewPop.clear();

	// 产生新一代的所有基因组
	while (vecNewPop.size() < popSize)
	{
		// 转盘随机抽出两个基因
		Genome mum
			= GetChromoRoulette();
		Genome dad
			= GetChromoRoulette();

		// 创建两个子代基因组
		vector<double> baby1, baby2;

		// 先把它们分别设置为父方和母方的基因
		baby1 = mum.vecGenome;
		baby2 = dad.vecGenome;
		// 使子代基因发生基因突变
		Mutate(baby1);
		Mutate(baby2);

		// 把两个子代基因放到新的基因组容器里面
		vecNewPop.push_back(Genome(baby1, 0));
		vecNewPop.push_back(Genome(baby2, 0));

	}  // 子代产生完毕

	// 如果你设置的人口总数非单数的话，就会出现报错
	if (vecNewPop.size() != popSize)
	{
		cout << "error" << endl;
		return;
	}
}



Genome GenAlg::GetBestFitness()
{
	return fittestGenome;
}



double GenAlg::GetAverageFitness()
{
	return averageFitness;
}



void GenEngine::report(const int &genNum)
{
	cout << "第"
		<< genNum << "代"
		<< endl;
	cout << "最佳适应度: "
		<< bestFitness << endl;
	cout << "最佳适应度基因取值: "
		<< bestSearch << endl;
	cout << "平均适应度: "
		<< averageFitness << endl
		<< endl;
}



void GenEngine::OnStartGenAlg()
{
	// 产生随机数
	srand((unsigned)time(NULL));

	// 初始化遗传算法引擎
	genAlg.init(g_popsize, g_dMutationRate, g_dCrossoverRate, g_numGen, g_LeftPoint, g_RightPoint);

	// 清空种群容器
	m_population.clear();

	// 种群容器装进经过随机初始化的种群
	m_population = genAlg.vecPop;

	vector<double> input;
	double output;

	input.push_back(0);

	for (int Generation = 0; Generation <= g_Generation; Generation++)
	{
		// 里面是对每一条染色体进行操作
		for (int i = 0; i < g_popsize; i++)
		{
			input = m_population[i].vecGenome;

			// 为每一个个体做适应性评价，如之前说的，评价分数就是函数值。其
			// Function函数的作用是输入自身量返回函数值，读者可以参考其代码
			output = (double)curve.function(input);
			m_population[i].fitness
				= output;
		}

		// 由父代种群进化出子代种群
		genAlg.Epoch(m_population);

		bestSearch
			= genAlg.GetBestFitness().vecGenome[0];
		bestFitness
			= genAlg.GetBestFitness().fitness;
		averageFitness
			= genAlg.GetAverageFitness();
		report(Generation + 1);
	}
}