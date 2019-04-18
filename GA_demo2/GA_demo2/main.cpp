#include <iostream>
#include "Genetic.h"

using namespace std;

extern Individual Population[GROUP_SCALE + 1];

void main()
{
	int  Xnration;
	int  i;
	int  seed = 123456789;

	showTime();
	initGroup(seed);
	evaluate();
	selectBest();

	for (Xnration = 0; Xnration < MAX_GENS; Xnration++)
	{
		selector(seed);
		crossover(seed);
		mutate(seed);
		report(seed);
		evaluate();
		elitist();
	}

	cout																<<	endl;
	cout << "  Best member after "	<<	MAX_GENS	<<	"  Xnrations: "	<<	endl;
	cout																<<	endl;


	for (i = 0; i < N_VARS; i++)
	{
		cout << "  X(" << i + 1 << ") = " << Population[GROUP_SCALE].Xn[i] <<	endl;
	}
	cout															<< endl;
	cout << " Best Fitness = "	<<	Population[GROUP_SCALE].Fitness	<<	endl;

	showTime();
	system("pause");
}