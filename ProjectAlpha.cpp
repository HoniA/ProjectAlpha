//ProjectAlpha.cpp : Defines the entry point for the console application.
//

/////// //////////////////////////////////////////////////
//// ME 493
///////Project Alpha
//////////Honi Ahmadian
////////////////// Worked with Sierra Gonzales, Max Pullman, Kevin Gang

#include "stdafx.h"
#include <iostream>
#include <random>
#include <assert.h>
#include <time.h>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <fstream>

using namespace std;

#define MeanRand (double)rand()/RAND_MAX*100.0 // mean falls between 0 and 100
#define StdDevRand (double)rand()/RAND_MAX // std dev falls between 0 and 1
#define LYRAND (double)rand()/RAND_MAX // between 0 and 1. Same as StdDevRand but created again to show application

class arm {
public:
	double mean;
	double stddev;
	double value;

	void init();
	void setMean(); //to give randomly assigned mean for normal distribution
	void setStddev(); // to give randomly assigned std dev for noraml distribution
	double pull(); //Manually pulls arm, returns reward
	void updateValue(double newV);
};

void arm::init() {
	mean = -1;
	stddev = -1;
	value = 0;
}

void arm::setMean() {
	//random number
	double m = MeanRand;
	//set to mean of arm
	mean = m;
}

void arm::setStddev() {
	//random number
	double sd = StdDevRand;
	//set to std dev of arm
	stddev = sd;
}

double generateGaussianNoise(double mu, double sigma);

double arm::pull() {
	//initialize variables
	double reward;

	//creates normal distribution from mean and std dev values
	reward = generateGaussianNoise(mean, stddev);

	//returns value using normal distribution
	return reward;
}

void arm::updateValue(double newV){
	double alpha = 0.1;
	value = newV*alpha + value*(1-alpha);
}

vector<arm> createMAB(int n);

void qLearner(vector<arm> V, int numArms, int numPulls, int runs);

void TestA();

void TestB();

int main()
{
	srand(time(NULL));
	
	//initialize variables;
	int n; //number of arms
	int N; // number of pulls
	int numRuns = 30; //number of statistical runs

	//ask user for number of arms
	cout << "Give number of arms for MAB (positive integers only): ";
	cin >> n;

	vector<arm> MAB = createMAB(n);

	cout << endl << "Give number of pulls for MAB (positive integers only): ";
	cin >> N; //Number of pulls, arbitrary number 

	// run q-learner
	qLearner(MAB, n, N, numRuns);

	//run Test A (includes data for learning curve)
	TestA();

	//run Test B (includes data for action curve)
	TestB();

	system("pause");
	return 0;
}

double generateGaussianNoise(double mu, double sigma) //found on wikipedia: Box-Mueller Transformation
{
	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

vector<arm> createMAB(int numArms) {

	//create vector
	vector<arm> MAB;

	//new arm for every spot in vector
	//randomly assigned mean and std dev for each
	for (int i = 0; i < numArms; i++)
	{
		arm single_arm;
		single_arm.init();
		single_arm.setMean();
		single_arm.setStddev();

		cout << "Arm " << i + 1 << "\t" << "Mean: " << single_arm.mean << endl;

		MAB.push_back(single_arm);
	}

	return MAB;
};

void qLearner(vector<arm> V, int numArms, int numPulls, int numRuns)
{
	double epsilon = 0.1; //greedy value
	double reward;
	int decisionIndex = 0;


	//want to output results to text file for documentation
	ofstream fout;
	fout.clear();
	fout.open("qLearner.txt");

	fout << "Run" << "\t" << "Pull" << "\t" << "Arm" << endl;
	// need to update value of each arm before next iteration
	// do this N times (for loop)
	for (int run = 0; run < numRuns; run++)
	{
		for (int j = 0; j < numPulls; j++)
		{
			double greedyIndex = LYRAND;

			// chooses randomly epsilon percent (random integer + if statement)
			if (greedyIndex <= epsilon) //10% of time because epsilon = 0.1
			{
				//choose randomly
				decisionIndex = rand() % numArms; //should pick int between 0 and n-1 which is every arm
				reward = V.at(decisionIndex).pull(); // pull randomly picked arm
				V.at(decisionIndex).updateValue(reward); //update value of arm from previously generated reward
			}

			// chooses best value 1 - epsilon percent (random integer + if statement)
			else // 90% of time
			{
				//choose best value
				//assume highest reward is first arm
				int highestIndex = 0;

				for (int k = 0; k < numArms; k++)
				{
					// loop through rest of MAB to see if any arms have higher reward
					if (V.at(k).value > V.at(highestIndex).value)
					{
						// if they do, assign position that arm to the highest index 
						highestIndex = k;
					}
				}
				// pull arm with highest reward and update value
				reward = V.at(highestIndex).pull();
				V.at(highestIndex).updateValue(reward);
				decisionIndex = highestIndex;
			}

			fout << run+1 <<"\t" << j+1 << "\t" << decisionIndex + 1 << "\t"; //decision index plus one to show which arm is being pulled but starting list from 1 instead of 0
			for (int i = 0; i < numArms; i++)
			{
				fout << V.at(i).value << "\t"; // output current value of every arm
			}
			fout << endl;

		}
	}

	fout.close();
}

void TestA()
{
	int N = 1000;
	int test = 0;

	ofstream fout;
	fout.clear();
	fout.open("LearningCurve.txt");
	//outputting running average to txt file to be used for learning curve plot for report

	//create arm
	arm testArm;
	testArm.init();
	testArm.setMean();
	testArm.setStddev();

	vector<double> values;

	//pull arm N times, store value in vector
	for (int i = 0; i < N; i++)
	{
		double reward = testArm.pull();
		values.push_back(reward);
		testArm.updateValue(reward);
	}

	double runningTotal = 0.0;
	double runningAverage = 0.0;

	fout << 0 << "\t" << 0 << endl;

	// find average
	for (int j = 0; j < N; j++)
	{
		runningTotal = runningTotal + values.at(j);
		runningAverage = runningTotal / (j + 1); //takes average of previous rewards after each pull
		fout << j + 1 << "\t" << runningAverage << endl;
	}

	fout.close();

	//double average = runningTotal / N;

	if (testArm.value > 0.9*testArm.mean && testArm.value < 1.1*testArm.mean)
	{
		test = 1;
	}

	cout << "Test A:" << endl;

	cout << "Final Value: " << testArm.value << endl;
	cout << "Mean: " << testArm.mean << endl;
	
	assert(test == 1);

	cout << "Test A Completed" << endl;
}


void TestB() 
{
	ofstream fout; // output to file for action curve
	fout.clear();
	fout.open("Action_Curve.txt");

	// first line is all zeros
	fout << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;

	//one arm has much more desirable mean and standard deviation numbers
	//q-learner pulls arms N times 
	//should pull best arm significantly higher percentage of time

	//create three arms, hard code meand and std. dev for one to make it clearly better.

	int n = 2;
	vector<arm> MAB;

	for (int i = 0; i < n; i++)
	{
		arm single_arm;
		single_arm.init();
		single_arm.setMean();
		single_arm.setStddev();

		MAB.push_back(single_arm);
	}

	arm best_arm;
	best_arm.init();
	best_arm.mean = 500;
	best_arm.stddev = 1;

	MAB.push_back(best_arm);

	// q-learner
	// iterate once for every pull and once every time the best_arm is pulled 
	int N = 1000;
	double pulls = 0;
	double bestPull = 0; //third arm
	double arm2Pull = 0;
	double arm1Pull = 0;
	double epsilon = 0.25;
	double reward= 0.0;
	int decisionIndex = -1;
	int test = 0;

	for (int j = 0; j < N; j++)
	{
		double greedyIndex = LYRAND;

		// chooses randomly epsilon percent (random integer + if statement)
		if (greedyIndex <= epsilon) //10% of time because epsilon = 0.1
		{
			//choose randomly
			decisionIndex = rand() % n+1; //should pick int between 0 and n+1 which is every arm
			reward = MAB.at(decisionIndex).pull();
			MAB.at(decisionIndex).updateValue(reward);
		}

		// chooses best value 1 - epsilon percent (random integer + if statement)
		else // 90% of time
		{
			//choose best value
			//assume highest reward is first arm
			decisionIndex = 0;

			for (int k = 0; k < n+1; k++)
			{
				if (MAB.at(k).value > MAB.at(decisionIndex).value)
				{
					decisionIndex = k;
				}
			}
			reward = MAB.at(decisionIndex).pull();
			MAB.at(decisionIndex).updateValue(reward);

		}

		//iterate total pulls
		pulls = pulls+1;

		//iterate number of pulls for chosen arm
		if (decisionIndex == 2)
		{
			bestPull = bestPull+1;
		}
		else if (decisionIndex == 1)
		{
			arm2Pull = arm2Pull + 1;
		}
		else
		{
			arm1Pull = arm1Pull + 1;
		}

		// send probabilities to txt file
		fout << pulls << "\t" << arm1Pull / pulls << "\t" << arm2Pull / pulls << "\t" << bestPull / pulls << endl;
	}

	fout.close();

	cout << "Percentage of time clearly best arm pulled: " << bestPull / pulls * 100 << " %" <<endl;

	// should pull best arm 1-epsilon * 100 percent of the time; test convergence (i.e. within 5% of that value)
	if (bestPull / pulls >= 0.95*((1 - epsilon) + epsilon / (n + 1)) && bestPull / pulls <= 1.05*((1 - epsilon) + epsilon / (n + 1)))
	{
		test = 1;
	}
	assert(test == 1); 

	cout << "Test B Completed" << endl;
}