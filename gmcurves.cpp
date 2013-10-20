#include "Angel.h"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

#define NUM_POINTS 10
#define RANGE_MIN -0.75
#define RANGE_MAX 0.75
#define RANGE_CHANGE 0.05
#define WIN_HEIGHT 500
#define WIN_WIDTH 500

using namespace std;

int valI, valN, leftButtonDown = 0, indexSpot, numControlPoints = 0;
double valueX, valueY;
double differenceValue, tempDifferenceValue;
double sizeOfAxis, sizeOfControlPoints, sizeOfControlLines, sizeOfBernsteinLines, sizeOfQuadLines, sizeOfCubicLines, sizeOfFourSchemeLines;
double *xValues;
double *yValues;
vec4 *controlPoints;
vec4 *bernsteinPoints;
vec4 *bernsteinLines;
vec4 *quadPoints;
vec4 *quadLines;
vec4 *cubicPoints;
vec4 *cubicLines;
vec4 *controlLines;
vec4 *arrayAxis;
vec4 *fourSchemePoints;
vec4 *fourSchemePoints2;
vec4 *fourSchemePoints3;
vec4 *fourSchemePoints4;
vec4 *fourSchemeLines;

vec4 vColor;
GLuint vColorID;

void readFile(string filename)
{
	ifstream infile;
	infile.open(filename.c_str());
	infile >> numControlPoints;
	xValues = new double[numControlPoints];
	yValues = new double[numControlPoints];
	for (int i = 0; i < numControlPoints; i++)
	{
		infile >> xValues[i] >> yValues[i];
	}
	infile.close();
}

double interpolate(double minA, double a, double maxA, double minB, double maxB)
{    
	double denominator = (maxA-minA);    
	return minB + ((a-minA) / denominator) * (maxB-minB);
}

double factorial(int f)
{
	double i = 0, fact = 1;

	if(f <= 1)
	{
		return(1);
	}
	else
	{
		for(i = 1; i <= f; i++)
		{
			fact = fact * i;
		}
		return(fact);
	}
}

double bernstein(int i, int n, double u)
{
	double answer, binomialCoefficient;

	binomialCoefficient = (factorial(n)) / (factorial(i) * factorial(n - i));
	answer = (binomialCoefficient * pow(u, i)) * (pow(1 - u, n - i));
	return answer;
}

void createBernstein()
{
	double valU = 0;

	bernsteinPoints = new vec4[NUM_POINTS];
	for (int i = 0; i < NUM_POINTS; i++)
	{
		double bernsteinValueX = 0, bernsteinValueY = 0;
		bernsteinValueX = (xValues[0] * (bernstein(0, 3, valU))) + (xValues[1] * (bernstein(1, 3, valU))) + (xValues[2] * (bernstein(2, 3, valU))) + (xValues[3] * (bernstein(3, 3, valU)));
		bernsteinValueY = (yValues[0] * (bernstein(0, 3, valU))) + (yValues[1] * (bernstein(1, 3, valU))) + (yValues[2] * (bernstein(2, 3, valU))) + (yValues[3] * (bernstein(3, 3, valU)));
		bernsteinPoints[i] = vec4(bernsteinValueX, bernsteinValueY, 0, 1);
		valU += 1.0 / (NUM_POINTS - 1);
	}
	bernsteinLines = new vec4[(2 * NUM_POINTS) - 2];
	for (int i = 0; i < (2 * NUM_POINTS) - 2; i++)
	{
		if (i < NUM_POINTS)
		{
			bernsteinLines[i] = bernsteinPoints[i];
		}
		else
		{
			bernsteinLines[i] = bernsteinPoints[(i + 1) - NUM_POINTS];
		}
	}
}

void createQuadBSpline()
{
	double valU = 0;
	double val1 = 0, val2 = 0, val3 = 0;

	quadPoints = new vec4[NUM_POINTS];
	for (int i = 0; i < NUM_POINTS; i++)
	{
		double quadValueX = 0, quadValueY = 0;
		val1 = (((1.0/2.0) * ((pow(valU, 2)) - (2 * valU) + 1)) * (xValues[0]));
		val2 = (((1.0/2.0) * ((-2 * (pow(valU, 2))) + (2 * valU) + 1)) * (xValues[1]));
		val3 = (((1.0/2.0) * (pow(valU, 2))) * (xValues[2]));
		quadValueX = val1 + val2 + val3;
		
		val1 = (((1.0/2.0) * ((pow(valU, 2)) - (2 * valU) + 1)) * (yValues[0]));
		val2 = (((1.0/2.0) * ((-2 * (pow(valU, 2))) + (2 * valU) + 1)) * (yValues[1]));
		val3 = (((1.0/2.0) * (pow(valU, 2))) * (yValues[2]));
		quadValueY = val1 + val2 + val3;

		quadPoints[i] = vec4(quadValueX, quadValueY, 0, 1);
		valU += 1.0 / (NUM_POINTS - 1);
	}
	quadLines = new vec4[(2 * NUM_POINTS) - 2];
	for (int i = 0; i < (2 * NUM_POINTS) - 2; i++)
	{
		if (i < NUM_POINTS)
		{
			quadLines[i] = quadPoints[i];
		}
		else
		{
			quadLines[i] = quadPoints[(i + 1) - NUM_POINTS];
		}
	}
}

void createCubicBSpline()
{
	double valU = 0;
	double val1 = 0, val2 = 0, val3 = 0, val4 = 0;
	int cubicCount = 0;
	cubicPoints = new vec4[numControlPoints * NUM_POINTS];
	for (int j = 0; j < numControlPoints; j++)
	{
		valU = 0;
		for (int i = 0; i < NUM_POINTS; i++)
		{
			if (j == 0)
			{
				double cubicValueX = 0, cubicValueY = 0;
				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (xValues[0]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (xValues[1]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (xValues[2]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (xValues[3]));
				cubicValueX = val1 + val2 + val3 + val4;

				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (yValues[0]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (yValues[1]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (yValues[2]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (yValues[3]));
				cubicValueY = val1 + val2 + val3 + val4;

				cubicPoints[cubicCount] = vec4(cubicValueX, cubicValueY, 0, 1);
				valU += 1.0 / (NUM_POINTS - 1);
				cubicCount++;
			}
			else if (j == (numControlPoints - 1))
			{
				double cubicValueX = 0, cubicValueY = 0;
				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (xValues[j]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (xValues[0]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (xValues[1]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (xValues[2]));
				cubicValueX = val1 + val2 + val3 + val4;

				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (yValues[j]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (yValues[0]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (yValues[1]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (yValues[2]));
				cubicValueY = val1 + val2 + val3 + val4;

				cubicPoints[cubicCount] = vec4(cubicValueX, cubicValueY, 0, 1);
				valU += 1.0 / (NUM_POINTS - 1);
				cubicCount++;
			}
			else if (j == (numControlPoints - 2))
			{
				double cubicValueX = 0, cubicValueY = 0;
				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (xValues[j]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (xValues[j + 1]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (xValues[0]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (xValues[1]));
				cubicValueX = val1 + val2 + val3 + val4;

				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (yValues[j]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (yValues[j + 1]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (yValues[0]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (yValues[1]));
				cubicValueY = val1 + val2 + val3 + val4;

				cubicPoints[cubicCount] = vec4(cubicValueX, cubicValueY, 0, 1);
				valU += 1.0 / (NUM_POINTS - 1);
				cubicCount++;
			}
			else if (j == (numControlPoints - 3))
			{
				double cubicValueX = 0, cubicValueY = 0;
				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (xValues[j]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (xValues[j + 1]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (xValues[j + 2]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (xValues[0]));
				cubicValueX = val1 + val2 + val3 + val4;

				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (yValues[j]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (yValues[j + 1]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (yValues[j + 2]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (yValues[0]));
				cubicValueY = val1 + val2 + val3 + val4;

				cubicPoints[cubicCount] = vec4(cubicValueX, cubicValueY, 0, 1);
				valU += 1.0 / (NUM_POINTS - 1);
				cubicCount++;
			}
			else
			{
				double cubicValueX = 0, cubicValueY = 0;
				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (xValues[j]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (xValues[j + 1]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (xValues[j + 2]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (xValues[j + 3]));
				cubicValueX = val1 + val2 + val3 + val4;

				val1 = (((1.0/6.0) * ((-1 * pow(valU, 3)) + (3 * pow(valU, 2)) - (3 * valU) + 1)) * (yValues[j]));
				val2 = (((1.0/6.0) * ((3 * pow(valU, 3)) - (6 * pow(valU, 2)) + 4)) * (yValues[j + 1]));
				val3 = (((1.0/6.0) * ((-3 * pow(valU, 3)) + (3 * pow(valU, 2)) + (3 * valU) + 1)) * (yValues[j + 2]));
				val4 = (((1.0/6.0) * (pow(valU, 3))) * (yValues[j + 3]));
				cubicValueY = val1 + val2 + val3 + val4;

				cubicPoints[cubicCount] = vec4(cubicValueX, cubicValueY, 0, 1);
				valU += 1.0 / (NUM_POINTS - 1);
				cubicCount++;
			}
		}
	}

	cubicLines = new vec4[(2 * (numControlPoints * NUM_POINTS))];
	for (int i = 0; i < (2 * (numControlPoints * NUM_POINTS)) - 2; i++)
	{
		if (i < (numControlPoints * NUM_POINTS))
		{
			cubicLines[i] = cubicPoints[i];
		}
		else
		{
			cubicLines[i] = cubicPoints[(i + 1) - (numControlPoints * NUM_POINTS)];
		}
	}
}

void createFourScheme()
{
	int fourEnd = numControlPoints;
	//First
	fourSchemePoints = new vec4[2 * numControlPoints];
	for (int i = 0; i < fourEnd; i++)
	{
			fourSchemePoints[2 * i] = vec4(xValues[i], yValues[i], 0, 1);
			
			double fourSchemeValueX = 0, fourSchemeValueY = 0;
			if (i == 0)
			{
				fourSchemeValueX += (xValues[fourEnd - 1] * (-1.0 / 16.0));
				fourSchemeValueY += (yValues[fourEnd - 1] * (-1.0 / 16.0));
				fourSchemeValueX += (xValues[0] * (9.0 / 16.0));
				fourSchemeValueY += (yValues[0] * (9.0 / 16.0));
				fourSchemeValueX += (xValues[1] * (9.0 / 16.0));
				fourSchemeValueY += (yValues[1] * (9.0 / 16.0));
				fourSchemeValueX += (xValues[2] * (-1.0 / 16.0));
				fourSchemeValueY += (yValues[2] * (-1.0 / 16.0));
				fourSchemePoints[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else if (i == (fourEnd - 1))
			{
				fourSchemeValueX += (xValues[i - 1] * (-1.0 / 16.0));
				fourSchemeValueY += (yValues[i - 1] * (-1.0 / 16.0));
				fourSchemeValueX += (xValues[i] * (9.0 / 16.0));
				fourSchemeValueY += (yValues[i] * (9.0 / 16.0));
				fourSchemeValueX += (xValues[0] * (9.0 / 16.0));
				fourSchemeValueY += (yValues[0] * (9.0 / 16.0));
				fourSchemeValueX += (xValues[1] * (-1.0 / 16.0));
				fourSchemeValueY += (yValues[1] * (-1.0 / 16.0));
				fourSchemePoints[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else if (i == (fourEnd - 2))
			{
				fourSchemeValueX += (xValues[i - 1] * (-1.0 / 16.0));
				fourSchemeValueY += (yValues[i - 1] * (-1.0 / 16.0));
				fourSchemeValueX += (xValues[i] * (9.0 / 16.0));
				fourSchemeValueY += (yValues[i] * (9.0 / 16.0));
				fourSchemeValueX += (xValues[i + 1] * (9.0 / 16.0));
				fourSchemeValueY += (yValues[i + 1] * (9.0 / 16.0));
				fourSchemeValueX += (xValues[0] * (-1.0 / 16.0));
				fourSchemeValueY += (yValues[0] * (-1.0 / 16.0));
				fourSchemePoints[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else
			{
				fourSchemeValueX += (xValues[i - 1] * (-1.0 / 16.0));
				fourSchemeValueY += (yValues[i - 1] * (-1.0 / 16.0));
				fourSchemeValueX += (xValues[i] * (9.0 / 16.0));
				fourSchemeValueY += (yValues[i] * (9.0 / 16.0));
				fourSchemeValueX += (xValues[i + 1] * (9.0 / 16.0));
				fourSchemeValueY += (yValues[i + 1] * (9.0 / 16.0));
				fourSchemeValueX += (xValues[i + 2] * (-1.0 / 16.0));
				fourSchemeValueY += (yValues[i + 2] * (-1.0 / 16.0));
				fourSchemePoints[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
	}

	//Second
	fourEnd = 2 * numControlPoints;
	fourSchemePoints2 = new vec4[4 * numControlPoints];
	for (int i = 0; i < fourEnd; i++)
	{
		fourSchemePoints2[2 * i] = vec4(fourSchemePoints[i].x, fourSchemePoints[i].y, 0, 1);
			
			double fourSchemeValueX = 0, fourSchemeValueY = 0;
			if (i == 0)
			{
				fourSchemeValueX += (fourSchemePoints[fourEnd - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[fourEnd - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[0].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[0].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[2].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[2].y * (-1.0 / 16.0));
				fourSchemePoints2[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else if (i == (fourEnd - 1))
			{
				fourSchemeValueX += (fourSchemePoints[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[0].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[0].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[1].y * (-1.0 / 16.0));
				fourSchemePoints2[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else if (i == (fourEnd - 2))
			{
				fourSchemeValueX += (fourSchemePoints[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[i + 1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i + 1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[0].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[0].y * (-1.0 / 16.0));
				fourSchemePoints2[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else
			{
				fourSchemeValueX += (fourSchemePoints[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[i + 1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i + 1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints[i + 2].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints[i + 2].y * (-1.0 / 16.0));
				fourSchemePoints2[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
	}
	
	//Third
	fourEnd = 4 * numControlPoints;
	fourSchemePoints3 = new vec4[8 * numControlPoints];
	for (int i = 0; i < fourEnd; i++)
	{
			fourSchemePoints3[2 * i] = vec4(fourSchemePoints2[i].x, fourSchemePoints2[i].y, 0, 1);
			
			double fourSchemeValueX = 0, fourSchemeValueY = 0;
			if (i == 0)
			{
				fourSchemeValueX += (fourSchemePoints2[fourEnd - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[fourEnd - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[0].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[0].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[2].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[2].y * (-1.0 / 16.0));
				fourSchemePoints3[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else if (i == (fourEnd - 1))
			{
				fourSchemeValueX += (fourSchemePoints2[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[0].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[0].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[1].y * (-1.0 / 16.0));
				fourSchemePoints3[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else if (i == (fourEnd - 2))
			{
				fourSchemeValueX += (fourSchemePoints2[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[i + 1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i + 1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[0].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[0].y * (-1.0 / 16.0));
				fourSchemePoints3[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else
			{
				fourSchemeValueX += (fourSchemePoints2[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[i + 1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i + 1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints2[i + 2].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints2[i + 2].y * (-1.0 / 16.0));
				fourSchemePoints3[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
	}

	//Four
	fourEnd = 8 * numControlPoints;
	fourSchemePoints4 = new vec4[16 * numControlPoints];
	for (int i = 0; i < fourEnd; i++)
	{
			fourSchemePoints4[2 * i] = vec4(fourSchemePoints3[i].x, fourSchemePoints3[i].y, 0, 1);
			
			double fourSchemeValueX = 0, fourSchemeValueY = 0;
			if (i == 0)
			{
				fourSchemeValueX += (fourSchemePoints3[fourEnd - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[fourEnd - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[0].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[0].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[2].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[2].y * (-1.0 / 16.0));
				fourSchemePoints4[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else if (i == (fourEnd - 1))
			{
				fourSchemeValueX += (fourSchemePoints3[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[0].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[0].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[1].y * (-1.0 / 16.0));
				fourSchemePoints4[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else if (i == (fourEnd - 2))
			{
				fourSchemeValueX += (fourSchemePoints3[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[i + 1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i + 1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[0].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[0].y * (-1.0 / 16.0));
				fourSchemePoints4[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
			else
			{
				fourSchemeValueX += (fourSchemePoints3[i - 1].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i - 1].y * (-1.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[i].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[i + 1].x * (9.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i + 1].y * (9.0 / 16.0));
				fourSchemeValueX += (fourSchemePoints3[i + 2].x * (-1.0 / 16.0));
				fourSchemeValueY += (fourSchemePoints3[i + 2].y * (-1.0 / 16.0));
				fourSchemePoints4[2 * i + 1] = vec4(fourSchemeValueX, fourSchemeValueY, 0, 1);
			}
	}
	
	fourEnd = 32 * numControlPoints;
	fourSchemeLines = new vec4[32 * numControlPoints];
	int count = 0, spot = 1;
	for (int i = 0; i < fourEnd - 2; i+=2)
	{
		if (i == fourEnd - 4)
		{
			fourSchemeLines[i] = fourSchemePoints4[count];
			fourSchemeLines[i + 1] = fourSchemePoints4[0];
			count++;
			spot++;
		}
		else
		{
			fourSchemeLines[i] = fourSchemePoints4[count];
			fourSchemeLines[i + 1] = fourSchemePoints4[spot];
			count++;
			spot++;
		}
	}
	
}

void createControls()
{
	controlPoints = new vec4[numControlPoints];
	for (int i = 0; i < numControlPoints; i++)
	{
		controlPoints[i] = vec4(xValues[i], yValues[i], 0, 1);
	}
	controlLines = new vec4[2 * numControlPoints];
	for (int i = 0; i < 2 * numControlPoints; i++)
	{
		if (i < numControlPoints)
		{
			controlLines[i] = vec4(xValues[i], yValues[i], 0, 1);
		}
		else if (i == (2 * numControlPoints) - 1)
		{
			controlLines[i] = vec4(xValues[0], yValues[0], 0, 1);
		}
		else
		{
			controlLines[i] = controlLines[(i + 1) - numControlPoints];
		}
	}
}

void createAxis()
{
	arrayAxis = new vec4[4];
	// X
	arrayAxis[0][0] = RANGE_MIN;
	arrayAxis[0][1] = 0;
	arrayAxis[0][2] = 0;
	arrayAxis[0][3] = 1;
	arrayAxis[1][0] = RANGE_MAX;
	arrayAxis[1][1] = 0;
	arrayAxis[1][2] = 0;
	arrayAxis[1][3] = 1;
	// Y
	arrayAxis[2][0] = 0;
	arrayAxis[2][1] = RANGE_MIN;
	arrayAxis[2][2] = 0;
	arrayAxis[2][3] = 1;
	arrayAxis[3][0] = 0;
	arrayAxis[3][1] = RANGE_MAX;
	arrayAxis[3][2] = 0;
	arrayAxis[3][3] = 1;
}

void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		leftButtonDown = 1;
		valueX = interpolate(0, x, WIN_WIDTH, -1, 1);
		valueY = interpolate(0, y, WIN_HEIGHT, 1, -1);
		for (int i = 0; i < numControlPoints; i++)
		{
			if (i == 0)
			{
				differenceValue = pow(valueX - xValues[i], 2) + pow(valueY - yValues[i], 2);
				indexSpot = i;
			}
			else
			{
				tempDifferenceValue = pow(valueX - xValues[i], 2) + pow(valueY - yValues[i], 2);
				if (differenceValue > tempDifferenceValue)
				{
					differenceValue = tempDifferenceValue;
					indexSpot = i;
				}
			}
		}
	}
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		leftButtonDown = 0;
	}
}

void mouseMotion(int x, int y)
{
	if (leftButtonDown == 1)
	{
		xValues[indexSpot] = interpolate(0, x, WIN_WIDTH, -1, 1);
		yValues[indexSpot] = interpolate(0, y, WIN_HEIGHT, 1, -1);
		createControls();
		createBernstein();
		createQuadBSpline();
		createCubicBSpline();
		createFourScheme();
		glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis, sizeOfControlPoints, controlPoints);
		glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints, sizeOfControlLines, controlLines);
		glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines, sizeOfBernsteinLines, bernsteinLines);
		glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines + sizeOfBernsteinLines, sizeOfQuadLines, quadLines);
		glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines + sizeOfBernsteinLines + sizeOfQuadLines, sizeOfCubicLines, cubicLines);
		glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines + sizeOfBernsteinLines + sizeOfQuadLines + sizeOfCubicLines, sizeOfFourSchemeLines, fourSchemeLines);
	}
	glutPostRedisplay();
}

void keyboard(unsigned char key, int width, int height)
{
	switch( key ) 
	{
	case 033:  // Escape key
	case 'q': case 'Q':
	    exit( EXIT_SUCCESS );
	    break;
	}
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//Axis
	vColor = vec4(0, 0, 0, 0);
	glUniform4fv(vColorID, 1, vColor);
	glDrawArrays(GL_LINES, 0, 4);
	//Control Points
	vColor = vec4(0, 0, 1, 1);
	glUniform4fv(vColorID, 1, vColor);
	glDrawArrays(GL_POINTS, 4, numControlPoints);
	//Control Box
	vColor = vec4(0, 0, 1, 1);
	glUniform4fv(vColorID, 1, vColor);
	glDrawArrays(GL_LINES, 4 + numControlPoints, (2 * numControlPoints));
	//Bernstein Curve
	vColor = vec4(1, 0, 0, 1);
	glUniform4fv(vColorID, 1, vColor);
	glDrawArrays(GL_LINES, (4 + numControlPoints) + ((2 * numControlPoints)), (2 * NUM_POINTS) - 2);
	//Quad Curve
	vColor = vec4(0, 1, 0, 1);
	glUniform4fv(vColorID, 1, vColor);
	glDrawArrays(GL_LINES, (4 + numControlPoints) + ((2 * numControlPoints)) + ((2 * NUM_POINTS) - 2), (2 * NUM_POINTS) - 2);
	//Cubic Curve
	vColor = vec4(0, 1, 1, 1);
	glUniform4fv(vColorID, 1, vColor);
	glDrawArrays(GL_LINES, (4 + numControlPoints) + ((2 * numControlPoints)) + ((2 * NUM_POINTS) - 2) + ((2 * NUM_POINTS) - 2), (2 * (numControlPoints * NUM_POINTS)));
	//Four Curve
	vColor = vec4(1, 0, 1, 1);
	glUniform4fv(vColorID, 1, vColor);
	glDrawArrays(GL_LINES, (4 + numControlPoints) + ((2 * numControlPoints)) + ((2 * NUM_POINTS) - 2) + ((2 * NUM_POINTS) - 2) + ((2 * (numControlPoints * NUM_POINTS))), (32 * numControlPoints));
    glFlush();
}

void init()
{   
    // Vertex array object
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Buffer object
    GLuint buffer;
    glGenBuffers(1, &buffer );
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
	
    // Empty buffer
	glBufferData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines + sizeOfBernsteinLines + sizeOfQuadLines + sizeOfCubicLines + sizeOfFourSchemeLines, NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeOfAxis, arrayAxis);
	glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis, sizeOfControlPoints, controlPoints);
	glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints, sizeOfControlLines, controlLines);
	glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines, sizeOfBernsteinLines, bernsteinLines);
	glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines + sizeOfBernsteinLines, sizeOfQuadLines, quadLines);
	glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines + sizeOfBernsteinLines + sizeOfQuadLines, sizeOfCubicLines, cubicLines);
	glBufferSubData(GL_ARRAY_BUFFER, sizeOfAxis + sizeOfControlPoints + sizeOfControlLines + sizeOfBernsteinLines + sizeOfQuadLines + sizeOfCubicLines, sizeOfFourSchemeLines, fourSchemeLines);

    // Load shaders and use the resulting shader program
    GLuint program = InitShader("vshader.glsl", "fshader.glsl");
    glUseProgram(program);

	vColorID = glGetUniformLocation(program, "vColor");

    // Initialize the vertex position attribute from the vertex shader    
    GLuint vPosition = glGetAttribLocation(program, "vPosition");
    glEnableVertexAttribArray(vPosition);
    glVertexAttribPointer(vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));

    glEnable(GL_DEPTH_TEST);
    glClearColor(1.0, 1.0, 1.0, 1.0);
	glPointSize(5);
	glLineWidth(2.5);
}

int main(int argc, char **argv)
{
	string filename = "ControlPoints.txt";
	readFile(filename);
	
	createAxis();
	createControls();
	createBernstein();
	createQuadBSpline();
	createCubicBSpline();
	createFourScheme();

	sizeOfAxis = 4 * sizeof(vec4);
	sizeOfControlPoints = numControlPoints * sizeof(vec4);
	sizeOfControlLines = (2 * numControlPoints) * sizeof(vec4);
	sizeOfBernsteinLines = (2 * NUM_POINTS - 2) * sizeof(vec4);
	sizeOfQuadLines = (2 * NUM_POINTS - 2) * sizeof(vec4);
	sizeOfCubicLines = (2 * (numControlPoints * NUM_POINTS)) * sizeof(vec4);
	sizeOfFourSchemeLines = (32 * numControlPoints) * sizeof(vec4);

	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(WIN_WIDTH, WIN_HEIGHT);
    glutCreateWindow("Geometric Modeling");

	glewExperimental = GL_TRUE;
	glewInit();
	init();	

    glutDisplayFunc(display);
    glutMouseFunc(mouse);
	glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);

    glutMainLoop();
	return 0;
}