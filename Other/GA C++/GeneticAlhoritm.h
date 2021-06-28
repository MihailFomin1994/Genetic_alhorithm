#ifndef __GA_H__
#define __GA_H__
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <random>
#include <iomanip>
#include <vector>
#include <cmath>
#include <set>




//��������� �����
struct Being
{
	unsigned long long numFenotype;//����� ��������
	unsigned long long int numInPopulation;//����� ����� � ���������
	double X;
	std::string chromosome;//��������� �����
	//double Y;
	double fx;//�������� ������� �������
	double Fx;//�������� ������� �����������������
};

template <class T>
void PrintVector(std::vector<T>&prV)
{
	for (auto it = prV.begin(); it != prV.end(); it++)
	{
		std::cout << *it << std::endl;
	}
}


//������� Sign
template <class Value>
int Sign(Value val)
{
	if (val == 0) return 0;
	if (val > 0) return 1;
	else return -1;
}


//����� �������� � ������
void ToGrathic(std::string str);


//��������� �����
struct Being;

//����� ������������� ���������
class GeneticAlgoritm
{
private:
	double Pc;//����������� �����������
	double Pm;//����������� �������
	double dF;//�������� ����������. ���� ���������� ����������������� ��������� ������ ��, �������������
	unsigned long long int MAX;//������������ ���������� ��������.
	int populationCount;//���������
	int genCount;//���������� ����� � ������ ������ �����
	double xMin, xMax;//������� � �������� �� �
	unsigned long long int D;//������������ ���������� ��������� �� X
	double stepX; //��� �� X
	unsigned long long numPopulation;//����� ������� ���������
	std::vector<Being> population;//���������
	double prevFs = 100;//���������� ���������� �������� ���������� ��� ���������

	std::ofstream out_population;



/* ��������������� ������� */
//������� ������� X � ���������������� ��� ������ �������
public:

	//������ � ����
	void WriteToFile(std::string str);

	//����� ����� �� �����
	void printBeing(Being being);

	//������ ��������� � ����
	void WritePopulation();
	
	//������ Fs
	double GetFsPopulation(std::vector<Being>& BeingsInPopulation);

	//������� X �� ������ �������
	double ReturnXToNumber(unsigned long long numberX);
	
	//������� �����������������
	double ResultFitness(double fx);

	//���������� ��������� ����� � ��������� �� 0 �� N
	unsigned long long int returnTheProbabiliryInt(unsigned long long N);

	//������� ���������� double � ��������� �� 0 �� 1
	double returnTheProbabilityDouble();


	//������� ������� �������� �������� �������
	double ResultFunction(double X);

	//������� �������� �� ���������� � ��������
	unsigned long long int DecimalToBinary(unsigned long long decimal);

	//������� �������� �� �������� ������ � ���������� �����
	unsigned long long int BinaryToDecimal(std::string binaryStr);

	//������� �������� ��������� ����� � ������
	std::string BinaryToBinaryStr(unsigned long long int binary);
	

	/* ----------------------- */

	/* ������� �����������, �������, �������� */
	//������������ ��������� ���������	
	void GetStartPopulation();

	//����������� ��������� �������
	void Breding();
	//������� ��������
	void Mutation(std::vector<std::string>& populationGen);
	//�������� �������
	void Reduction_Roulette(std::vector<std::string>& populationGen);
	
	//������� ���������
	void Stop();
	
	//����� �������� �� ������� �� ������ ���������
	void AddValuesToGraprhic();

	GeneticAlgoritm();
	~GeneticAlgoritm()
	{
		out_population.close();
	}
};

#endif // __GA_H__