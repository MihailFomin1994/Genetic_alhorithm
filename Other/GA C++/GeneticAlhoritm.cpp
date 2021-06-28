#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <random>
#include <iomanip>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>
#include "GeneticAlhoritm.h"


std::random_device rd;
std::mt19937 gen;
auto rnd = std::default_random_engine{};

//����� ����� �� �����
void GeneticAlgoritm::printBeing(Being being)
{
	std::cout << std::fixed;
	std::cout.precision(4);
	std::cout << "Num in population: " << being.numInPopulation << "\tX = " << being.X << "\tChromosome:" << being.chromosome << "\tf(x) = " << being.fx << "\tF(x) = " << being.Fx << std::endl;
}

//������� �����������������
double GeneticAlgoritm::ResultFitness(double fx)
{
	double Fx = fx + 45;//���� ������� �����������������
	return Fx;
}


//������� ������� �������� �������� �������
double GeneticAlgoritm::ResultFunction(double X)
{
	double result = -1 * exp(-3.4 * X) * Sign(sin(43 * X));
	return result;
}

//������� �������� ��������� ����� � ������
std::string GeneticAlgoritm::BinaryToBinaryStr(unsigned long long int binary)
{
	std::string binaryStr = std::to_string(binary);
	if (binaryStr.length() < static_cast<unsigned int>(genCount))
	{
		int countZero = genCount - binaryStr.length();
		std::string zeroStr = "";
		for (int i = 0; i < countZero; i++)
		{
			zeroStr += '0';
		}
		binaryStr = zeroStr + binaryStr;
	}
	return binaryStr;
}


//������� �������� �� �������� ������ � ���������� �����
unsigned long long int GeneticAlgoritm::BinaryToDecimal(std::string binaryStr)
{
	//������� ���������� ���������� ����
	binaryStr.erase(0, std::min(binaryStr.find_first_not_of('0'), binaryStr.size() - 1));
	if (binaryStr.length() == 0) binaryStr = "0";//���� ����� ����� ������� �� �����
	int value;
	size_t length;
	value = 0;
	length = binaryStr.length();
	for (size_t i = 0; i < length; i++)
	{
		value |= (binaryStr[i] == '1') ? (1 << (length - i - 1)) : 0;
	}
	return value;
}



//������� �������� �� ���������� � ��������
unsigned long long int GeneticAlgoritm::DecimalToBinary(unsigned long long decimal)
{
	unsigned long long int binary = 0, k = 1;
	while (decimal)
	{
		binary += (decimal % 2) * k;
		k *= 10;
		decimal /= 2;
	}
	return binary;
}


//���������� ��������� ����� � ��������� �� 0 �� N
unsigned long long int GeneticAlgoritm::returnTheProbabiliryInt(unsigned long long N)
{
	N = N - 1;
	std::uniform_int_distribution<> dis(0, N);
	return dis(gen);
}

//������� ���������� double � ��������� �� 0 �� 1
double GeneticAlgoritm::returnTheProbabilityDouble()
{
	
	std::uniform_real_distribution<> dis(0.01, 0.99);


	double tmp = dis(gen);

	return tmp;
}





//������ Fs
double GeneticAlgoritm::GetFsPopulation(std::vector<Being>& BeingsInPopulation)
{
	double sum = 0;
	for (auto it = BeingsInPopulation.begin(); it != BeingsInPopulation.end(); it++)
	{
		sum += it->Fx;
	}
	return sum;
}

//������ � ����
void GeneticAlgoritm::WriteToFile(std::string str)
{
	if (out_population.is_open())
	{
		out_population << str << std::endl;
	}
	
}
//����� ��������� �� �������
void AddValuesToGraprhic(std::string str)
{
	std::ofstream out;
	out.open("Graphic_f.txt", std::ios::app);
	if (out.is_open())
	{
		out << str << std::endl;
	}
	out.close();
}


//����� �������� � ������
void ToGrathic(std::string str)
{
	std::ofstream out;
	out.open("Grathic.txt", std::ios::app);
	if (out.is_open())
	{
		out << str << std::endl;
	}
	out.close();
}

//������� X �� ������ �������
double GeneticAlgoritm::ReturnXToNumber(unsigned long long numberX)
{
	return xMin + stepX * numberX;
}

//������ ��������� � ����
void GeneticAlgoritm::WritePopulation()
{
	WriteToFile("Num population: " + std::to_string(numPopulation));
	for (auto it = population.begin(); it != population.end(); it++)
	{
		std::string fx = std::to_string(it->fx);
		std::string tmp = "\t";
		if (fx.size() < 9)
		{
			tmp += "\t";
			
		}
		fx = fx + tmp;
		std::string str = "Num in population: " + std::to_string(it->numInPopulation) + "\tX = " + std::to_string(it->X) + "\tChromosome:" + it->chromosome + "\tf(x) = " + fx + "F(x) = " + std::to_string(it->Fx);
		WriteToFile(str);
	}

}


//����������� �� ���������
GeneticAlgoritm::GeneticAlgoritm()
{
	dF = 0.001;
	MAX = 100;


	numPopulation = 0;
	xMin = -1;
	xMax = 1;
	genCount = 20;
	D = pow(2, genCount);
	populationCount = 50;
	stepX = (xMax - xMin) / (static_cast<double>(D) - 1);
	Pc = 0.85;
	Pm = 0.15;

	out_population.open("Populations.txt", std::ios::app);

	//�����������
	GetStartPopulation();

}

//������������ ��������� ���������	
void GeneticAlgoritm::GetStartPopulation()
{
	std::set<unsigned long long> startPopulation;//������������ ��� ��������
	unsigned long long i = 0;
	while (i < populationCount)
	{
		long long num = returnTheProbabiliryInt(D);
		if (startPopulation.find(num) == startPopulation.end())
		{
			i++;
			startPopulation.insert(num);
			double X = ReturnXToNumber(num);
			double fx = ResultFunction(X);
			unsigned long long bin_hromosome = DecimalToBinary(num);
			std::string str_hromosom = BinaryToBinaryStr(bin_hromosome);
			Being being;
			being.X = X;
			being.chromosome = str_hromosom;
			being.fx = fx;
			being.numFenotype = num;
			being.numInPopulation = i;
			being.Fx = ResultFitness(fx);
			population.push_back(being);
		}
	}
	startPopulation.clear();
	numPopulation++;
	WritePopulation();//������� ������ ��������� � ����. ����� �������������� � ������� � ������. 
	Breding();//�����������
}

//����������� ��������� �������
void GeneticAlgoritm::Breding()
{
	//��� �������� ������������� ������� �� ��������� ������ � ������� � ������
	std::vector<std::string>populationGens;
	for (auto it = population.begin(); it != population.end(); it++)
	{
		populationGens.push_back(it->chromosome);
	}
	
	std::shuffle(populationGens.begin(), populationGens.end(), rnd);//������������� �������


	std::vector<std::string> newPopulation;
	//�����������. ����� ��������������� 2 ��������. ����� ���� ����������� ����������� ���� �������������, ���������� ��������, ���� ���� - ����� ���������
	for (auto it = populationGens.begin(); it != populationGens.end(); it = it + 2)
	{
		auto mather = it;
		auto father = it + 1;
		double P = returnTheProbabilityDouble();

		if (P < Pc)//���� ��������� ����� ������ ����������� �����������
		{
			int crossingoverPoint = static_cast<int>(returnTheProbabiliryInt(genCount));

			if ((crossingoverPoint != 0) && (crossingoverPoint != (genCount - 1)))
			{
				newPopulation.push_back(mather->substr(0, crossingoverPoint) + father->substr(crossingoverPoint, father->size() - 1));
				newPopulation.push_back(father->substr(0, crossingoverPoint) + mather->substr(crossingoverPoint, father->size() - 1));
			}
			else
			{
				newPopulation.push_back(*mather);
				newPopulation.push_back(*father);
			}
		}
		else
		{
			newPopulation.push_back(*mather);
			newPopulation.push_back(*father);
		}
	}
	//������������ ������
	populationGens.clear();
	population.clear();
	populationGens.shrink_to_fit();
	population.shrink_to_fit();

	//����� �������� �������
	Mutation(newPopulation);
}




//������� ��������
void GeneticAlgoritm::Mutation(std::vector<std::string>& populationGen)
{
	std::vector<std::string> mutationGen;
	for (auto code = populationGen.begin(); code != populationGen.end(); code++)
	{
		double P = returnTheProbabilityDouble();
		if (P < Pm)//���� ������� ����������
		{
			int allele = static_cast<int>(returnTheProbabiliryInt(genCount));
			std::string gens = *code;
			if (gens[allele] == '1')
			{
				gens[allele] = '0';
			}
			else
			{
				gens[allele] = '1';
			}
			
			mutationGen.push_back(gens);
		}
	}
	populationGen.insert(populationGen.end(), mutationGen.begin(), mutationGen.end());
	mutationGen.clear();
	mutationGen.shrink_to_fit();
	Reduction_Roulette(populationGen);//������ �������� "�������"
}

//�������� �������
void GeneticAlgoritm::Reduction_Roulette(std::vector<std::string>& populationGen)
{
	std::vector<Being> allBeingsInPopulation;
	//���������� ���� ������ � ���������, ��������� ���


	int i = 0;
	for (auto it = populationGen.begin(); it != populationGen.end(); it++)
	{
		i++;//���������
		Being* being = new Being;
		being->chromosome = *it;
		being->numFenotype = BinaryToDecimal(*it);
		being->X = ReturnXToNumber(being->numFenotype);

		//����� � ��������� ���� ���������� �� �������� (�� �����������)
		being->numInPopulation = i;
		/// ///////////////////////////////////////////////////////

		being->fx = ResultFunction(being->X);
		being->Fx = ResultFitness(being->fx);
		allBeingsInPopulation.push_back(*being);
	}
	int numInPopulation = 1;
	//���� � �������. 
	while (population.size() < populationCount)
	{
		double Fs = GetFsPopulation(allBeingsInPopulation);
		double max = 0;
		std::vector<double> roulette;
		for (auto it = allBeingsInPopulation.begin(); it != allBeingsInPopulation.end(); it++)
		{
			max += it->Fx / Fs;
			roulette.push_back(max);
		}
		double P = returnTheProbabilityDouble();//�������� ��������� �������� �� 0 �� 1
		int winter = 0;
		for (int i = 0; i < roulette.size(); i++)
		{
			if (roulette[i] > P)
			{
				allBeingsInPopulation[i].numInPopulation = numInPopulation;
				population.push_back(allBeingsInPopulation[i]);//���������� � ������� ��������� (��� ������) �������� ����� 
				allBeingsInPopulation.erase(allBeingsInPopulation.begin() + i);
				numInPopulation++;
				break;
			}
		}
		roulette.clear();
		roulette.shrink_to_fit();
		
	}
	allBeingsInPopulation.clear();
	allBeingsInPopulation.shrink_to_fit();
	numPopulation++;//����������� ����� ���������
	WritePopulation();


	populationGen.clear();
	
	populationGen.shrink_to_fit();

	Stop();
}

//������� ��������� ��� ���������� �������
void GeneticAlgoritm::Stop()
{
	double Fs = GetFsPopulation(population);
	double population_dF = (Fs - prevFs) / Fs;
	if (numPopulation >= MAX || population_dF<dF)
	{
		std::cout << numPopulation << std::endl;
		system("pause");
		std::cout << "Maximux" << std::endl;
	}
	else
	{
		prevFs = Fs;
		Breding();
	}
}
