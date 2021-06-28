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




//Структура особи
struct Being
{
	unsigned long long numFenotype;//Номер фенотипа
	unsigned long long int numInPopulation;//Номер особи в популяции
	double X;
	std::string chromosome;//Хромосома особи
	//double Y;
	double fx;//Значение целевой функции
	double Fx;//Значение функции приспособленности
};

template <class T>
void PrintVector(std::vector<T>&prV)
{
	for (auto it = prV.begin(); it != prV.end(); it++)
	{
		std::cout << *it << std::endl;
	}
}


//Функция Sign
template <class Value>
int Sign(Value val)
{
	if (val == 0) return 0;
	if (val > 0) return 1;
	else return -1;
}


//Вывод значений в график
void ToGrathic(std::string str);


//Структура особи
struct Being;

//Класс генетического алгоритма
class GeneticAlgoritm
{
private:
	double Pc;//Вероятность скрещивания
	double Pm;//Вероятность мутации
	double dF;//Величина приращения. Если приращение приспособленности популяции меньше ее, останавливаем
	unsigned long long int MAX;//Максимальное количество итераций.
	int populationCount;//Популяция
	int genCount;//Количество генов в геноме каждой особи
	double xMin, xMax;//Минимум и максимум по Х
	unsigned long long int D;//Максимальное количество координат по X
	double stepX; //Шаг по X
	unsigned long long numPopulation;//Номер текущей популяции
	std::vector<Being> population;//Популяция
	double prevFs = 100;//Запоминаем предыдущее значение приращения для остановки

	std::ofstream out_population;



/* ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ */
//Фукнция расчета X и соответствующего ему номера позиции
public:

	//Запись в файл
	void WriteToFile(std::string str);

	//Вывод особи на экран
	void printBeing(Being being);

	//Запись популяции в файл
	void WritePopulation();
	
	//Расчет Fs
	double GetFsPopulation(std::vector<Being>& BeingsInPopulation);

	//Считает X по номеру позиции
	double ReturnXToNumber(unsigned long long numberX);
	
	//Функция приспособленности
	double ResultFitness(double fx);

	//Возвращает случайное число в интервале от 0 до N
	unsigned long long int returnTheProbabiliryInt(unsigned long long N);

	//Возврат случайного double в интервале от 0 до 1
	double returnTheProbabilityDouble();


	//Функция расчета значения заданной функции
	double ResultFunction(double X);

	//Функция перевода из десятичной в двоичную
	unsigned long long int DecimalToBinary(unsigned long long decimal);

	//Функция перевода из двоичной строки в десятичное число
	unsigned long long int BinaryToDecimal(std::string binaryStr);

	//Функция перевода двоичного числа в строку
	std::string BinaryToBinaryStr(unsigned long long int binary);
	

	/* ----------------------- */

	/* ФУНКЦИИ СКРЕЩИВАНИЯ, МУТАЦИИ, РЕДУКЦИИ */
	//Формирование начальной популяции	
	void GetStartPopulation();

	//Скрещивание случайным отбором
	void Breding();
	//Мутация точечная
	void Mutation(std::vector<std::string>& populationGen);
	//Редукция рулетка
	void Reduction_Roulette(std::vector<std::string>& populationGen);
	
	//Функция остановки
	void Stop();
	
	//Вывод значений на графики на основе популяций
	void AddValuesToGraprhic();

	GeneticAlgoritm();
	~GeneticAlgoritm()
	{
		out_population.close();
	}
};

#endif // __GA_H__