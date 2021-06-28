using System;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.Diagnostics;
using System.Windows.Forms.DataVisualization.Charting;
using System.ComponentModel;
using System.Drawing;
using Genetic_Algorihm;


namespace GA
{
    public static class RandomProvider
    {
        
        private static int seed = Environment.TickCount;
        private static ThreadLocal<Random> randomWrapper = new ThreadLocal<Random>(() =>
        new Random(Interlocked.Increment(ref seed)));
        public static Random GetThreadRandom()
        {
            return randomWrapper.Value;
        }
    }
    //Считает целевую функцию
    public static class Function
    {
        //Считает целевую функцию f(x,y)
        public static double Get_fxy(double x, double y = 0)
        {
            //f(x,y) = 0.1*x+0.1*y-4*cos(0.8*x)+4*cos(0.8*y)+8
            return 0.1 * x + 0.1 * y - 4 * Math.Cos(0.8 * x) + 4 * Math.Cos(0.8 * y) + 8;
        }
        //Считает функцию приспособленности F(x,y)
        public static double Get_Fxy(double x, double y)
        {
            //F(x,y) = f(x,y)+1;
            return 1 + Get_fxy(x, y);
        }
    }

    public struct Specimen
    {
        public string X;
        public string Y;
    }



    public class GeneticAlgorihm
    {

        private int N = 20;//Количество генов в хромосоме
        private int populationSize;//Размер популяции
        private int D;//Количество фенотипов
        private double Pc, Pm;//Вероятность скрещивания и мутации
        private double minSearch;
        private double maxSearch;
        private double step;//Шаг в области поиска
        private double dF_lim;//Приращение приспособленности
        private double iterationMAX;//Максимальное количество итераций

        private int populationNumber = 0;//Номер текущей популяции
        private List<Specimen> population = new List<Specimen>();//Текущая популяция
        private double prevFS = 10000;//Для подсчета приращения. Может быть любое число, лишь бы не останавливало на первой итерации

        private Random rnd = RandomProvider.GetThreadRandom();//Определяет рандом
        private bool twoVariable = false;//Использование функции двух переменных
        private Chart graphic_cel;


        private Color calculateColor(double min, double z, double max)
        {
            int c = (int)Math.Round(0 + (z - min) * (255 - 0) / (max - min));
            return Color.FromArgb(255, c, c);
        }


        //Конструктор
        public GeneticAlgorihm(Chart chart1, int N = 20, int populationSize = 50, double Pc = 0.85, double Pm = 0.15, double minSearch = -1, double maxSearch = 1, double dF_lim = 0.01, int iterationMAX = 100)
        {
            this.N = N;
            this.populationSize = populationSize;
            D = Convert.ToInt32(Math.Pow(2, N));
            this.Pc = Pc;
            this.Pm = Pm;
            this.minSearch = minSearch;
            this.maxSearch = maxSearch;
            this.step = (Math.Abs(maxSearch) + Math.Abs(minSearch)) / this.D;
            this.dF_lim = dF_lim;
            this.iterationMAX = iterationMAX;



            //Инициализация графиков
            this.graphic_cel = chart1;

            graphic_cel.Legends["Legend1"].Enabled = false;
            graphic_cel.ChartAreas["ChartArea1"].Area3DStyle.Enable3D = true;
            graphic_cel.ChartAreas["ChartArea1"].AxisX.Enabled = AxisEnabled.True;
            graphic_cel.ChartAreas["ChartArea1"].AxisY.Enabled = AxisEnabled.True;
            graphic_cel.ChartAreas["ChartArea1"].Area3DStyle.Inclination = 30;

            graphic_cel.Series.RemoveAt(0);


            
            double x, y, z;
            double x1 = -6, x2 = 6, y1 = -6, y2 = 6;
            int countPoints = 10;
            double xstep = 12 / (countPoints);
            double ystep = 12 / (countPoints);
            double max = Function.Get_fxy(x1, y1), min = max;
            Color Color1;
            for (x = x1; x <= x2; x += xstep)
                for (y = y1; y <= y2; y += ystep)
                {
                    z = Function.Get_fxy(x, y);
                    if (z > max) max = z;
                    else if (z < min) min = z;
                }
            
            for (int i = 0; i < countPoints; i++)
            {
                Series Series1 = new Series();
                Series1.ChartType = SeriesChartType.Spline;
                graphic_cel.Series.Add(Series1);
                for (int j = 0; j < countPoints; j++)
                {
                    x = x1 + i * xstep;
                    y = y1 + j * ystep;
                    z = Function.Get_fxy(x, y);
                    Color1 = calculateColor(min, z, max);
                    //graphic_cel.Series[i].Points.AddXY((double)j, z);
                    graphic_cel.Series[i].Points.AddXY(x, y, z);
                    
                    graphic_cel.Series[i].Points[j].Color = Color1;
                    graphic_cel.Series[i].Points[j].BackSecondaryColor = Color1;
                }



            }
            


            




            InitialPopulation();
        }

        //Считает целевую функцию f(x,y)
        public double Get_fxy(Specimen specimen)
        {
            double x = Convert.ToDouble(Chrimosome_Decoding(specimen.X)) * this.step;
            double y = Convert.ToDouble(Chrimosome_Decoding(specimen.Y)) * this.step;
            return Function.Get_fxy(x, y);
        }
        //Считает функцию приспособленности F(x,y)
        public double Get_Fxy(Specimen specimen)
        {
            double x = Convert.ToDouble(Chrimosome_Decoding(specimen.X)) * this.step;
            double y = Convert.ToDouble(Chrimosome_Decoding(specimen.Y)) * this.step;
            return Function.Get_Fxy(x, y);
        }
        //Считает суммарную приспосабливаемость
        public double Get_FS(List<Specimen> population)
        {
            double Fs = 0;
            foreach(var sp in population)
            {
                Fs = Get_Fxy(sp);
            }
            return Fs;
        }


        //Инициализация начальной популяции (она одинаковая для всех)
        public void InitialPopulation()
        {
            List<int> zeroPopulationX = new List<int>();
            List<int> zeroPopulationY = new List<int>();
            zeroPopulationX = GetRandomNumbers(this.N, this.D);
            if(twoVariable)
            {
                zeroPopulationY = GetRandomNumbers(this.N, this.D);
            }
            else
            {
                zeroPopulationY = Enumerable.Repeat(0, zeroPopulationX.Count()).ToList();
            }
            //Заполнение X
            int count = 0;
            foreach(int num in zeroPopulationX)
            {
                Specimen sp = new Specimen();
                sp.X = Chromosome_Coding(num);
                sp.Y = Chromosome_Coding(zeroPopulationY[count]);
                population.Add(sp);
                count++;
            }
            
        }

        //Скрещивание
        public void Crossbreeding()
        {
            Dictionary<Specimen, Specimen> parentPairs = FormatioOfParentPairs_Random();//Используется случайный отбор

            //Далее для каждой пары смотрим вероятность скрещивания
            List<Specimen> newPopulation = new List<Specimen>();//Формировние новой популяции
            int crossingPoint;
            double P;
            foreach (var pairs in parentPairs)
            {
                P = GetRandom_Double();
                Specimen ch1, ch2;
                if(P<this.Pc)//Тогда создаем потомков
                {
                    crossingPoint = GetRandom_Int(max:this.N - 1);
                    string child1_X = GetChild_1(pairs.Key.X, pairs.Value.X, crossingPoint);
                    string child2_X = GetChild_1(pairs.Key.X, pairs.Value.X, crossingPoint);
                    string child1_Y, child2_Y;

                    if (twoVariable)
                    {
                        crossingPoint = GetRandom_Int(max: this.N - 1);
                        child1_Y = GetChild_1(pairs.Key.Y, pairs.Value.Y, crossingPoint);
                        child2_Y = GetChild_1(pairs.Key.Y, pairs.Value.Y, crossingPoint);
                    }
                    else
                    {
                        child1_Y = child2_Y = FirstZeroBinary("0");
                    }
                    ch1.X = child1_X;
                    ch1.Y = child1_Y;
                    ch2.X = child2_X;
                    ch2.Y = child2_Y;
                }
                //Иначе вносим самих родителей
                ch1 = pairs.Key;
                ch2 = pairs.Value;
            }

        }
        //Формирование родительских пар - случайный отбор
        public Dictionary<Specimen, Specimen> FormatioOfParentPairs_Random()
        {
            Dictionary<Specimen, Specimen> parentPairs = new Dictionary<Specimen, Specimen>();
            List<int> selection = new List<int>();
            selection = GetRandomNumbers(population.Count(), population.Count() - 1);//Заполнение массива случайными неповторяющимися числами
            for (int i = 0; i < selection.Count(); i+=2)
            {
                parentPairs.Add(population[i], population[i + 1]);
            }
            return parentPairs;
        }

        //Формирование родительских пар - лучшие с лучшими
        public Dictionary<Specimen, Specimen> FormatioOfParentPairs_BestWithBest()
        {
            //Не очень понимаю, как именно реализовать "лучшие с лучшими". В дальнейшем можно исправить
            Dictionary<Specimen, Specimen> parentPairs = new Dictionary<Specimen, Specimen>();
            SortedDictionary<double, Specimen> bestParent = new SortedDictionary<double, Specimen>();
            foreach(Specimen specimen in population)//Считаем значение функции приспособленности для каждой особи популяции и сортируем в порядке возрастания
            {
                bestParent.Add(Get_Fxy(specimen), specimen);
            }
            for(int i = bestParent.Count-1; i>=0; i-=2)//Заполняем родительские пары, используя сортированный ранее словарь с конца
            {
                parentPairs.Add(bestParent.ElementAt(i).Value, bestParent.ElementAt(i - 1).Value);
            }
            return parentPairs;
        }


        //Мутация
        public void Mutation(List<Specimen> newPopulation)
        {
            List<Specimen> mutPopulation = newPopulation;
            double P = GetRandom_Double();
            for(int i = 0; i<newPopulation.Count(); i++)
            {
                if(P<this.Pm)
                {
                    var tmp = newPopulation[i];
                    //Используется точечная мутация
                    tmp.X = Mutation_Point(newPopulation[i].X);
                    if (twoVariable) tmp.Y = Mutation_Point(newPopulation[i].Y);
                    mutPopulation.Add(tmp);
                }
            }
            newPopulation.Clear();
        }
        //Мутация - точечная
        public string Mutation_Point(string chromosome)
        {
            int allele = GetRandom_Int(max: this.N - 1);
            if(chromosome[allele] == '0')
            {
                chromosome = chromosome.Remove(allele, 1).Insert(allele, "1");
            }
            else
            {
                chromosome = chromosome.Remove(allele, 1).Insert(allele, "0");
            }
            return chromosome;
        }
        //Мутация - инверсия
        public string Mutation_Inversion(string chromosome)
        {
            List<int> alleles = GetRandomNumbers(2, this.N - 1);
            var gen1 = chromosome[alleles[0]].ToString();
            var gen2 = chromosome[alleles[1]].ToString();
            chromosome = chromosome.Remove(alleles[0], 1).Insert(alleles[0], gen1);
            chromosome = chromosome.Remove(alleles[1], 1).Insert(alleles[1], gen2);
            return chromosome;
        }


        //Редукция
        public void Reduction(List<Specimen> mutPopulation)
        {
            this.population.Clear();//Очищаем популяцию
            Reduction_Roulette(mutPopulation);//Редукция - рулетка

        }
        //Редукция - турнир
        public void Reduction_Tournament(List<Specimen> mutPopulation)
        {
            while(population.Count()<this.populationSize)
            {
                List<int> rivals = GetRandomNumbers(2, mutPopulation.Count() - 1);
                double rivals1_F = Get_Fxy(mutPopulation[rivals[0]]);
                double rivals2_F = Get_Fxy(mutPopulation[rivals[1]]);
                if(rivals1_F > rivals2_F)
                {
                    population.Add(mutPopulation[rivals[0]]);
                    mutPopulation.RemoveAt(rivals[0]);
                }
                else
                {
                    population.Add(mutPopulation[rivals[1]]);
                    mutPopulation.RemoveAt(rivals[1]);
                }
            }
        }

        //Редукция - рулетка
        public void Reduction_Roulette(List<Specimen> mutPopulation)
        {
            
            while(population.Count() < this.populationSize)
            {
                List<double> roulette = new List<double>();
                double begin = 0;
                double Fs = Get_FS(mutPopulation);
                foreach (var sp in population)
                {
                    roulette.Add(begin);
                    begin += Get_Fxy(sp)/Fs;

                }
                double P = GetRandom_Double();
                for(int i = 0; i<population.Count(); i++)
                {
                    if(roulette[i] > P)
                    {
                        population.Add(mutPopulation[i]);
                        mutPopulation.RemoveAt(i);
                        break;
                    }
                }
            }
        }

        //Остановка
        public void Stop()
        {
            PopulationOnGraph();
            this.populationNumber++;
            double Fs = Get_FS(this.population);
            double dF = (Fs - prevFS) / Fs;
            if(this.populationNumber < this.iterationMAX &&  dF >= this.dF_lim)
            {
                Crossbreeding();
            }
        }


        //Кодирование хромосомы
        private string Chromosome_Coding(int dec)
        {
            //Кодирование двоичным кодом
            return DecimalToBinary(dec);
        }

        //Декодирование хромосомы
        private int Chrimosome_Decoding(string code)
        {
            return BinaryToDecmal(code);
        }


        //Возвращает первого ребенка
        public string GetChild_1(string mather, string father, int crossingPoint)
        {
            return mather.Substring(0, crossingPoint) + father.Substring(crossingPoint);
        }
        //Возвращает второго ребенка
        public string GetChild_2(string mather, string father, int crossingPoint)
        {
            return father.Substring(0, crossingPoint) + mather.Substring(crossingPoint);
        }

        //Выдает разные целочисленные значения установленной размерности в диапазоне от 0 до max
        private List<int> GetRandomNumbers(int count, int max)
        {
            //Не самый лучший алгоритм. Можно оптимизировать
            List<int> randomNumbers = new List<int>();
            int number;
            while(randomNumbers.Count()<count)
            {
                number = GetRandom_Int(max: max);
                if(!randomNumbers.Contains(number))
                {
                    randomNumbers.Add(number);
                    Debug.WriteLine(number);
                }
            }
            return randomNumbers;
        }

        

        //Из десятичного в двоичное (строку) установленной длинны
        public string DecimalToBinary(int Decimal)
        {
            return FirstZeroBinary(Convert.ToString(Decimal, 2));
        }
        
        //Из двоичного в десятичное
        public int BinaryToDecmal(string Binary)
        {
            return Convert.ToInt32(Binary, 2);
        }
        
        //Заполнение лидирующими нулями
        public string FirstZeroBinary(string binary)
        {
            int zeroCount = N - binary.Length;
            if (zeroCount < 0) return binary;
            string firstZero = "";
            while(firstZero.Length < zeroCount)
            {
                firstZero += '0';
            }
            return firstZero + binary;
        }
        
        //Возвращает рандомное число с плавающей запятой (по умолчанию от 0 до 1)
        public double GetRandom_Double(double min = 0.0, double max = 1.0)
        {
            return rnd.NextDouble() * (max - min) + min;
        }
        
        //Возвращает рандомное целое число
        public int GetRandom_Int(int min = 0, int max = 100)
        {
            return rnd.Next(min, max);
        }
        
        //Перевод десятичного числа в код Грея(строка) установленной длинны
        public string DecimalToGrey(int dec)
        {
            int result = dec ^ (dec >> 1);
            return Convert.ToString(FirstZeroBinary(Convert.ToString(dec ^ (dec >> 1), 2)));
        }
        
        //Переводит код Грея в десятичное число
        public int GreyToDecimal(string greyStr)
        {
            int grey = Convert.ToInt32(greyStr, 2);
            int mask = grey;
            while (mask != 0)
            {
                mask >>= 1;
                grey ^= mask;
            }
            return grey;
        }


        
       

        

        
        


        

        //Вывод популяции на график
        public void PopulationOnGraph()
        {

        }
        
    }
}
