import time
import sympy.combinatorics.graycode as b_g
import numpy as np
import random
from numpy import sign, sin, cos, sqrt, exp

func = "x"
prFunc = "f(x)"



def f(x, y = 0):
    global func
    return eval(func)
    
def F(x, y = 0):
    global prFunc
    return eval(prFunc)


#Генератор случайных чисел
def GetRandom_Double():
    return random.random()

def GetRandom_Int(max):
    return random.randint(0, max)

#Возвращает list неповторяющихся случайных номеров
def ReturnRandomAr(countEl, maxValue):
    randomAr = list()
    sizeAr = countEl - 1
    if(countEl == maxValue):
        lenAr = [i for i in range(0, countEl)]
        while sizeAr >= 0:
            i = GetRandom_Int(sizeAr)
            randomAr.append(lenAr[i])
            del lenAr[i]
            sizeAr = sizeAr - 1
    
    else:
        value = GetRandom_Int(maxValue)
        
        while(sizeAr >= 0):
            if(value not in randomAr):
                randomAr.append(value)
                sizeAr = sizeAr - 1
            value = GetRandom_Int(maxValue)
    return randomAr


class GA(object):
    #Конструктор
    def __init__(self, f = "x**2", F = "f(x)+5", N = 20, popCount = 50, Pc = 0.85, Pm = 0.15, dF = 0.01, 
    MAX_ITERATION = 100, MIN_lim = -1.0, MAX_lim = 1.0, codeGrey = False, parenRandom = False, reductionRoulete = False):
        try:
            #Инициализация устанавливаемых параметров
            global func, prFunc
            func = str(f)                           #Целевая функция
            prFunc = str(F)                         #Проверочная функция
            
            self.N = int(N)                         #Количество бит кодирования хромосомы
            self.popCount = int(popCount)           #Особей в популяции (кратно 2)
            self.Pc = float(Pc)                     #Вероятность скрещивания
            self.Pm = float(Pm)                     #Вероятность мутации
            self.dF = float(dF)                     #Минимальное приращение приспособленности
            self.MAX_ITERATION = int(MAX_ITERATION) #Максимальное количество популяций
            self.MIN = float(MIN_lim)               #Минимум диапазона поиска экстремума
            self.MAX = float(MAX_lim)               #Максимальный диапазон поиска экстремума
            self.grey = codeGrey                    #Да - код Грея, нет - бинарный код
            self.parenRandom = parenRandom          #Да - случайный отбор, нет - лучшие с лучшими
            self.reductionRoulete = reductionRoulete#Да - редукция Рулетка, нет - турнир
            
            #Инициализация используемых параметров
            self.population = list()
            self.numIter = 0
            self.R = 2**N
            self.step = (self.MAX-self.MIN)/(self.R - 1)
            self.twoVariable = False
            self.Fsn_1 = -10000
            self.countP = 0
        except:
            print("Неверный формат ввода данных")
        self.CheckInputParametres()
    
    #region Функции
    #Возвращает значение переменной по хромосоме
    def Value(self, chromosome):
        return self.MIN + float(self.Decoding(chromosome) * self.step)
    #Возвращает результат приспособленности
    def Fx(self, chromosome):
        if(self.twoVariable == False):
            x = self.MIN + self.Decoding(chromosome) * self.step
            return F(x)
        else:
            x = self.MIN + self.Decoding(chromosome[0]) * self.step
            y = self.MIN + self.Decoding(chromosome[1]) * self.step
            return F(x, y)
    #Возвращает значение целевой функции по хромосоме
    def fxy(self, chromosome):
        if(self.twoVariable == False):
            x = self.Value(chromosome)
            return f(x)
        else:
            x = self.Value(chromosome[0])
            y = self.Value(chromosome[1])
            return f(x, y)
   
    #Возвращает координаты особи с максимальной приспособленностью
    def ReturnMaximumFitness(self):
        numIndiv = len(self.population)
        populationFx = list()
        maxChrom = ""
        max_Fxy = -100000
        for chrom in self.population:
            Fxy = self.Fx(chrom)
            if Fxy > max_Fxy:
                maxChrom = chrom
                max_Fxy = Fxy
        coordinate = []
        if self.twoVariable == False:
            x = self.Value(maxChrom)
            fx = f(x)
            coordinate = [x, fx]
        else:
            x = self.Value(maxChrom[0])
            y = self.Value(maxChrom[1])
            fxy = f(x, y)
            coordinate = [x, y, fxy]
        return coordinate
    #endregion

    #region Проверка входных параметров
    def CheckInputParametres(self):
        #Функция одной или 2 переменных
        global func, prFunc
        if 'y' in func:
            self.twoVariable = True
        else:
            self.twoVariable = False
        check_prFunk = False
        if 'y' in prFunc:
            check_prFunk = True
        else:
            check_prFunk = False
        #Управление вводом отдано управляющему классу!!!
        """
        if self.twoVariable != check_prFunk: 
            print("Несоотвествие целевой и проверочной функции")
            exit()
        #Четное ли количество особей в популяции
        if self.N%2 != 0:
            print("Нечетное количество особей в популяции")
            exit()
        #Невозмодность создать установленное количество особей в популяции, используя данную кодировку
        if (self.R < self.popCount) or (self.popCount < 2):
            print("Невозможно создать популяцию из представленных фенотипов или популяция = 0")
            exit()
        #Проверка интервалов
        if self.MIN >= self.MAX:
            print("Неправильный ввод интервалов")
            exit()
        """

    #endregion
    
    #region Кодирование/Декодирование
    #Возвращает хромосому но номеру фенотипа
    def Coding(self, decimal):
        if(self.grey == False):
            chromosome = self.DecimalToBinary(decimal)
        else:
            chromosome = self.DecimalToGrey(decimal)
        return self.ChromosomeToLen(chromosome)
    #Возвращает номер фенотипа по хромосоме
    def Decoding(self, chromosome):
        if(self.grey == False):
            decimal = self.BinaryToDecimal(chromosome)
        else:
            decimal = self.GreyToDecimal(chromosome)
        return decimal
    #endregion

    #region Способы кодирования и декодирования
    #Приведение хромосомы к установленному размеру
    def ChromosomeToLen(self, chromosome):
        return chromosome.rjust(self.N, '0')

    #Из номера фенотипа в двоичную хромосому
    def DecimalToBinary(self, decimal):
        binary = bin(decimal).replace("0b", "")
        return binary
    #Из двоичной хромосомы к номеру фенотипа
    def BinaryToDecimal(self, chromosome):
        decimal = int(chromosome, 2)
        return decimal

    #Из номера фенотипа в хромосому Грея
    def DecimalToGrey(self, decimal):
        binary = self.DecimalToBinary(decimal)
        grey = b_g.bin_to_gray(binary)
        return grey
    #Из хромосомы Грея к номеру фенотипа
    def GreyToDecimal(self, chromosome):
        binary = b_g.gray_to_bin(chromosome)
        decimal = self.BinaryToDecimal(binary)
        return decimal
    #endregion

    #region Формирование начальной популяции
    def CreateStartPopulation(self):
        startPopulationX = ReturnRandomAr(self.popCount, self.R)
        if(self.twoVariable == False):
            for number in startPopulationX:
                self.population.append(self.Coding(number))
        else:
            startPopulationY = ReturnRandomAr(self.popCount, self.R)
            for i in range(self.popCount):
                self.population.append([self.Coding(startPopulationX[i]), self.Coding(startPopulationY[i])])
    #endregion

    #region Скрещивание
    #Случайное скрещивание
    def Crossbreing_random(self):
        
        parentPairs_index = ReturnRandomAr(self.popCount, self.popCount)
        parentPairs = list()
        for i in range(self.popCount):
            parentPairs.append(self.population[parentPairs_index[i]])
        return parentPairs
    #Скрещивание лучшие с лучшими
    def Crossbreing_winTowin(self):
        parentPairs = list()
        for i in range(self.popCount):
            Fxi = self.Fx(self.population[i])
            parentPairs.append([Fxi, self.population[i]])
        parentPairs.sort(key=lambda i:i[0], reverse=True)
        parentPairs = [el[1] for el in parentPairs]
        return parentPairs
    #Скрещивание
    def Crossbreing(self):
        parentPairs = list()
        if(self.parenRandom == True): parentPairs = self.Crossbreing_random()
        else: parentPairs = self.Crossbreing_winTowin()

        parent_mathers = parentPairs[::2]
        parent_fathers = parentPairs[1::2]

        self.population.clear()
        for i in range(int(self.popCount/2)):
            P = GetRandom_Double()
            if(P<self.Pc):
                if self.twoVariable == False:
                    crossPoint = GetRandom_Int(self.N)
                    self.population.append(parent_mathers[i][:crossPoint:] + parent_fathers[i][crossPoint::])
                    self.population.append(parent_fathers[i][:crossPoint:] + parent_mathers[i][crossPoint::])
                else:
                    crossPointX = GetRandom_Int(self.N)
                    crossPointY = GetRandom_Int(self.N)
                    self.population.append([parent_mathers[i][0][:crossPointX:] + parent_fathers[i][0][crossPointX::], parent_mathers[i][0][:crossPointY:] + parent_fathers[i][0][crossPointY::]])
                    self.population.append([parent_mathers[i][0][crossPointX::] + parent_fathers[i][0][:crossPointX:], parent_mathers[i][0][crossPointY::] + parent_fathers[i][0][:crossPointY:]])
            else:
                self.population.append(parent_mathers[i])
                self.population.append(parent_fathers[i])
    #endregion

    #region Мутация
    #Точечная мутация
    def Mutation_point(self, chromosome):
        alleleX = GetRandom_Int(self.N - 1)
        if(self.twoVariable == False):
            if chromosome[alleleX] == '0':
                chromosome = chromosome[:alleleX] + '1' + chromosome[alleleX+1:]
            else:
                chromosome = chromosome[:alleleX] + '0' + chromosome[alleleX+1:]
        else:
            alleleY = GetRandom_Int(self.N - 1)
            if(chromosome[0][alleleX] == '0'):
                chromosome[0] = chromosome[0][:alleleX] + '1' + chromosome[0][alleleX+1:]
            else:
                chromosome[0] = chromosome[0][:alleleX] + '0' + chromosome[0][alleleX+1:]
            
            if(chromosome[1][alleleY] == '0'):
                chromosome[1] = chromosome[1][:alleleY] + '1' + chromosome[1][alleleY+1:]
            else:
                chromosome[1] = chromosome[1][:alleleY] + '0' + chromosome[1][alleleY+1:]
        self.population.append(chromosome)

    #Мутация - инверсия
    def Mutation_invers(self, chromosome):
        allelesX = ReturnRandomAr(2, self.N-1)
        allelesX.sort()
        if(self.twoVariable == False):
            chromosome = chromosome[:allelesX[0]] + chromosome[allelesX[1]] + chromosome[allelesX[0] + 1: allelesX[1]] + chromosome[allelesX[0]] + chromosome[allelesX[1] + 1:]
        else:
            allelesY = ReturnRandomAr(2, self.N-1)
            allelesY.sort()
            chromosome[0] = chromosome[0][:allelesX[0]] + chromosome[0][allelesX[1]] + chromosome[0][allelesX[0] + 1: allelesX[1]] + chromosome[0][allelesX[0]] + chromosome[0][allelesX[1] + 1:]
            chromosome[1] = chromosome[1][:allelesY[0]] + chromosome[1][allelesY[1]] + chromosome[1][allelesY[0] + 1: allelesY[1]] + chromosome[1][allelesY[0]] + chromosome[1][allelesY[1] + 1:]
        self.population.append(chromosome)
    
    #Мутация
    def Mutatin(self):
        newPopulation = self.population
        for chromosome in newPopulation:
            P = GetRandom_Double()
            if(P<self.Pm):
                if self.grey == True: self.Mutation_invers(chromosome)
                else: self.Mutation_point(chromosome)

    
    #endregion

    #region Отбор
    #Суммарная приспособленность
    def Fs(self):
        result = 0
        for gen in self.population:
            result = result + self.Fx(gen)
        return result

    
           
    
    #Возвращает list для рулетки
    def GetList_roulet(self):
        value = 0
        Fs = self.Fs()
        roulet = list()
        count_el = len(self.population)
        for i in range(count_el):
            value = value + self.Fx(self.population[i])/Fs
            roulet.append(value)
        return roulet

    #Отбор - рулетка
    def Selection_roulet(self):
        if len(self.population) == self.popCount: return
        newPopulation = list()
        max_population = len(self.population)
        stop = 0
        while stop < self.popCount:
            roulet = self.GetList_roulet()
            value_r = GetRandom_Double()
            for i in range(max_population):
                if value_r <= roulet[i]:
                    newPopulation.append(self.population.pop(i))
                    
                    break
            stop = stop + 1
            max_population = max_population - 1
        self.population = newPopulation
    
    #Отбор - турнир
    def Selection_tourney(self):
        if len(self.population) == self.popCount: return
        newPopulation = list()
        stop = 0
        while stop < self.popCount:
            partis = ReturnRandomAr(2, len(self.population) - 1)
            party_1 = self.Fx(self.population[partis[0]])
            party_2 = self.Fx(self.population[partis[1]])
            if party_1 > party_2:
                newPopulation.append(self.population[partis[0]])
                del self.population[partis[1]]
            else:
                newPopulation.append(self.population[partis[0]])
                del self.population[partis[1]]
            stop = stop + 1
        self.population = newPopulation
        
    #Отбор
    def Selection(self):
        if self.reductionRoulete == True: self.Selection_roulet()
        else: self.Selection_tourney()
        self.countP = self.countP + 1
    #endregion
    
    #Возвращает популяцию в виде набора координат
    def GetPopulationToGraph(self):
        coordinates = list()
        if self.twoVariable == False:
            coordinatesX = list()
            coordinatesY = list()
            for chrom in self.population:
                x = self.Value(chrom)
                y = f(x)
                coordinatesX.append(x)
                coordinatesY.append(y)
            coordinates.append(coordinatesX)
            coordinates.append(coordinatesY)
        else:
            coordinatesX = list()
            coordinatesY = list()
            coordinatesZ = list()
            for chrom in self.population:
                x = self.Value(chrom[0])
                y = self.Value(chrom[1])
                z = f(x, y)
                coordinatesX.append(x)
                coordinatesY.append(y)
                coordinatesZ.append(z)
            coordinates.append(coordinatesX)
            coordinates.append(coordinatesY)
            coordinates.append(coordinatesZ)
        coordinates = np.array(coordinates)
        return coordinates

    #Возвращает популяцию по шаблону отчета
    #| Номер особи в популяции | x | (y) | хромосома | f(x, y) | F(x, y) |
    def GetPopulationToReport(self):
        report = list()
        for i in range(len(self.population)):
            if self.twoVariable == False:
                report.append(str(i) + " | " + str(round(self.Value(self.population[i]), 3)) + " | " +  self.population[i] + " | " + str(round(self.fxy(self.population[i]), 3)) + " | " + str(round(self.Fx(self.population[i]), 3)))
            else:
                report.append(str(i) + " | " + str(round(self.Value(self.population[i][0]), 3)) + " | " +  str(round(self.Value(self.population[i][1]), 3)) + " | " +  
                self.population[i][0] + " " + self.population[i][1] + " | " + str(round(self.fxy(self.population[i]), 3)) + " | " + str(round(self.Fx(self.population[i]), 3)))
        return report

    #Тестовая функция
    def Run(self):
        self.CreateStartPopulation()
        while self.countP < self.MAX_ITERATION:
            self.Crossbreing()
            self.Mutatin()
            self.Selection()
            self.Graph()
            Fsn = self.Fs()
            dFn = (Fsn - self.Fsn_1)/Fsn
            if dFn < self.dF:
                break
            self.Fsn_1 = Fsn
        #print("--- %s seconds ---" % (time.time() - self.start))
        print("Finish on iteration = " + str(self.countP))

    




