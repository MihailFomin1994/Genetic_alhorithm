import sys
from PyQt5 import QtWidgets
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QLabel, QGridLayout, QPushButton, QVBoxLayout
import GUI
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib import cm
import GA
from mpl_toolkits.mplot3d import Axes3D, proj3d
import matplotlib.animation as animation
import time
import random

pointCount = 500


class GeneticAlgoriphm(QtWidgets.QMainWindow, GUI.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.OnStart()
        self.RunGA.pressed.connect(self.OnGA)


    #При запуске программы
    def OnStart(self):
        #Виджет целевой функции
        self.objLayout = QtWidgets.QVBoxLayout()
        self.objWidget.setLayout(self.objLayout)
        
        #Виджет проверочной функции
        self.fitLayout = QtWidgets.QVBoxLayout()
        self.fitWidget.setLayout(self.fitLayout)
        self.CheckInputParameters()
        self.CreateGraphics()

        X = np.linspace(0, 100, pointCount)
        Y = np.random.rand(pointCount) + 50
        
        #Запосление графика приспособленности при старте программы
        self.fit_graph, = self.axes_fit.plot(X, Y, c = 'blue')

    #Функция запуска
    def OnGA(self):
        self.ClearGraphics()
        self.CheckInputParameters()
        if self.statusRun == True:
            self.CreateGraphics()
            self.StartGA()
        else:
            print("Неверные вхожные параметры!")

    #Очистка виджетов графиков
    def ClearGraphics(self):
        self.valueExtremum.setText("none")
        #Очистка виджета целевой функции
        while self.objLayout.count():
            child = self.objLayout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        #Очистка виджета функции приспособленности
        while self.fitLayout.count():
            child = self.fitLayout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
    
    #Проверка параметров
    def CheckInputParameters(self):
        #Просто проверяем на пустоту инпуты и корректность ввода в них
        self.Value_func = ""
        self.Value_fitFunc = ""
        self.Value_N = 0
        self.Value_popCount = 0
        self.Value_Pc = 0
        self.Value_Pm = 0
        self.Value_minL = 0
        self.Value_maxL = 0
        self.Value_greyCode = False
        self.Value_parentRandom = False
        self.Value_reductionRoulette = False
        self.Value_MAX = 0
        self.Value_dF = 0
        
        try:
            #region Проверка функций
            self.Plot3D = False
            global func, prFunc
            x = random.random()
            y = random.random()
            #Проверка целевой функции
            if len(self.objFunc.text()) == 0:
                GA.func = "x**2 + y**2 + 2"
            else:
                GA.func = str(self.objFunc.text())
                GA.func = GA.func.lower()
                self.Plot3D = False
                if 'y' in GA.func:
                    self.Plot3D = True
                    GA.f(x, y)
                else:
                    self.Plot3D = False
                    GA.f(x)
            self.Value_func =  GA.func
            #Проверка проверочной функции
            if len(self.fitFunc.text()) == 0:
                GA.prFunc = "f(x, y) + 3"
            else:
                GA.prFunc = str(self.fitFunc.text())
                GA.prFunc = GA.prFunc.lower()
                fitCheck3D = False
                #Совпадение количества переменных целевой и проверочной функций
                if 'y' in GA.prFunc:
                    fitCheck3D = True
                    GA.F(x, y)
                else:
                    fitCheck3D = False
                    GA.F(x)
                if self.Plot3D != fitCheck3D:
                    print("Несоотвествие целевой и проверочной функции")
                    self.statusRun = False
                else:
                    self.statusRun = True
            self.Value_fitFunc = GA.prFunc
            #endregion
            
            #region Проверка других параметров
            #Количество генов в хромосоме
            if len(self.N.text()) == 0:
                self.Value_N = 20
            else:
                self.Value_N = int(self.N.text())
                if self.Value_N<2:
                    print("Неверно указанно количество кодирования хромосомы")
                    self.statusRun = False
            
            #Количество особей в популяции
            if len(self.popCount.text()) == 0:
                self.Value_popCount = 20
            else:
                #Невозможность при заданных параметрах создать популяцию
                R = 2**self.Value_N
                self.Value_popCount = int(self.popCount.text())
                if (self.Value_popCount%2 != 0) or (self.Value_popCount < 2) or (R < self.Value_popCount):
                    print("Невозможно создать популяцию из представленных фенотипов или популяция = 0")
                    self.statusRun = False
            
            #Вероятности скрещивания и мутации
            if len(self.Pc.text()) == 0:
                self.Value_Pc = 0.85
            else:
                self.Value_Pc = float(self.Pc.text())
                if self.Value_Pc<0 or self.Value_Pc>1:
                    print("Вероятность скрещивания указана неверно")
                    self.statusRun = False
            
            if len(self.Pm.text()) == 0:
                self.Value_Pm = 0.85
            else:
                self.Value_Pm = float(self.Pm.text())
                if self.Value_Pm<0 or self.Value_Pm>1:
                    print("Вероятность мутации указана неверно")
                    self.statusRun = False
            
            #Проверка диапазона поиска
            if len(self.minL.text()) == 0:
                self.Value_minL = 0
            else:
                self.Value_minL = float(self.minL.text())
            
            if len(self.maxL.text()) == 0:
                self.Value_maxL = 0
            else:
                self.Value_maxL = float(self.maxL.text())

            if self.Value_maxL <= self.Value_minL:
                print("Неправильно указанны диапазоны поиска")
                self.statusRun = False
            
            #Проверка радиобоксов на указанное значение
            self.Value_greyCode = False
            self.Value_parentRandom = True
            self.Value_reductionRoulette = True
            
            if self.radio_binCode.isChecked() == False:
                self.Value_greyCode = True
            
            if self.radio_parentRandom.isChecked() == False:
                self.Value_parentRandom = False
            
            if self.radio_roulete.isChecked() == False:
                self.Value_reductionRoulette = False
            #Проверка количества итераций
            if len(self.MAX.text()) == 0:
                self.Value_MAX = 100
            else:
                self.Value_MAX = int(self.MAX.text())
                if self.Value_MAX < 1:
                    print("Неверное количество итераций")
                    self.statusRun = False
            
            #Приращение приспособленности
            self.Value_dF = float(self.dF.text())

            
            #endregion
            
        except inputErrore:
            self.statusRun = False
            
            

    #Создание виджетов графиков
    def CreateGraphics(self):
        #Виджет целевой функции
        self.obj_fig = Figure()
        self.obj_canvas = FigureCanvas(self.obj_fig)
        self.objLayout.addWidget(self.obj_canvas)
        
        #Виджет проверочной функции
        self.fit_fig = Figure()
        self.fit_canvas = FigureCanvas(self.fit_fig)
        
        self.fitLayout.addWidget(self.fit_canvas)
        self.axes_fit = self.fit_canvas.figure.subplots()

        self.inctimentData = list()
        self.fitFuncData = list()

        self.fit_graph, = self.axes_fit.plot([], [], c = 'blue')

        if self.Plot3D == True:
            self.CreateFuncGaraph_3D()
        else:
            self.CreateFuncGaraph_2D()
    
    #Создание графика поверхности целевой функции
    def CreateFuncGaraph_3D(self):
        #График целевой функции
        self.axes_obj = self.obj_fig.add_subplot(projection='3d')
        #Используется точность
        X = np.linspace(float(self.Value_minL), float(self.Value_maxL), int(pointCount))
        Y = np.linspace(float(self.Value_minL), float(self.Value_maxL), int(pointCount))
        X, Y = np.meshgrid(X, Y)
        Z = GA.f(X, Y)
        self.axes_obj.plot_surface(X, Y, Z, cmap = cm.coolwarm, alpha=0.7, antialiased=True)
        self.points = self.axes_obj.scatter([], [], [], c = 'black')
            
    #Создание графика целевой функции
    def CreateFuncGaraph_2D(self):
        self.axes_obj = self.obj_fig.subplots()
        X = np.linspace(float(self.Value_minL), float(self.Value_maxL), int(pointCount))
        Y = GA.f(X)
        self.axes_obj.plot(X, Y, c='orange')

        self.axes_obj.set_xlim(self.Value_minL, self.Value_maxL)

        #self.axes_obj.set_ylim(self.minL - 1, self.maxL + 1)
        self.points = self.axes_obj.scatter([], [], c = 'black')

    #Обновление особей популяции 3D
    def updatePlot_3D(self, X, Y, Z):
        self.points._offsets3d = (X, Y, Z)
        self.points.figure.canvas.draw() 

    #Обновление особей популяции 2D
    def updatePlot_2D(self, X, Y):
        self.points.set_offsets(np.c_[X, Y])
        self.points.figure.canvas.draw()
    
    #Обновление функции приспособленности
    def updateFitPlot(self):
        self.inctimentData.append(self.ga.countP)
        self.fitFuncData.append(self.ga.Fs())
        X = np.array(self.inctimentData)
        Y = np.array(self.fitFuncData)
        minLX = np.min(X)
        maxLX = np.max(X)
        minLY = np.min(Y)
        maxLY = np.max(Y)
        if (minLX != maxLX):
            self.axes_fit.set_xlim(minLX, maxLX)
            self.axes_fit.set_ylim(minLY, maxLY)
        self.fit_graph.set_xdata(X)
        self.fit_graph.set_ydata(Y)
        self.fit_graph.figure.canvas.draw()
    
    #Обновление на графике популяции
    def UpdatePopulation(self):
        coordinates = self.ga.GetPopulationToGraph()
        #Обновление популяции на графике целевой функции
        if self.Plot3D == True:
            self.updatePlot_3D(coordinates[0], coordinates[1], coordinates[2])
        else:
            self.updatePlot_2D(coordinates[0], coordinates[1])
        #Обновление функции приспособленности
        self.updateFitPlot()
        #Обновление счетчика номера популяции
        self.countPopulation.setText(str(self.ga.countP))

    #Начало работы алгоритма
    def StartGA(self):
        #Создание экземпляра генетического алгоритма
        self.ga = GA.GA(self.Value_func, self.Value_fitFunc, self.Value_N, self.Value_popCount, 
        self.Value_Pc, self.Value_Pm, self.Value_dF, self.Value_MAX, self.Value_minL, 
        self.Value_maxL, self.Value_greyCode, self.Value_parentRandom, self.Value_reductionRoulette)
        #Начальная популяция
        self.ga.CreateStartPopulation()
        #Вывод начальной популяции на экран
        self.UpdatePopulation()
        print("Population №0")
        self.OnWindow()

        self.timer = self.fit_canvas.new_timer(50)
        self.timer.add_callback(self.Run)
        self.timer.start()
            
    #Запуск поиска экстремума
    def Run(self):
        Fsn = self.ga.Fs()
        dFn = (Fsn - self.ga.Fsn_1)/Fsn
        if (self.ga.countP >= self.Value_MAX) or (dFn < self.Value_dF):
            self.timer.stop()
            coordinates = self.ga.ReturnMaximumFitness()
            if self.Plot3D == False:
                self.valueExtremum.setText("(" + str(round(coordinates[0], 3)) + "; " + str(round(coordinates[1], 3)) + " )")
            else:
                self.valueExtremum.setText("(" + str(round(coordinates[0], 3)) + "; " + str(round(coordinates[1], 3)) + "; " + str(round(coordinates[2], 3)) + " )")
            #Вывод последней популяции на экран
            print("Population №" + str(self.ga.countP))
            self.OnWindow()
        else:
            self.ga.Crossbreing()
            self.ga.Mutatin()
            self.ga.Selection()
            self.UpdatePopulation()
            self.ga.Fsn_1 = Fsn
            #Вывод популяции №3 на экран
            if(self.ga.countP == 2):
                print("Population №" + str(self.ga.countP))
                self.OnWindow()

    #Вывод на экран
    def OnWindow(self):
        repotr = self.ga.GetPopulationToReport()
        for el in repotr:
            print(el + '\n')


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = GeneticAlgoriphm()
    window.show()
    app.exec_()

if __name__ == '__main__':
    main()