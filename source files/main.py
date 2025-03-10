import pandas as pd
import matplotlib.pyplot as plt
from Matrix import Matrix
import math

def readCoordinates(name_file):
    df = pd.read_csv(name_file)
    x, y = df['Dystans (m)'].to_list(), df['Wysokość (m)'].to_list()
    return x, y

def plotHeightProfile(x, y, interpolatePointsX, interpolatePointsY, interpolateX, interpolateY, title):
    plt.plot(x, y)

    plt.scatter(interpolatePointsX, interpolatePointsY, c='red', s=100)

    plt.plot(interpolateX, interpolateY)

    plt.xlabel("Dystans [m]")
    plt.ylabel("Wysokość [m n.p.m.]")
    plt.title(title)
    plt.legend(["Oryginalna funkcja", "Węzły interpolacji", "Interpolacja funkcji"])
    plt.show()


def choosePoints(x, y, numPoints):
    diff = int(len(x) / (numPoints-1))
    interpolatePointsX, interpolatePointsY = [], []
    interpolatePointsX.append(x[0])
    interpolatePointsY.append(y[0])
    index = diff
    for i in range(numPoints-2):
        interpolatePointsX.append(x[index])
        interpolatePointsY.append(y[index])
        index += diff
    interpolatePointsX.append(x[-1])
    interpolatePointsY.append(y[-1])
    return interpolatePointsX, interpolatePointsY

def LagrangeInterpolate(point, interpolatePointsX, interpolatePointsY, x):
    interpolateX, interpolateY = [], []

    for i in range(len(x)):
        value = 0
        for j in range(point):
            denominator, numerator = 1, 1
            for k in range(point):
                if (k != j):
                    numerator = numerator * (x[i] - interpolatePointsX[k])
                    denominator = denominator * (interpolatePointsX[j] - interpolatePointsX[k])
            current_value = numerator / denominator
            value += current_value * interpolatePointsY[j]

        interpolateX.append(x[i])
        interpolateY.append(value)

    return interpolateX, interpolateY

def chooseChebyshevPoints(x, y, point):
    interpolatePointsX, interpolatePointsY = [], []

    interpolatePointsX.append(x[0])
    interpolatePointsY.append(y[0])

    for i in range(point-1, 1, -1):
        index = int((len(x)-1)/2 + (len(x)-1)/2 * math.cos((2*i-1) * math.pi / (2*point)))
        interpolatePointsX.append(x[index])
        interpolatePointsY.append(y[index])

    interpolatePointsX.append(x[-1])
    interpolatePointsY.append(y[-1])

    return interpolatePointsX, interpolatePointsY


def splajnyInterpolate(point, interpolatePointsX, interpolatePointsY, x):
    coefficients = Matrix((point-1)*4, (point-1)*4)
    values = Matrix(1, (point-1)*4)

    index, row = 0, 0
    for i in range(point-1):
        coefficients.setValueInPosition(row, index, 1)
        values.setValueInPosition(row, 0, interpolatePointsY[i])
        h = interpolatePointsX[i+1] - interpolatePointsX[i]
        row += 1
        coefficients.setValueInPosition(row, index, 1)
        coefficients.setValueInPosition(row, index+1, h)
        coefficients.setValueInPosition(row, index+2, h*h)
        coefficients.setValueInPosition(row, index+3, h*h*h)
        values.setValueInPosition(row, 0, interpolatePointsY[i+1])

        index += 4
        row += 1

    index = 0
    for i in range(point-2):
        h = interpolatePointsX[i+1] - interpolatePointsX[i]
        coefficients.setValueInPosition(row, index+1, 1)
        coefficients.setValueInPosition(row, index + 2, 2*h)
        coefficients.setValueInPosition(row, index + 3, 3*h*h)
        coefficients.setValueInPosition(row, index + 5, -1)

        row += 1

        coefficients.setValueInPosition(row, index + 2, 2)
        coefficients.setValueInPosition(row, index + 3, 6 * h)
        coefficients.setValueInPosition(row, index + 6, -2)

        row += 1
        index += 4

    coefficients.setValueInPosition(row, 2, 1)
    row += 1
    coefficients.setValueInPosition(row, index + 2, 2)
    h = interpolatePointsX[-1] - interpolatePointsX[-2]
    coefficients.setValueInPosition(row, index + 3, 6 * h)

    resultCoefficients = coefficients.faktoryzacjaLU(coefficients, values)

    index, interpolateX, interpolateY, offset = 0, [], [], 0
    for i in range(len(x)):
        interpolateX.append(x[i])
        h = x[i] - interpolatePointsX[offset]
        value = resultCoefficients.getValueInPosition(index, 0)
        value += resultCoefficients.getValueInPosition(index+1, 0) * h
        value += resultCoefficients.getValueInPosition(index+2, 0) * h * h
        value += resultCoefficients.getValueInPosition(index+3, 0) * h * h * h
        interpolateY.append(value)
        if x[i] == interpolatePointsX[offset + 1]:
            index += 4
            offset += 1

    return interpolateX, interpolateY

def doInterpolate(file_name, title):
    x, y = readCoordinates(file_name)

    points = [5, 10, 15, 20]


    # interpolacja metodą Lagrange'a dla punktów rozmieszczonych równomiernie
    for point in points:
        interpolatePointsX, interpolatePointsY = choosePoints(x, y, point)
        interpolateX, interpolateY = LagrangeInterpolate(point, interpolatePointsX, interpolatePointsY, x)
        plotHeightProfile(x, y, interpolatePointsX, interpolatePointsY, interpolateX, interpolateY, f'Lagrange {title} - {point} węzły rozmieszczone równomiernie')


    # interpolacja metodą Lagrange'a dla punktów rozmieszczonych nierównomiernie
    for point in points:
        interpolatePointsX, interpolatePointsY = chooseChebyshevPoints(x, y, point)
        interpolateX, interpolateY = LagrangeInterpolate(point, interpolatePointsX, interpolatePointsY, x)
        plotHeightProfile(x, y, interpolatePointsX, interpolatePointsY, interpolateX, interpolateY, f'Lagrange {title} - {point} węzły Czebyszewa')


    points = [5, 10, 15, 40]
    
    # interpolacja splajnami dla punktów rozmieszczonych równomiernie
    for point in points:
        interpolatePointsX, interpolatePointsY = choosePoints(x, y, point)
        interpolateX, interpolateY = splajnyInterpolate(point, interpolatePointsX, interpolatePointsY, x)
        plotHeightProfile(x, y, interpolatePointsX, interpolatePointsY, interpolateX, interpolateY, f'Splajny {title} - {point} węzły rozmieszczone równomiernie')



doInterpolate("MountEverest.csv", "Mount Everest")
doInterpolate("SpacerniakGdansk.csv", "Spacerniak Gdańsk")
doInterpolate("WielkiKanionKolorado.csv", "Wielki Kanion Kolorado")



