
class Matrix:

    def __init__(self, width, length):
        self._width = width
        self._length = length
        self._values = [[0 for _ in range(width)] for _ in range(length)]

    def getWidth(self):
        return self._width

    def getLength(self):
        return self._length

    def getValueInPosition(self, row, column):
        return self._values[row][column]

    def setValueInPosition(self, row, column, value):
        self._values[row][column] = value

    def printMatrix(self):
        for w in range(self._length):
            print(self._values[w])


    def _pivoting(self, L, U, index, b):
        value, indexChange = abs(U.getValueInPosition(index, index)), index

        for i in range(index+1, U.getLength()):
            if abs(U.getValueInPosition(i, index)) > value:
                value, indexChange = abs(U.getValueInPosition(i, index)), i

        if indexChange != index:
            for j in range(U.getWidth()):
                if j >= index:
                    pom = U.getValueInPosition(index, j)
                    U.setValueInPosition(index, j, U.getValueInPosition(indexChange, j))
                    U.setValueInPosition(indexChange, j, pom)
                else:
                    pom = L.getValueInPosition(index, j)
                    L.setValueInPosition(index, j, L.getValueInPosition(indexChange, j))
                    L.setValueInPosition(indexChange, j, pom)

            pom = b.getValueInPosition(index, 0)
            b.setValueInPosition(index, 0, b.getValueInPosition(indexChange, 0))
            b.setValueInPosition(indexChange, 0, pom)


    def _createLU(self, A, L, U, b):
        for j in range(L.getLength()):
            for i in range(L.getWidth()):
                if i == j:
                    L.setValueInPosition(i, j, 1)

        for j in range(A.getWidth()):

            self._pivoting(L, U, j, b)

            for i in range(j+1, A.getLength()):
                value = U.getValueInPosition(i, j) / U.getValueInPosition(j,j)
                L.setValueInPosition(i, j, value)

                for k in range(j, A.getWidth()):
                    value = U.getValueInPosition(i, k) - L.getValueInPosition(i, j) * U.getValueInPosition(j, k)
                    U.setValueInPosition(i, k, value)


    def _solveY(self, L, b):
        y = Matrix(1, b.getLength())

        for i in range(y.getLength()):
            value = b.getValueInPosition(i, 0)
            for k in range(i):
                value -= L.getValueInPosition(i, k) * y.getValueInPosition(k, 0)
            y.setValueInPosition(i, 0, value)

        return y

    def _solveX(self, U, y):
        x = Matrix(1, y.getLength())

        for i in range(x.getLength()-1, -1, -1):
            value = y.getValueInPosition(i, 0)
            for k in range(y.getLength()-1, i, -1):
                value -= U.getValueInPosition(i, k) * x.getValueInPosition(k, 0)
            x.setValueInPosition(i, 0, value/U.getValueInPosition(i,i))

        return x



    def faktoryzacjaLU(self, A, b):
        L, U = Matrix(A.getWidth(), A.getLength()), A

        self._createLU(A, L, U, b)

        y = self._solveY(L, b)

        return self._solveX(U, y)
