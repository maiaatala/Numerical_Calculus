from sympy import *
from sympy import init_printing
from icecream import ic
import random

init_printing(use_unicode=True)
""" x, y = symbols("x y")
expression = sin(x) * exp(x)
pprint(expression)
# pprint(integrate(diff(expression)))
currentXValue = N(diff(expression), 4, subs={x: 1})
ic(currentXValue)
pprint(N(expression, 20, subs={x: currentXValue})) """


class BasicInfo:
    eq = acceptableError = None


class NewtonRaphson(BasicInfo):
    currentXValue = diffeq = None

    def __init__(self, *args):
        self.eq, self.currentXValue, self.acceptableError = args


class Weighted(BasicInfo):
    leftXValue = rightXValue = None

    def __init__(self, *args):
        self.eq, self.leftXValue, self.rightXValue, self.acceptableError = args


class Secante(BasicInfo):
    pastXValue = currentXValue = None

    def __init__(self, *args):
        self.eq, self.pastXValue, self.currentXValue, self.acceptableError = args


def SecanteMethod(eq, pastXValue, currentXValue, accpError, iteration):
    if iteration > 100:
        print("iteration limit reached")
        return None, iteration

    if abs(currentXValue - pastXValue) <= accpError:
        return random.choice([pastXValue, currentXValue]), iteration

    if abs((pastYValue := N(eq, 20, subs={x: pastXValue}))) <= accpError:
        return pastXValue, iteration
    if abs((currentYValue := N(eq, 20, subs={x: currentXValue}))) <= accpError:
        return currentXValue, iteration

    iteration = iteration + 1
    tempXValue = currentXValue - (currentYValue / (currentYValue - pastYValue)) * (
        currentXValue - pastXValue
    )

    ic(iteration, tempXValue, N(eq, 20, subs={x: tempXValue}))

    return SecanteMethod(eq, currentXValue, tempXValue, accpError, iteration)


def NewtonMethod(eq, diffeq, currentXValue, accpError, iteration):
    YValue = N(eq, 20, subs={x: currentXValue})

    # ic(iteration, currentXValue, YValue)
    if abs(YValue) <= accpError:
        return currentXValue, iteration

    iteration = iteration + 1
    DiffYValue = N(diffeq, 20, subs={x: currentXValue})
    currentXValue = currentXValue - YValue / DiffYValue

    if iteration > 100:
        print("iteration limit reached")
        return None, iteration

    return NewtonMethod(eq, diffeq, currentXValue, accpError, iteration)


def WeightedMethod(eq, leftXValue, rightXValue, accpError, iteration):
    if iteration > 100:
        print("iteration limit reached")
        return None, iteration

    if abs(rightXValue - leftXValue) <= accpError:
        return random.choice([leftXValue, rightXValue]), iteration

    if abs((leftYValue := N(eq, 20, subs={x: leftXValue}))) <= accpError:
        return leftXValue, iteration
    if abs((rightYValue := N(eq, 20, subs={x: rightXValue}))) <= accpError:
        return rightXValue, iteration

    iteration = iteration + 1
    tempXValue = (leftXValue * rightYValue - rightXValue * leftYValue) / (
        rightYValue - leftYValue
    )
    tempYValue = N(eq, 20, subs={x: tempXValue})
    # ic(leftXValue, rightXValue, leftYValue, rightYValue, tempXValue, tempYValue)

    if leftYValue > 0:
        if tempYValue > 0:
            return WeightedMethod(eq, tempXValue, rightXValue, accpError, iteration)
        else:
            return WeightedMethod(eq, leftXValue, tempXValue, accpError, iteration)
    else:
        if tempYValue < 0:
            return WeightedMethod(eq, tempXValue, rightXValue, accpError, iteration)
        else:
            return WeightedMethod(eq, leftXValue, tempXValue, accpError, iteration)


def BisecMethod(eq, leftXValue, rightXValue, accpError, iteration):
    if iteration > 100:
        print("iteration limit reached")
        return None

    if abs(rightXValue - leftXValue) <= accpError:
        return random.choice([leftXValue, rightXValue]), iteration

    if abs((leftYValue := N(eq, 20, subs={x: leftXValue}))) <= accpError:
        return leftXValue
    if abs((rightYValue := N(eq, 20, subs={x: rightXValue}))) <= accpError:
        return rightXValue

    iteration = iteration + 1
    tempXValue = (leftXValue + rightXValue) / 2
    tempYValue = N(eq, 20, subs={x: tempXValue})
    # ic(leftXValue, rightXValue, leftYValue, rightYValue, tempXValue, tempYValue)

    if leftYValue > 0:
        if tempYValue > 0:
            return WeightedMethod(eq, tempXValue, rightXValue, accpError, iteration)
        else:
            return WeightedMethod(eq, leftXValue, tempXValue, accpError, iteration)
    else:
        if tempYValue < 0:
            return WeightedMethod(eq, tempXValue, rightXValue, accpError, iteration)
        else:
            return WeightedMethod(eq, leftXValue, tempXValue, accpError, iteration)


x = symbols("x")
method = 3
foundXValue = None
k = 0

# newton
if method == 2:
    # metodoNewton = NewtonRaphson(sin(x) + ln(x), 0.55, 1 / 100000)
    metodoNewton = NewtonRaphson(2 * pow(x, 3) - cos(x + 1) - 3, -1, 1 / 10000000)
    pprint(metodoNewton.eq)
    metodoNewton.diffeq = diff(metodoNewton.eq)

    foundXValue, k = NewtonMethod(
        metodoNewton.eq,
        metodoNewton.diffeq,
        metodoNewton.currentXValue,
        metodoNewton.acceptableError,
        k,
    )

    print(f"The x solve value is: {foundXValue}\nFound in {k} tries")

# weighted
if method == 1:
    # metodoBisec = Weighted(ln(x) + x, 0.001, 1, 1 / 100000)
    metodoWeighted = Weighted(2 * pow(x, 3) - cos(x + 1) - 3, -1, 2, 1 / 100000)
    pprint(metodoWeighted.eq)

    foundXValue, k = WeightedMethod(
        metodoWeighted.eq,
        metodoWeighted.leftXValue,
        metodoWeighted.rightXValue,
        metodoWeighted.acceptableError,
        0,
    )

    print(f"The x solve value is: {foundXValue}\nFound in {k} tries")

# bisec
if method == 0:
    # metodoBisec = Weighted(ln(x) + x, 0.001, 1, 1 / 100000)
    metodoBisec = Weighted(2 * pow(x, 3) - cos(x + 1) - 3, -1, 2, 1 / 100000)
    pprint(metodoBisec.eq)

    foundXValue, k = BisecMethod(
        metodoBisec.eq,
        metodoBisec.leftXValue,
        metodoBisec.rightXValue,
        metodoBisec.acceptableError,
        0,
    )

    print(f"The x solve value is: {foundXValue}\nFound in {k} tries")

# secante
if method == 3:
    # metodoBisec = Weighted(ln(x) + x, 0.001, 1, 1 / 100000)
    metodoSecante = Secante(2 * pow(x, 3) - cos(x + 1) - 3, 2, -1, 1 / 100000)
    pprint(metodoSecante.eq)

    foundXValue, k = SecanteMethod(
        metodoSecante.eq,
        metodoSecante.pastXValue,
        metodoSecante.currentXValue,
        metodoSecante.acceptableError,
        0,
    )

    print(f"The x solve value is: {foundXValue}\nFound in {k} tries")
