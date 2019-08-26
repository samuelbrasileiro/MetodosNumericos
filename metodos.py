from math import *
from cmath import *
from sympy import *
from sympy.parsing.sympy_parser import parse_expr

arq = open('entrada.txt', 'r')
entrada = arq.read()
metodos = entrada.splitlines()

def Euler(metodo):
    metodo = metodo.split(' ')
    print(metodo[0], metodo[1],metodo[2],metodo[3],metodo[4])
    #y0, t0, h, passos, f

    #DECLARACAO
    y0 = float(metodo[0])
    t0 = float(metodo[1])
    h = float(metodo[2])
    passos = int(metodo[3],10)
    funcao = parse_expr(metodo[4])
    y = symbols('y')
    t = symbols('t')

    def convert(Y,T):
        f = funcao.subs(y, Y)
        f = f.subs(t,T)
        return f
    #CALCULO DO EULER
    for i in range(passos):
        f = convert(y0,t0)
        y0 = y0 + h*f
        t0 = t0 + h
        print ('i =', i+1, '\ty =', y0)


Euler(metodos[0])