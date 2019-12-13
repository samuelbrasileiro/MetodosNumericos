from math import *
from cmath import *
from re import *
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import matplotlib.pyplot as plt

    
def Euler(metodo):
    
    #y0, t0, h, passos, f

    #DECLARACAO
    y0 = float(metodo[0])
    y = [y0]
    t0 = float(metodo[1])
    h = float(metodo[2])
    passos = int(metodo[3],10)
    funcao = parse_expr(metodo[4])
    y_ = symbols('y')
    t_ = symbols('t')

    def convert(Y,T):
        f = funcao.subs(y_, Y)
        f = f.subs(t_,T)
        return f
    #CALCULO DO EULER
    for i in range(passos):
        f = convert(y0,t0)
        y0 = y0 + h*f
        t0 = t0 + h
        
        y.append(y0)
    
    return y

def Euler_Inverso(metodo):

    #DECLARACAO
    y0 = float(metodo[0])
    y = [y0]
    t0 = float(metodo[1])
    h = float(metodo[2])
    passos = int(metodo[3],10)
    funcao = parse_expr(metodo[4])
    y_ = symbols('y')
    t_ = symbols('t')
    
    def convert(Y,T):
        f = funcao.subs(y_, Y)
        f = f.subs(t_,T)
        return f
    #CALCULO DO EULER
    for i in range(passos):
        f = convert(y0,t0)
        f2 = convert(y0+h*f,t0+h)

        y0 = y0 + h*f2
        t0 = t0 + h

        y.append(y0)

    return y

def Euler_Aprimorado(metodo):

    #DECLARACAO
    y0 = float(metodo[0])
    y = [y0]
    t0 = float(metodo[1])
    h = float(metodo[2])
    passos = int(metodo[3],10)
    funcao = parse_expr(metodo[4])
    y_ = symbols('y')
    t_ = symbols('t')
    
    def convert(Y,T):
        f = funcao.subs(y_, Y)
        f = f.subs(t_,T)
        return f
    #CALCULO DO EULER
    for i in range(passos):
        f1 = convert(y0,t0)
        f2 = convert(y0+h*f1,t0+h)

        y0 = y0 + h*(f1+f2)/2
        t0 = t0 + h
        
        y.append(y0)

    return y

def Runge_Kutta(metodo):
    
    #DECLARACAO
    y0 = float(metodo[0])
    y = [y0]
    t0 = float(metodo[1])
    h = float(metodo[2])
    passos = int(metodo[3],10)
    funcao = parse_expr(metodo[4])
    y_ = symbols('y')
    t_ = symbols('t')
    
    def convert(Y,T):
        f = funcao.subs(y_, Y)
        f = f.subs(t_,T)
        return f
    #CALCULO DO EULER
    for i in range(passos):
        k1 = convert(y0,t0)
        k2 = convert(y0+h*k1/2,t0+h/2)
        k3 = convert(y0+h*k2/2,t0+h/2)
        k4 = convert(y0+h*k3,t0+h)

        y0 = y0 + h*(k1+2*k2+2*k3+k4)/6
        t0 = t0 + h
        y.append(y0)
        
    return y

def Adam_Bashforth(nome, metodo):
    
    #DECLARACAO
    
    ordem = int(metodo[-1],10)
    
    funcao = parse_expr(metodo[-2])
    y_ = symbols('y')
    t_ = symbols('t')
    
    passos = int(metodo[-3],10)
    h = float(metodo[-4])
    t0 = float(metodo[-5])
    numPassos = str(ordem-1)
    if(search('euler$', nome)):
        saida.write('Método Adam-Bashforth por Euler')
        y0 = Euler([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('euler_inverso$', nome)):
        saida.write('Método Adam-Bashforth por Euler Inverso')
        y0 = Euler_Inverso([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('euler_aprimorado$', nome)):
        saida.write('Método Adam-Bashforth por Euler Aprimorado')
        y0 = Euler_Aprimorado([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('runge_kutta$', nome)):
        saida.write('Método Adam-Bashforth por Runge Kutta')
        y0 = Runge_Kutta([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    else: #Com os pontos determinados
        if(search('adam_bashforth$',nome)):
            saida.write('Método Adam-Bashforth')
        y0 = metodo[:-5]
        
    for i in range(len(y0)):
        y0[i-1] = float(y0[i-1])
    
    def convert(Y,T):
        f = funcao.subs(y_, Y)
        f = f.subs(t_,T)
        return f
    #CALCULO DO EULER
    
    t0 = t0 + (ordem - 1)*h
    for i in range( passos - ordem + 1):
        if ordem == 2:
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            y0.append(y0[-1] + h*(f0*3/2 - f1*1/2))
            
        elif ordem == 3:
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            y0.append(y0[-1] + h*(f0*23/12 - f1*4/3 + f2*5/12))
            
        elif ordem == 4:
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            y0.append(y0[-1] + h*(f0*55/24 - f1*59/24 + f2*37/24 - f3*3/8))
            
        elif ordem == 5:
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            f4 = convert(y0[-5],t0-4*h)
            y0.append(y0[-1] + h*(f0*1901/720 - f1*1387/360 + f2*109/30 - f3*637/360 + f4*251/720))
            
        elif ordem == 6:
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            f4 = convert(y0[-5],t0-4*h)
            f5 = convert(y0[-6],t0-5*h)
            y0.append(y0[-1] + h*(f0*4277/1440 - f1*2641/480 + f2*4991/720 - f3*3649/720 + f4*959/480 - f5*95/288))
            
        elif ordem == 7:
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            f4 = convert(y0[-5],t0-4*h)
            f5 = convert(y0[-6],t0-5*h)
            f6 = convert(y0[-7],t0-6*h)
            y0.append(y0[-1] + h*(f0*198721/60480 - f1*18637/2520 + f2*235183/20160 - f3*10754/945 + f4*135713/20160 - f5*5603/2520 + f6*19087/60480))
        
        elif ordem == 8:
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            f4 = convert(y0[-5],t0-4*h)
            f5 = convert(y0[-6],t0-5*h)
            f6 = convert(y0[-7],t0-6*h)
            f7 = convert(y0[-8],t0-7*h)
            y0.append(y0[-1] + h*(f0*16083/4480 - f1*1152169/120960 + f2*242653/13440 - f3*296053/13440 + f4*2102243/120960 - f5*115747/13440 + f6*32863/13440 - f7*5257/17280))
        
        
        t0 = t0 + h
    
    return y0


def Adam_Moulton(nome, metodo):
    
    #DECLARACAO
    
    ordem = int(metodo[-1],10)
    
    funcao = parse_expr(metodo[-2])
    y_ = symbols('y')
    t_ = symbols('t')
    
    passos = int(metodo[-3],10)
    h = float(metodo[-4])
    t0 = float(metodo[-5])
    numPassos = str(ordem-2)
    if(search('euler$', nome)):
        saida.write('Método Adam-Moulton por Euler')
        y0 = Euler([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('euler_inverso$', nome)):
        saida.write('Método Adam-Moulton por Euler Inverso')
        y0 = Euler_Inverso([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('euler_aprimorado$', nome)):
        saida.write('Método Adam-Moulton por Euler Aprimorado')
        y0 = Euler_Aprimorado([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('runge_kutta$', nome)):
        saida.write('Método Adam-Moulton por Runge Kutta')
        y0 = Runge_Kutta([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    else: #Com os pontos determinados
        saida.write('Método Adam-Moulton')
        y0 = metodo[:-5]
    for i in range(len(y0)):
        y0[i-1] = float(y0[i-1])

    def convert(Y,T):
        f = funcao.subs(y_, Y)
        f = f.subs(t_,T)
        return f
    #CALCULO    
    t0 = t0 + (len(y0)-1)*h 
    for i in range( passos - ordem + 2):
        if ordem == 2:

            bashforth = Euler([str(i) for i in [y0[-1], t0, h, 1, metodo[-2]]]) #euler = bashforth 1
            fp = convert(bashforth[-1],t0+h)
            f0 = convert(y0[-1],t0)

            y0.append(y0[-1] + h*(fp*1/2 + f0*1/2))
            
        elif ordem == 3:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-2], y0[-1], t0-h, h, ordem - 1, metodo[-2], ordem - 1]])
            fp = convert(bashforth[-1],t0+h)
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            y0.append(y0[-1] + h*(fp*5/12 + f0*2/3 - f1*1/12))
            
        elif ordem == 4:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-3], y0[-2], y0[-1], t0-2*h, h, ordem - 1, metodo[-2], ordem - 1]])
            fp = convert(bashforth[-1],t0+h)
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            y0.append(y0[-1] + h*(fp*3/8 + f0*19/24 - f1*5/24 + f2*1/24))
            
        elif ordem == 5:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-4], y0[-3], y0[-2], y0[-1], t0 - 3*h, h, ordem - 1, metodo[-2], ordem - 1]])
            fp = convert(bashforth[-1],t0+h)
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            y0.append(y0[-1] + h*(fp*251/720 + f0*323/720 - f1*11/30 + f2*53/360 - f3*19/720))
            
        elif ordem == 6:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-5], y0[-4], y0[-3], y0[-2], y0[-1], t0 - 4*h, h, ordem - 1, metodo[-2], ordem - 1]])
            fp = convert(bashforth[-1],t0+h)
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            f4 = convert(y0[-5],t0-4*h)
            y0.append(y0[-1] + h*(fp*95/288 + f0*1427/1440 - f1*133/240 + f2*241/720 - f3*173/1440 + f4*3/160))
            
        elif ordem == 7:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-6], y0[-5], y0[-4], y0[-3], y0[-2], y0[-1], t0 - 5*h, h, ordem - 1, metodo[-2], ordem - 1]])            
            fp = convert(bashforth[-1],t0+h)
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            f4 = convert(y0[-5],t0-4*h)
            f5 = convert(y0[-6],t0-5*h)
            y0.append(y0[-1] + h*(fp*19087/60480 + f0*2713/2520 - f1*15487/20160 + f2*586/945 - f3*6737/20160 + f4*263/2520 - f5*863/60480))
        
        elif ordem == 8:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-7], y0[-6], y0[-5], y0[-4], y0[-3], y0[-2], y0[-1], t0 - 6*h, h, ordem - 1, metodo[-2], ordem - 1]])            
            fp = convert(bashforth[-1],t0+h)
            f0 = convert(y0[-1],t0)
            f1 = convert(y0[-2],t0-h)
            f2 = convert(y0[-3],t0-2*h)
            f3 = convert(y0[-4],t0-3*h)
            f4 = convert(y0[-5],t0-4*h)
            f5 = convert(y0[-6],t0-5*h)
            f6 = convert(y0[-7],t0-6*h)
            y0.append(y0[-1] + h*(fp*5257/17280 + f0*139849/120960 - f1*4511/4480 + f2*123133/120960 - f3*88547/120960 + f4*1537/4480 - f5*11351/120960 + f6*275/24192))
        
        t0 = t0 + h

    return y0


def Formula_Inversa(nome, metodo):
    
    #DECLARACAO
    
    ordem = int(metodo[-1],10)
    
    funcao = parse_expr(metodo[-2])
    y_ = symbols('y')
    t_ = symbols('t')
    
    passos = int(metodo[-3],10)
    h = float(metodo[-4])
    t0 = float(metodo[-5])
    numPassos = str(ordem-1)
    if(search('euler$', nome)):
        saida.write('Método Fórmula Inversa de Diferenciação por Euler')
        y0 = Euler([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('euler_inverso$', nome)):
        saida.write('Método Fórmula Inversa de Diferenciação por Euler Inverso')
        y0 = Euler_Inverso([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('euler_aprimorado$', nome)):
        saida.write('Método Fórmula Inversa de Diferenciação por Euler Aprimorado')
        y0 = Euler_Aprimorado([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    elif(search('runge_kutta$', nome)):
        saida.write('Método Fórmula Inversa de Diferenciação por Runge Kutta')
        y0 = Runge_Kutta([metodo[0],metodo[1], metodo[2], numPassos, metodo[-2]])
    else: #Com os pontos determinados
        saida.write('Método Fórmula Inversa de Diferenciação')
        y0 = metodo[:-5]
    for i in range(len(y0)):
        y0[i-1] = float(y0[i-1])

    def convert(Y,T):
        f = funcao.subs(y_, Y)
        f = f.subs(t_,T)
        return f
    #CALCULO    
    t0 = t0 + (len(y0)-1)*h 
    for i in range( passos - ordem + 1):
        if ordem == 2:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-2], y0[-1], t0-h, h, ordem, metodo[-2], ordem]])
            fp = convert(bashforth[-1],t0+h)

            y0.append(h*fp*3/2 + y0[-1]*4/3 - y0[-2]*1/3)
            
        elif ordem == 3:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-3], y0[-2], y0[-1], t0-2*h, h, ordem - 1, metodo[-2], ordem]])
            fp = convert(bashforth[-1],t0+h)

            y0.append(h*fp*6/11 + y0[-1]*18/11 - y0[-2]*9/11 + y0[-3]*2/11)
            
        elif ordem == 4:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-4], y0[-3], y0[-2], y0[-1], t0-3*h, h, ordem, metodo[-2], ordem]])
            fp = convert(bashforth[-1],t0+h)

            y0.append(h*fp*12/25 + y0[-1]*48/25 - y0[-2]*36/25 + y0[-3]*16/25 - y0[-4]*3/25)
            
        elif ordem == 5:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-5], y0[-4], y0[-3], y0[-2], y0[-1], t0-4*h, h, ordem, metodo[-2], ordem]])
            fp = convert(bashforth[-1],t0+h)
            
            y0.append(h*fp*60/137 + y0[-1]*300/137 - y0[-2]*300/137 + y0[-3]*200/137 - y0[-4]*75/137 + y0[-5]*12/137)
            
        elif ordem == 6:
            bashforth = Adam_Bashforth('',[str(i) for i in [y0[-6], y0[-5], y0[-4], y0[-3], y0[-2], y0[-1], t0-5*h, h, ordem, metodo[-2], ordem]])
            fp = convert(bashforth[-1],t0+h)
            
            y0.append(h*fp*60/147 + y0[-1]*360/147 - y0[-2]*450/147 + y0[-3]*400/147 - y0[-4]*225/147 + y0[-5]*72/147 - y0[-6]*10/147)
            
        
        t0 = t0 + h

    return y0

def Vetor_Ideal(metodo):
    #y0, t0, h, passos, f

    #DECLARACAO
    y0 = float(metodo[0])
    y = [y0]
    t0 = float(metodo[1])
    h = float(metodo[2])
    passos = int(metodo[3],10)
    funcao = parse_expr(metodo[4])
    y_ = symbols('y')
    t_ = symbols('t')

    def convert(Y,T):
        f = funcao.subs(y_, Y)
        f = f.subs(t_,T)
        return f
    #CALCULO DO EULER
    
    for i in range(passos):
        f = convert(y0,t0 + h*(i+1))
        y0 = f
        y.append(y0)
    
    return y

#MAIN
arq = open('entrada.txt', 'r')
entrada = arq.read()
saida = open('saida.txt','w+')
casos = entrada.splitlines()
casos.pop(0)
count = 0
ideal = []
for metodos in casos:
    count = count + 1
    metodos = metodos.split(' ')
    y = []
    t = []
    h  = metodos[-3]
    t0 = metodos[-4]
    passos = metodos[-2]
    print(metodos[0])
    if metodos[0] == 'erro':
        ideal = Vetor_Ideal(metodos[1:])
        print("HAHAHA")
        continue
    if metodos[0] == 'euler':
        saida.write('Método de Euler\n')
        y = Euler(metodos[1:])
    elif metodos[0] == 'euler_inverso':
        saida.write('Método de Euler Inverso\n')
        y = Euler_Inverso(metodos[1:])
    elif metodos[0] == 'euler_aprimorado':
        saida.write('Método de Euler Aprimorado\n')
        y = Euler_Aprimorado(metodos[1:])
    elif metodos[0] == 'runge_kutta':
        saida.write('Método de Range-Kutta\n')
        y = Runge_Kutta(metodos[1:])
    else:

        h  = metodos[-4]
        t0 = metodos[-5]
        passos = metodos[-3]
        
        if search('^adam_bashforth',metodos[0]):
            y = Adam_Bashforth(metodos[0], metodos[1:])    
        elif search('^adam_moulton',metodos[0]):
            y = Adam_Moulton(metodos[0], metodos[1:])
        elif search('^formula_inversa',metodos[0]):
            y = Formula_Inversa(metodos[0], metodos[1:])
        else:
            saida.write('Não existe\n')
            continue

        saida.write('\nOrdem = %s\n' %metodos[-1])
    #descobrir erro
    # for i in range(int(passos)):
        
    #     erro = (ideal[i]-y[i])/ideal[i]
    #     if(erro<0):
    #         erro = -erro
    #     print('ideal[ '+str(i)+'] = '+str(ideal[i]))
    #     print('y[ '+str(i)+'] = '+str(y[i]))
    #     print('erro[ '+str(i)+'] = '+ str(erro))
    erro = (ideal[int(passos)]-y[int(passos)])/ideal[int(passos)]
    if(erro<0):
        erro = -erro
    print('ideal[', int(passos), '] = ', "{:.6}".format(ideal[int(passos)]))
    print('y[', int(passos), '] = ', "{:.6}".format(y[int(passos)]))
    print('erro[', int(passos), '] = ', "{:.4%}".format(erro/100))

    t = [float(t0)]
    for i in range(int(passos)):
        t.append(t[-1] + float(h))
    
    plt.plot(t,y, label=('(' + str(count) + ') ' + metodos[0]))
    
    saida.write('y(%f) = ' %float(t0))
    saida.write('%f\n' %float(metodos[1]))
    saida.write('h = %f\n' %float(h))
    if(count%5==0):
        
        plt.ylabel('Quantidade de sal(kg)')
        plt.xlabel('Tempo(min)')
        plt.legend(loc='upper left',ncol = 2, fontsize = 'xx-small')
        plt.show()
    for i in range(int(passos)+1):
        saida.write ('%d\t'%i)
        saida.write('%f\n' %y[i])
    saida.write('\n')
plt.ylabel('Quantidade de sal(kg)')
plt.xlabel('Tempo(min)')
plt.legend(loc='upper left',ncol = 2, fontsize = 'xx-small')
plt.show()

arq.close()
saida.close()
