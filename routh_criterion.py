import numpy as np
from PIL import Image as im

class Poly:
    def __init__(self, grau):
        self.grau = grau
        self.necessary_condition = 0
        self.coefs = 0
        self.non_null_coefs = 0
        self.impar_coefs = 0
        self.n_impar_coefs = 0
        self.par_coefs = 0
        self.n_par_coefs = 0

        if self.grau <= 0 or not isinstance(self.grau, int):
            raise ValueError('O grau do polinômio deve ser inteiro e maior que zero')

        def input_coef():
            def print_mold():
                mold = ''
                for i in range(grau):
                    mold = mold + f'a_{grau - i}*s^{grau - i} + '
                mold = mold + f'a_{0}*s^{0}'
                print(mold)
            print_mold()
            coefs = []
            for i in range(grau+1):
                a = float(input(f'a_{grau - i} = '))
                coefs.append(a)
            return coefs
        self.coefs = np.array(input_coef())
        def non_null_array (array):
            return array[array != 0]
        self.non_null_coefs = non_null_array(self.coefs)
        def impar_coefs ():
            impar_array = np.arange(self.grau + 1)
            impar_array = impar_array[(impar_array%2 == 0)]
            self.impar_coefs = non_null_array(np.take(self.coefs, [impar_array]))
        impar_coefs()
        def par_coefs():
            par_array = np.arange(self.grau + 1)
            par_array = par_array[(par_array%2 != 0)]
            self.par_coefs = non_null_array(np.take(self.coefs, [par_array]))
        par_coefs()

        

    def print_poly(self):
        poly = ''
        for i in range(self.grau):
            poly = poly + f'{self.coefs[i]}s^{self.grau - i} + '
        poly = poly + f'{self.coefs[self.grau]}'
        print(poly)

    def check_necessary_condition (self):
        if (np.all(self.coefs >= 0)):
            self.necessary_condition = True
        else:
            self.necessary_condition = False


class Routh_Table:
    def __init__(self, poly):
        self.poly = poly
        self.sufficient_condition = 0
    
    def b (square_matrix):
        inv_det = (square_matrix[0,1]*square_matrix[1,0]) - (square_matrix[0,0]*square_matrix[0,0])
        return 1/inv_det

    def make_table (self):
        routh_table = np.zeros(shape=(poly.grau + 1, max(len(poly.par_coefs), len(poly.impar_coefs))), dtype=float)
        routh_table[0,0:len(poly.par_coefs)] = poly.par_coefs
        routh_table[1,0:len(poly.impar_coefs)] = poly.impar_coefs

poly = Poly(5)
poly.print_poly()
poly.check_necessary_condition()
print(poly.necessary_condition)