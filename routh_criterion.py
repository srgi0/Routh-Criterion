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
            impar_array = impar_array[(impar_array%2 != 0)]
            self.impar_coefs = non_null_array(np.take(self.coefs, [impar_array]))
        impar_coefs()
        def par_coefs():
            par_array = np.arange(self.grau + 1)
            par_array = par_array[(par_array%2 == 0)]
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
        self.routh_table = 0
        self.routh_table_simplified = 0
        self.sufficient_condition = 0
    
    def make_table (self, *simplified):
        def calc (square_matrix):
            inv_det = (square_matrix[0,1]*square_matrix[1,0]) - (square_matrix[0,0]*square_matrix[1,1])
            return inv_det/square_matrix[1,0]

        max_column = max(len(poly.par_coefs), len(poly.impar_coefs))
        pre_routh_table = np.zeros(shape=(poly.grau + 1, max_column), dtype=float)
        pre_routh_table[0,0:len(poly.par_coefs)] = poly.par_coefs
        pre_routh_table[1,0:len(poly.impar_coefs)] = poly.impar_coefs
        
        if (simplified.count(True) > 0):
            self.routh_table_simplified = pre_routh_table.copy()
            for row in range(2, poly.grau+1):
                stable_column = self.routh_table_simplified[(row-2):2+(row-2), [0]]
                for column in range(max_column-1):
                    variant_column = self.routh_table_simplified[(row-2):2+(row-2), [column+1]]
                    square_matrix = np.concatenate((stable_column,variant_column), axis=1)
                    self.routh_table_simplified[row,column] = calc(square_matrix)

                if (np.all(self.routh_table_simplified[row,:]%2 == 0)):
                    self.routh_table_simplified[row,:] = self.routh_table_simplified[row,:]/2
        else:
            self.routh_table = pre_routh_table.copy()
            for row in range(2, poly.grau+1):
                stable_column = self.routh_table[(row-2):2+(row-2), [0]]
                for column in range(max_column-1):
                    variant_column = self.routh_table[(row-2):2+(row-2), [column+1]]
                    square_matrix = np.concatenate((stable_column,variant_column), axis=1)
                    self.routh_table[row,column] = calc(square_matrix)
        

    def print_table(self, *simplified):
        if (simplified.count(True) > 0):
            print('Simplified')
            print(self.routh_table_simplified)
        else:
            print('Not simplified')
            print(self.routh_table)

    def check_sufficient_condition (self):
        if (np.all(self.routh_table[:,0] > 0)):
            self.sufficient_condition = True
        else:
            self.sufficient_condition = False


poly = Poly(4)
poly.print_poly()
poly.check_necessary_condition()

routh_table = Routh_Table(poly)
routh_table.make_table(True)
print(routh_table.print_table(True))
