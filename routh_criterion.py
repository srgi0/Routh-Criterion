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

class RouthTable:
    def __init__(self, Poly):
        self.Poly = Poly
        self.routh_table = 0
        self.routh_table_simplified = 0
        self.sufficient_condition = 0
        self.special_case_1 = 0
        self.special_case_2 = 0
        self.routh_table_rows = 0
        self.routh_table_columns = 0
    
    def make_table (self, *simplified):
        def calc (square_matrix):
            inv_det = (square_matrix[0,1]*square_matrix[1,0]) - (square_matrix[0,0]*square_matrix[1,1])
            return inv_det/square_matrix[1,0]

        def simplify_row (row_from_matrix):
            row_from_matrix_simplified = row_from_matrix.copy()
            for div in (2,3,5,7,11):
                    while (np.all(row_from_matrix_simplified%div == 0)):
                        row_from_matrix_simplified = row_from_matrix_simplified/div
            return row_from_matrix_simplified

        self.routh_table_rows = self.Poly.grau + 1
        self.routh_table_columns = max(len(self.Poly.par_coefs), len(self.Poly.impar_coefs))

        pre_routh_table = np.zeros(shape=(self.Poly.grau + 1, self.routh_table_columns), dtype=float)
        pre_routh_table[0,0:len(self.Poly.par_coefs)] = self.Poly.par_coefs
        pre_routh_table[1,0:len(self.Poly.impar_coefs)] = self.Poly.impar_coefs
        
        if (simplified.count(True) > 0):
            self.routh_table_simplified = pre_routh_table.copy()
            self.routh_table_simplified[0,:] = simplify_row(self.routh_table_simplified[0,:])
            self.routh_table_simplified[1,:] = simplify_row(self.routh_table_simplified[1,:])

            for row in range(2, self.Poly.grau+1):
                stable_column = self.routh_table_simplified[(row-2):2+(row-2), [0]]
                for column in range(self.routh_table_columns-1):
                    variant_column = self.routh_table_simplified[(row-2):2+(row-2), [column+1]]
                    square_matrix = np.concatenate((stable_column,variant_column), axis=1)
                    self.routh_table_simplified[row,column] = calc(square_matrix)
                self.routh_table_simplified[row,:] = simplify_row(self.routh_table_simplified[row,:])

        else:
            self.routh_table = pre_routh_table.copy()
            for row in range(2, self.Poly.grau+1):
                stable_column = self.routh_table[(row-2):2+(row-2), [0]]
                for column in range(self.routh_table_columns-1):
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
    

    def check_special_case_1 (self):
        if (np.any(self.routh_table[:,0] == 0)):
            self.special_case_1 = True
        else:
            self.special_case_1 = False       

    def check_special_case_2 (self):
        for row in range(self.routh_table_rows):
            if (np.all(self.routh_table[row,:] == 0)):
                self.special_case_2 = True
            else:
                self.special_case_2 = False



poly = Poly(4)
poly.print_poly()
poly.check_necessary_condition()

routh_table = RouthTable(poly)
routh_table.make_table(True)
print(routh_table.print_table(True))

