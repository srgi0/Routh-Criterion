import numpy as np
from PIL import Image as im
import sys

diferential_zero = sys.float_info.min

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
        # a partir daqui tudo numpy.ndarray
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


    def check_necessary_condition (self):
        if (np.all(self.coefs >= 0)):
            self.necessary_condition = True
        else:
            self.necessary_condition = False

    def print_poly(self):
        poly = ''
        for i in range(self.grau):
            poly = poly + f'{self.coefs[i]}s^{self.grau - i} + '
        poly = poly + f'{self.coefs[self.grau]}'
        print(poly)
        

class RouthTable:
    def __init__(self, Poly):
        self.Poly = Poly
        self.routh_table = 0
        self.sufficient_condition_flag = 0
        self.special_case_1_flag = 0
        self.special_case_2_flag = 0
        self.routh_table_rows = 0
        self.routh_table_columns = 0
    
    def make_table (self, simplified):
        def calc (square_matrix):
            inv_det = (square_matrix[0,1]*square_matrix[1,0]) - (square_matrix[0,0]*square_matrix[1,1])
            return inv_det/square_matrix[1,0]

        def simplify_row (row_from_matrix):
            row_from_matrix_simplified = row_from_matrix.copy()
            for div in (2,3,5,7,11):
                    while (np.all(row_from_matrix_simplified%div == 0)):
                        row_from_matrix_simplified = row_from_matrix_simplified/div
            return row_from_matrix_simplified

        def check_special_cases (row_from_matrix):
            if ( self.routh_table[row_from_matrix,0] == 0 ) and \
            not ( np.all(self.routh_table[row_from_matrix,:] == 0) ) :
                self.special_case_1_flag = True
                self.special_case_2_flag = False
            else:
                self.special_case_1_flag = False
                self.special_case_2_flag = True
        
        def resolve_special_case_1 (row):
            self.routh_rable[row,0] = diferential_zero

        def resolve_special_case_2 (row):
            derivative = np.array((np.poly1d(self.routh_table[row-1,:]).deriv()))
            
            
            
            if derivative.size < self.routh_table_columns:
                derivative = np.concatenate((derivative, [0]), axis = 0)
                            
            self.routh_table[row,:] = derivative

        self.routh_table_rows = self.Poly.grau + 1
        self.routh_table_columns = max(len(self.Poly.par_coefs), len(self.Poly.impar_coefs))

        pre_routh_table = np.zeros(shape=(self.Poly.grau + 1, self.routh_table_columns), dtype=float)
        pre_routh_table[0,0:len(self.Poly.par_coefs)] = self.Poly.par_coefs
        pre_routh_table[1,0:len(self.Poly.impar_coefs)] = self.Poly.impar_coefs

        if (simplified == True):
            print('Simplified')
            self.routh_table = pre_routh_table.copy()
            self.routh_table[0,:] = simplify_row(self.routh_table[0,:])
            self.routh_table[1,:] = simplify_row(self.routh_table[1,:])

            for row in range(2, self.Poly.grau+1):
                stable_column = self.routh_table[(row-2):2+(row-2), [0]]
                for column in range(self.routh_table_columns-1):
                    variant_column = self.routh_table[(row-2):2+(row-2), [column+1]]
                    square_matrix = np.concatenate((stable_column,variant_column), axis=1)
                    self.routh_table[row,column] = calc(square_matrix)
                self.routh_table[row,:] = simplify_row(self.routh_table[row,:])

        else:
            print('Not simplified')
            self.routh_table = pre_routh_table.copy()
            for row in range(2, self.Poly.grau+1):
                stable_column = self.routh_table[(row-2):2+(row-2), [0]]
                for column in range(self.routh_table_columns-1):
                    variant_column = self.routh_table[(row-2):2+(row-2), [column+1]]
                    square_matrix = np.concatenate((stable_column,variant_column), axis=1)
                    self.routh_table[row,column] = calc(square_matrix)
                
                check_special_cases(row)
                if (self.special_case_1_flag):
                    resolve_special_case_1(row)
                elif (self.special_case_2_flag):
                    resolve_special_case_2(row)


    def check_sufficient_condition (self):
            if (np.all(self.routh_table[:,0] > 0)):
                self.sufficient_condition_flag = True
            else:
                self.sufficient_condition_flag = False

    



poly = Poly(5)
poly.print_poly()
poly.check_necessary_condition()
print(f'Necessary condition is {poly.necessary_condition}')

routh_table = RouthTable(poly)
routh_table.make_table(False)
print(routh_table.routh_table)
routh_table.check_sufficient_condition()
print(f'Special Case 1 is {routh_table.special_case_1_flag}')
print(f'Special Case 2 is {routh_table.special_case_2_flag}')
print(f'Sufficient condition is {routh_table.sufficient_condition_flag}')