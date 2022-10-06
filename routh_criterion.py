import numpy as np
from PIL import Image as im
import sys

diferential_zero = sys.float_info.min

class Poly:
    def __init__(self, grau):
        self.grau = grau
        self.size = 0
        self.necessary_condition = 0
        self.coefs_list = 0
        self.non_null_coefs = 0
        self.impar_coefs = 0
        self.n_impar_coefs = 0
        self.par_coefs = 0
        self.n_par_coefs = 0

        self.size = self.grau + 1

        if self.grau <= 0 or not isinstance(self.grau, int):
            raise ValueError('O grau do polinômio deve ser inteiro e maior que zero')

        def input_coef():
            def print_mold():
                mold = ''
                for i in range(self.grau):
                    mold = mold + f'a_{grau - i}*s^{grau - i} + '
                mold = mold + f'a_{0}*s^{0}'
                print(mold)
            print_mold()
            coefs = []
            for i in range(self.size):
                a = float(input(f'a_{self.grau - i} = '))
                coefs.append(a)
            return coefs
        # a partir daqui tudo numpy.ndarray
        self.coefs_list = np.array(input_coef())
        def non_null_array (array):
            return array[array != 0]
        self.non_null_coefs = non_null_array(self.coefs_list)
        def impar_coefs ():
            impar_array = np.arange(self.size)
            impar_array = impar_array[(impar_array%2 != 0)]
            self.impar_coefs = non_null_array(np.take(self.coefs_list, [impar_array]))
        impar_coefs()
        def par_coefs():
            par_array = np.arange(self.size)
            par_array = par_array[(par_array%2 == 0)]
            self.par_coefs = non_null_array(np.take(self.coefs_list, [par_array]))
        par_coefs()


    def check_necessary_condition (self):
        if (np.all(self.coefs_list >= 0)):
            self.necessary_condition = True
        else:
            self.necessary_condition = False

    def print_poly(self):
        poly = ''
        for i in range(self.grau):
            poly = poly + f'{self.coefs_list[i]}s^{self.grau - i} + '
        poly = poly + f'{self.coefs_list[self.grau]}'
        print(poly)
        

class RouthTable:
    def __init__(self, poly):
        self.poly = poly
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
                self.special_case_1_flag = 1
                self.special_case_2_flag = 0
            elif ( np.all(self.routh_table[row_from_matrix,:] == 0) ):
                self.special_case_1_flag = 0
                self.special_case_2_flag = 1
        
        def resolve_special_case_1 (row):
            self.routh_rable[row,0] = diferential_zero

        def resolve_special_case_2 (row_from_matrix):
            def replace_row(row_from_matrix):
                poly = np.zeros(shape=(1,self.poly.size), dtype=float)[0]
                column_routh_table = 0
                # considerando que somente pegara a segunda linha (ímpares)
                for column in range(self.poly.size):
                    if (column%2 == 0):
                        poly[column] = 0
                    else:
                        poly[column] = self.routh_table[row_from_matrix-1,column_routh_table]
                        column_routh_table =+ 1

                derivative = np.array(np.poly1d(poly).deriv())
                non_null_derivative = derivative[derivative != 0]

                if non_null_derivative.size < self.routh_table_columns:
                    non_null_derivative = np.concatenate((non_null_derivative, [0]), axis = 0)
                                
                self.routh_table[row_from_matrix,:] = non_null_derivative
            replace_row(row_from_matrix)
            
            # trocar self.routh_table[:,0] por self.routh_first_column
            def count_signal_change ():
                sign_change_counter = 0
                for row in range(self.poly.size - 1):
                    if np.sign(self.routh_table[row,0]) != np.sign(self.routh_table[row+1,0]):
                        sign_change_counter =+ 1
            
            


        self.routh_table_rows = self.poly.size
        self.routh_table_columns = max(len(self.poly.par_coefs), len(self.poly.impar_coefs))

        pre_routh_table = np.zeros(shape=(self.poly.grau + 1, self.routh_table_columns), dtype=float)
        pre_routh_table[0,0:len(self.poly.par_coefs)] = self.poly.par_coefs
        pre_routh_table[1,0:len(self.poly.impar_coefs)] = self.poly.impar_coefs

        if (simplified == True):
            print('Simplified')
            self.routh_table = pre_routh_table.copy()
            self.routh_table[0,:] = simplify_row(self.routh_table[0,:])
            self.routh_table[1,:] = simplify_row(self.routh_table[1,:])

            for row in range(2, self.poly.size):
                stable_column = self.routh_table[(row-2):2+(row-2), [0]]
                for column in range(self.routh_table_columns-1):
                    variant_column = self.routh_table[(row-2):2+(row-2), [column+1]]
                    square_matrix = np.concatenate((stable_column,variant_column), axis=1)
                    self.routh_table[row,column] = calc(square_matrix)
                self.routh_table[row,:] = simplify_row(self.routh_table[row,:])

                check_special_cases(row)
                if (self.special_case_1_flag):
                    print(f'Special Case 1 is {self.special_case_1_flag}')
                    resolve_special_case_1(row)
                    self.special_case_1_flag = 0
                elif (self.special_case_2_flag):
                    print(f'Special Case 2 is {self.special_case_2_flag}')
                    resolve_special_case_2(row)
                    self.special_case_2_flag = 0

        else:
            print('Not simplified')
            self.routh_table = pre_routh_table.copy()
            for row in range(2, self.poly.size):
                stable_column = self.routh_table[(row-2):2+(row-2), [0]]
                for column in range(self.routh_table_columns-1):
                    variant_column = self.routh_table[(row-2):2+(row-2), [column+1]]
                    square_matrix = np.concatenate((stable_column,variant_column), axis=1)
                    self.routh_table[row,column] = calc(square_matrix)
                
                check_special_cases(row)
                if (self.special_case_1_flag):
                    print(f'Special Case 1 is {self.special_case_1_flag}')
                    resolve_special_case_1(row)
                    self.special_case_1_flag = 0
                elif (self.special_case_2_flag):
                    print(f'Special Case 2 is {self.special_case_2_flag}')
                    resolve_special_case_2(row)
                    self.special_case_2_flag = 0

    def check_sufficient_condition (self):
        if (np.all(self.routh_table[:,0] > 0)):
            self.sufficient_condition_flag = 1
        else:
            self.sufficient_condition_flag = 0



if __name__ == '__main__':
    poly = Poly(5)
    poly.print_poly()
    poly.check_necessary_condition()
    print(f'Necessary condition is {poly.necessary_condition}')

    routh_table = RouthTable(poly)
    routh_table.make_table(False)
    print(routh_table.routh_table)
    routh_table.check_sufficient_condition()
    print(f'Sufficient condition is {routh_table.sufficient_condition_flag}')