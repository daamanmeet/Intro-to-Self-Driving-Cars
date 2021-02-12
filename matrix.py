import math
from math import sqrt
import numbers

def zeroes(height, width):
        """
        Creates a matrix of zeroes.
        """
        g = [[0.0 for _ in range(width)] for __ in range(height)]
        return Matrix(g)

def identity(n):
        """
        Creates a n x n identity matrix.
        """
        I = zeroes(n, n)
        for i in range(n):
            I.g[i][i] = 1.0
        return I

class Matrix(object):

    # Constructor
    def __init__(self, grid):
        self.g = grid
        self.h = len(grid)
        self.w = len(grid[0])

    #
    # Primary matrix math methods
    #############################
 
    def determinant(self):
        """
        Calculates the determinant of a 1x1 or 2x2 matrix.
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate determinant of non-square matrix.")
        if self.h > 2:
            raise(NotImplementedError, "Calculating determinant not implemented for matrices largerer than 2x2.")
        
        
        # For 1X1 matrix
        if self.h == 1:
            return self.g[0]
        
        # For 2X2 matrix
        elif self.h == 2:
            if (self.g[0][0]*self.g[1][1])==(self.g[0][1]*self.g[1][0]):
                raise ValueError('The matrix is not invertible.')
            else:            
                result = self.g[0][0]*self.g[1][1]-self.g[0][1]*self.g[1][0]
                return result    
                   

    def trace(self):
        """
        Calculates the trace of a matrix (sum of diagonal entries).
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate the trace of a non-square matrix.")

        t_sum = 0 
        for value in range(self.h):
            t_sum += self.g[value][value]
            
        return t_sum
            

    def inverse(self):
        """
        Calculates the inverse of a 1x1 or 2x2 Matrix.
        """
        if not self.is_square():
            raise(ValueError, "Non-square Matrix does not have an inverse.")
        if self.h > 2:
            raise(NotImplementedError, "inversion not implemented for matrices larger than 2x2.")

        # inverse storage matrix
        inverse =[]
        
        # For 1X1 matrix 
        if self.h == 1:
            inverse.append([1/self.g[0][0]])
        
        # For 2X2 matrix
        elif self.h == 2:
            a = self.g[0][0]
            b = self.g[0][1]
            c = self.g[1][0]
            d = self.g[1][1]
            
            d_factor = 1/self.determinant()
            inverse = [[d,-b],
                       [-c,a]]
            for i in range(len(inverse)):
                for j in range(len(inverse[0])):
                    inverse[i][j] = d_factor*inverse[i][j]
        
        return Matrix(inverse)
        

    def T(self):
        """
        Returns a transposed copy of this Matrix.
        """
        # Transpose storage matrix
        t_matrix = []
        for c in range(self.w):
            t_row = []
            for r in range(self.h):
                t_row.append(self.g[r][c])
            t_matrix.append(t_row)
        return Matrix(t_matrix)

    def is_square(self):
        return self.h == self.w
    
      
    #
    # Begin Operator Overloading
    ############################
    def __getitem__(self,idx):
        """
        Defines the behavior of using square brackets [] on instances
        of this class.

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > my_matrix[0]
          [1, 2]

        > my_matrix[0][0]
          1
        """
        return self.g[idx]

    def __repr__(self):
        """
        Defines the behavior of calling print on an instance of this class.
        """
        s = ""
        for row in self.g:
            s += " ".join(["{} ".format(x) for x in row])
            s += "\n"
        return s

    def __add__(self,other):
        """
        Defines the behavior of the + operator
        """
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be added if the dimensions are the same") 

        """
        #Append_method
        
        p_matrix = []
        for r in range(self.h):
            p_row = []
            for c in range(self.w):
                p_row.append((self.g[r][c])+(other.g[r][c]))
            p_matrix.append(p_row)
            
        return Matrix(p_matrix)
        
        """        
        #storage matrix
        p_matrix = zeroes(self.h,self.w)
        for r in range(self.h):
            for c in range(self.w):
                p_matrix[r][c] = (self.g[r][c])+(other.g[r][c])
            
        return p_matrix
       

    def __neg__(self):
        """
        Defines the behavior of - operator (NOT subtraction)

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > negative  = -my_matrix
        > print(negative)
          -1.0  -2.0
          -3.0  -4.0
        """
        #   
        # Scalar multiplication with -1
        #
               
        #create a storage matrix
        n_matrix = zeroes(self.h,self.w)
        for r in range(self.h):
            for c in range(self.w):
                n_matrix[r][c] = -1*self.g[r][c]
            
        return n_matrix

    def __sub__(self, other):
        """
        Defines the behavior of - operator (as subtraction)
        """
        #   
        # Matrix subtraction
        #
        
        #create a storage matrix
        s_matrix = zeroes(self.h,self.w)
        for r in range(self.h):
            for c in range(self.w):
                s_matrix[r][c] = ((self.g[r][c])-(other.g[r][c]))
            
        return s_matrix

    def __mul__(self, other):
        """
        Defines the behavior of * operator (matrix multiplication)
        """
        #   
        # Matrix multiplication
        #
        
        # In order to perform Matrix multiplication first
        #check if n_cols of matrixA is equal to n_rows of matrixB
        # 
        
        if(self.w!=other.h):
            raise(ValueError, "matrixA no. of cols should be equal to matrixB no. of rows for matrix mul")
        elif(self.w==other.h):
            #create a storage matrix 
            """
            # m_A(mXn) * m_B(nXp) = m_C(mXp)
            #
            """
            m_matrix = zeroes(self.h,other.w)
            
            # Traverse thorugh MatrixM row
            for i in range(self.h):
                ## Traverse thorugh MatrixA row and MatrixB column
                for j in range(other.w):
                    #Traverse thorugh MatrixA column and MatrixB row
                    for k in range(other.h):
                        #update storage matrixM[i][j] = storage_matrix[i][j] + (MatrixA[i][k] * MatrixB[k][j])
                        m_matrix[i][j] += self.g[i][k]*other.g[k][j]
                        
            return m_matrix
            
            
            
    def __rmul__(self, other):
        """
        Called when the thing on the left of the * is not a _.

        Example:

        > identity = Matrix([ [1,0], [0,1] ])
        > doubled  = 2 * identity
        > print(doubled)
          2.0  0.0
          0.0  2.0
        """
        
        if isinstance(other, numbers.Number):
            pass
            #   
            # scalar multipliaction with other(number)
            #
            
            #storage matrix
            r_matrix = zeroes(self.h,self.w)
            for r in range(self.h):
                for c in range(self.w):
                    r_matrix[r][c] = other*(self.g[r][c])
                    
            return r_matrix
            