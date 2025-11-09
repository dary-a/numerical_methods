import math

import numpy as np


class LinearAlgebraicSystems:
    def __init__(self, n, A, b):
        self.n = n  # matrix dimension
        self.A = A  # matrix
        self.b = b
        self.LU_decomposition()


    def LU_decomposition(self):
        self.L = np.eye(self.n)
        self.U = np.copy(self.A)
        self.p = np.arange(self.n)
        self.swaps = 0

        for k in range(self.n):  # K IS A LINE NUMBERS
            max_index_in_subcol = np.argmax(
                np.abs(self.U[self.p[k:], k]))  # index in minor, using p because some lines we swapped
            max_index_in_col = max_index_in_subcol + k
            if k != max_index_in_col:
                self.p[[max_index_in_col, k]] = self.p[[k, max_index_in_col]]
                self.L[[max_index_in_col, k]] = self.L[[k, max_index_in_col]]
                self.swaps += 1

            pivot = self.U[self.p[k], k]
            for i in range(k + 1, n):
                self.L[i, k] = self.U[self.p[i], k] / pivot  # L без перестановок
                self.U[self.p[i], k:] -= self.L[i, k] * self.U[self.p[k], k:]


    def determinant(self):
        self.determinant = np.prod(np.diag(self.U[self.p])) * (-1) ** self.swaps
        return self.determinant


    def linear_system_solution(self, A=None, b=None):
        self.A = A if A is not None else self.A
        self.b = b if b is not None else self.b
        y = np.zeros(self.n)
        self.x = np.zeros(self.n)

        for i in range(self.n):
            y[i] = self.b[self.p][i] - np.dot(self.L[i, :i], y[:i])
        for i in range(self.n - 1, -1, -1):
            self.x[i] = (y[i] - np.dot(self.U[self.p][i, i + 1:], self.x[i + 1:])) / self.U[self.p][i, i]
        return self.x


    def find_inverse_matrix(self):
        ones = np.eye(self.n)
        inverse_matrix = np.zeros((n, n))
        for i in range(self.n):
            inverse_matrix[:, i] = self.linear_system_solution(self.A, ones[i])
        return inverse_matrix


    def condition_number(self):
        matrix_norm = max(np.sum(np.abs(self.A), axis=1))
        inversed_matrix = self.find_inverse_matrix()
        inversed_matrix_norm = max(np.sum(np.abs(inversed_matrix), axis=1))
        return matrix_norm * inversed_matrix_norm

    def print_matrix(self):
        header = " ".join(f"{'a' + str(j + 1):>8}" for j in range(A.shape[1])) + " |     b"
        print(header)
        print("-" * len(header))
        for i in range(self.n):
            print(" ".join(f"{self.A[i, j]:8.3f}" for j in range(self.A.shape[1])), "|", f"{self.b[i]:8.3f}")
        print()


n = 4
A = np.random.randint(1, 10, size=(n, n)).astype(float)
b = np.random.randint(1, 10, size=n).astype(float)
s = LinearAlgebraicSystems(n, A, b)
s.print_matrix()
# print(f'My determinant: {s.determinant()}\nReal determinant: {np.linalg.det(A)}\n')
# print(f'My solution: {s.linear_system_solution()}\nReal solution: {np.linalg.solve(A, b)}')
# print(f'Inverse matrix for matrix A is \n{s.find_inverse_matrix()}')
print(f'Condition number for matrix A is {s.condition_number()}')
