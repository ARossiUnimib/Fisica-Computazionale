import numpy as np

def legendre_matrix(order):
    matrix = np.zeros((order, order))
    for i in range(order - 1):
        matrix[i, i + 1] = matrix[i + 1, i] = np.sqrt((i + 1) * (i + 2)) / (2 * i + 3)
    return matrix

def power_method(A, tol=1e-16, max_iter=10000000, alpha = 0):
    n = A.shape[0]
    x = np.random.rand(n)
    for _ in range(max_iter):
        A_shifted = A + alpha * np.identity(n)
        x_new = np.dot(A_shifted, x)
        
        x_new /= np.linalg.norm(x_new)
        
        eigenvalue = np.dot(x_new, np.dot(A, x_new)) / np.dot(x_new, x_new)
        
        if np.linalg.norm(x_new - x) < tol:
            return eigenvalue-alpha, x_new
        
        x = x_new
    
    return eigenvalue-alpha, x

def deflation_power_method(A, num_eigenvalues, tol=1e-6, max_iter=1000, alpha = 0):
    eigenvalues = []
    eigenvectors = []
    
    for _ in range(num_eigenvalues):
        eigenvalue, eigenvector = power_method(A, tol, max_iter, alpha)
        eigenvalues.append(eigenvalue)
        eigenvectors.append(eigenvector)
        
        A -= eigenvalue * np.outer(eigenvector, eigenvector)
    
    return eigenvalues, eigenvectors

order = 10
A = legendre_matrix(order)

num_eigenvalues = 10
eigenvalues, eigenvectors = deflation_power_method(A, num_eigenvalues, 0.5)

print("Eigenvalues of the Legendre polynomial (order 10):")
for i, eigenvalue in enumerate(eigenvalues):
    print(f"Eigenvalue {i + 1}: {eigenvalue}")

