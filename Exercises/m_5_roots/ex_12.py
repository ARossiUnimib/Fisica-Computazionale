import numpy as np

def legendre_matrix(order):
    # For the Legendre polynomial, the matrix is typically tridiagonal
    # Here we construct a matrix that represents the recurrence relation for Legendre polynomials
    matrix = np.zeros((order, order))
    for i in range(order - 1):
        matrix[i, i + 1] = matrix[i + 1, i] = np.sqrt((i + 1) * (i + 2)) / (2 * i + 3)
    return matrix

def power_method(A, tol=1e-6, max_iter=1000):
    n = A.shape[0]
    x = np.random.rand(n)
    for _ in range(max_iter):
        # Apply the matrix A to vector x
        x_new = np.dot(A, x)
        
        # Normalize the vector
        x_new /= np.linalg.norm(x_new)
        
        # Compute the Rayleigh quotient (approximates the eigenvalue)
        eigenvalue = np.dot(x_new, np.dot(A, x_new)) / np.dot(x_new, x_new)
        
        # Check for convergence
        if np.linalg.norm(x_new - x) < tol:
            return eigenvalue, x_new
        
        x = x_new
    
    return eigenvalue, x

def deflation_power_method(A, num_eigenvalues, tol=1e-6, max_iter=1000):
    eigenvalues = []
    eigenvectors = []
    
    for _ in range(num_eigenvalues):
        # Compute the dominant eigenvalue and eigenvector
        eigenvalue, eigenvector = power_method(A, tol, max_iter)
        eigenvalues.append(eigenvalue)
        eigenvectors.append(eigenvector)
        
        # Deflate the matrix by removing the eigenvalue contribution
        # A_new = A - eigenvalue * outer_product of eigenvector with itself
        A -= eigenvalue * np.outer(eigenvector, eigenvector)
    
    return eigenvalues, eigenvectors

# Define the matrix for the Legendre polynomial of order 10
order = 10
A = legendre_matrix(order)

# Calculate the eigenvalues using the deflation power method
num_eigenvalues = 10
eigenvalues, eigenvectors = deflation_power_method(A, num_eigenvalues)

# Output the results
print("Eigenvalues of the Legendre polynomial (order 10):")
for i, eigenvalue in enumerate(eigenvalues):
    print(f"Eigenvalue {i + 1}: {eigenvalue}")

