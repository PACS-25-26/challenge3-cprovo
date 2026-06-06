import numpy as np

def f(x, y):
    return 8 * np.pi**2 * np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y)

def exact(x, y):
    return np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y)

n = 50
h = 1.0 / (n - 1)
u = np.zeros((n, n))
u_new = np.zeros((n, n))

# Initialize with exact on boundary
for i in range(n):
    for j in range(n):
        if i == 0 or i == n - 1 or j == 0 or j == n - 1:
            u[i, j] = exact(i * h, j * h)

for it in range(100):
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            u_new[i, j] = 0.25 * (u[i - 1, j] + u[i + 1, j] + u[i, j - 1] + u[i, j + 1] + h * h * f(i * h, j * h))
    u[1:-1, 1:-1] = u_new[1:-1, 1:-1]

print("Error with +h^2 f: ", np.max(np.abs(u - exact(np.arange(n)[:, None] * h, np.arange(n)[None, :] * h))))

for it in range(100):
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            u_new[i, j] = 0.25 * (u[i - 1, j] + u[i + 1, j] + u[i, j - 1] + u[i, j + 1] - h * h * f(i * h, j * h))
    u[1:-1, 1:-1] = u_new[1:-1, 1:-1]

print("Error with -h^2 f: ", np.max(np.abs(u - exact(np.arange(n)[:, None] * h, np.arange(n)[None, :] * h))))
