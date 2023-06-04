from scipy.optimize import least_squares
from PIL import Image
import numpy as np
from numpy import *
import time

surface_resolution = (35, 35)  # x, y
surface_dimensions = (10, 10)  # mm

image_path = 'img/lena.png'
light_direction = (1, 0, -0.2)

target_path = "generated/lena.obj"

max_iterations = 20

def create_mesh():
    dx = surface_dimensions[0] / surface_resolution[0]
    dy = surface_dimensions[1] / surface_resolution[1]
    vertices = []
    triangles = []

    # grid vertices
    for y in range(0, surface_resolution[1] + 1):
        for x in range(0, surface_resolution[0] + 1):
            vertices.append([x * dx, y * dy, 0])

    # center vertices
    for y in range(0, surface_resolution[1]):
        for x in range(0, surface_resolution[0]):
            vertices.append([(x + 0.5) * dx, (y + 0.5) * dy, 0])

    # add polygons
    offset = (surface_resolution[0] + 1) * (surface_resolution[1] + 1)
    for y in range(0, surface_resolution[1]):
        for x in range(0, surface_resolution[0]):
            upper_left = x + y * (surface_resolution[0] + 1)
            upper_right = (x + 1) + y * (surface_resolution[0] + 1)
            lower_left = x + (y + 1) * (surface_resolution[0] + 1)
            lower_right = (x + 1) + (y + 1) * (surface_resolution[0] + 1)
            center = offset + x + y * (surface_resolution[0])
            triangles.append([upper_left, upper_right, center])
            triangles.append([lower_left, upper_left, center])
            triangles.append([lower_right, lower_left, center])
            triangles.append([upper_right, lower_right, center])

    return np.array(vertices), triangles

def calculate_L(x, y, vertices, lightNormal):
    p_upper_left = p(x, y + 1, vertices)
    p_upper_right = p(x + 1, y + 1, vertices)
    p_lower_left = p(x, y, vertices)
    p_lower_right = p(x + 1, y, vertices)
    p_center = p_c(x, y, vertices)
    
    normal_1 = np.cross(np.subtract(p_center, p_lower_left), np.subtract(p_center, p_lower_right))
    normal_2 = np.cross(np.subtract(p_center, p_lower_right), np.subtract(p_center, p_upper_right))
    normal_3 = np.cross(np.subtract(p_center, p_upper_right), np.subtract(p_center, p_upper_left))
    normal_4 = np.cross(np.subtract(p_center, p_upper_left), np.subtract(p_center, p_lower_left))
    
    n_1 = np.linalg.norm(normal_1)
    n_2 = np.linalg.norm(normal_2)
    n_3 = np.linalg.norm(normal_3)
    n_4 = np.linalg.norm(normal_4)

    L_x, L_y, L_z = np.multiply(lightNormal, -1.0)

    A = (0.5*(L_x + L_y)*(1/n_1 + 1/n_4), 0.5*(-L_x + L_y)*(1/n_1 + 1/n_2), 0.5*(-L_x - L_y)*(1/n_2 + 1/n_3), 0.5*(L_x - L_y)*(1/n_3 + 1/n_4), L_x*(1/n_2 - 1/n_4) + L_y*(1/n_1 - 1/n_3))
    B = (p_lower_left[2], p_lower_right[2], p_upper_right[2], p_upper_left[2], p_center[2])
    A_transposed = np.transpose(A)
    return np.dot(A_transposed, B) + L_z*(1/n_1 + 1/n_2 + 1/n_3 + 1/n_4)
    

def compress_gradients(gradients, alpha_d=1):
    # Apply the compression function to each gradient value
    compressed_gradients = np.sign(gradients) * (1 / alpha_d) * np.log(1 + alpha_d * np.abs(gradients))
    return compressed_gradients

def weighs(b, x, y):
    return 1

def p(x, y, vertices):
    return vertices[x + y * (surface_resolution[0] + 1)]

def h(x, y, vertices):
    return p(x, y, vertices)[2]

def p_c(x, y, vertices):
    offset = (surface_resolution[0] + 1) * (surface_resolution[1] + 1)
    return vertices[offset + x + y * (surface_resolution[0])]

def h_c(x, y, vertices):
    return p_c(x, y, vertices)[2]

def dh(x, y, vertices):
    alpha_h = 1
    return alpha_h*log(1+alpha_h*h(x, y, vertices))

def solve_heightmap(vertices, I0, light_dir, max_iter):

    # Compress the image gradients using the specified compression function (eq 9)
    D0_x = compress_gradients(np.gradient(I0, axis=0))
    D0_y = compress_gradients(np.gradient(I0, axis=1))

    m = surface_resolution[0]
    n = surface_resolution[1]

    def compute_residuals(heights):
        vertices[:, 2] = heights

        # Calculate the actual radiance values (eq 6)
        L = []
        for y in range(0, n):
            L.append([])
            for x in range(0, m):
                L[y].append(calculate_L(x, y, vertices, light_dir))

        # save the luminance values as an image for visualization (remove for performance gain)
        im = Image.fromarray(np.asarray(np.abs(np.multiply(L, 4))))
        im = im.convert('RGB')
        im.save("generated/luminance.png")

        # Calculate the gradients of radiance values (eq 11)
        L0_x, L0_y = np.gradient(L)

        # Calculate squared differences between (eq 12)
        E0_x = weighs(0,0,0) * np.square(L0_x - D0_x)
        E0_y = weighs(0,0,0) * np.square(L0_y - D0_y)

        # Calculate the global energy for fitting the gradients (eq 13)
        E0_g = np.sum(E0_x[:, :m-1]) + np.sum(E0_y[:n-1, :])

        # Minimize sum of squared Laplacians (eq 14)
        E_c = 0
        for x in range(0, m):
            for y in range(0, n):
                E_c = E_c + np.square(h(x, y, vertices) + h(x + 1, y, vertices) + h(x + 1, y + 1, vertices) + h(x, y + 1, vertices) + 4*h_c(x, y, vertices))
        
        # Second order smoothness of corner vertices of each pixel
        E_p = 0
        for x in range(0, m):
            for y in range(0, n):
                E_p = E_p + np.square(h(x, y, vertices) - h(x + 1, y, vertices) + h(x + 1, y + 1, vertices) - h(x, y + 1, vertices))
        
        # Ensures the resulting heights are close to desired heights
        E_h = 0
        for x in range(0, m):
            for y in range(0, n):
                E_h = E_h + np.square(h(x, y, vertices) - dh(x, y, vertices))
        
        # calculate the total energy 
        E = E0_g# + 5*E_c + 5*E_p + 1.0*E_h

        return [E]

    # Use least squares optimization to find the optimal heights
    result = least_squares(compute_residuals, vertices[:, 2], ftol=1e-6, xtol=1e-6, gtol=1e-6, verbose=2, max_nfev=max_iter)

    # Retrieve the optimized heights
    vertices[:, 2] = result.x

    return vertices

def export_3dmesh(vertices, triangles, filename):
    with open(filename, 'w') as thefile:
        for vertex in vertices:
            thefile.write("v {0} {1} {2}\n".format(vertex[0], vertex[2], -vertex[1]))

        for triangle in triangles:
            thefile.write("f {0} {1} {2}\n".format(triangle[0] + 1, triangle[1] + 1, triangle[2] + 1))

image = Image.open(image_path)
image = image.resize(surface_resolution)

# Convert the image to grayscale
gray_image = image.convert("L")

# Get the pixel values as a 2D array
pixel_values = list(gray_image.getdata())
width, height = gray_image.size
pixel_values = [pixel_values[i * width:(i + 1) * width] for i in range(height)]
pixel_values = np.multiply(pixel_values, 1)

# Example usage
vertices, triangles = create_mesh()

# Solve the heights of the vertces
start_time = time.time()
vertices = solve_heightmap(vertices, pixel_values, light_direction, max_iterations)
print("--- took %s seconds ---" % round((time.time() - start_time), 2))

# Generate obj file
export_3dmesh(vertices, triangles, target_path)
