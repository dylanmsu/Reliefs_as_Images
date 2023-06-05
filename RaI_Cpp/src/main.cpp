#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

//#include <eigen3/Eigen/Core>
//#include "ceres/ceres.h"
//#include "glog/logging.h"


#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

/*using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;*/

#define RES_X 200
#define RES_Y 200

#define MMPP 1 //mm per pixel

int generate_mesh(std::vector<std::vector<double>> &vertices, std::vector<std::vector<int>> &triangles) {
    double dx = RES_X*MMPP;
    double dy = RES_Y*MMPP;

    // grid vertices
    for (int y=0; y<RES_Y + 1; y++){
        for (int x=0; x<RES_X + 1; x++){
            vertices.push_back({x*dx, y*dy, 0.0f});
        }
    }

    // center vertices
    for (int y=0; y<RES_Y; y++){
        for (int x=0; x<RES_X; x++){
            vertices.push_back({(x + 0.5f)*dx, (y + 0.5f)*dy, 0.0f});
        }
    }

    int offset = (RES_X + 1) * (RES_Y + 1);
    for (int y=0; y<RES_Y; y++){
        for (int x=0; x<RES_X; x++){
            int upper_left = x + y * (RES_X + 1);
            int upper_right = (x + 1) + y * (RES_X + 1);
            int lower_left = x + (y + 1) * (RES_X + 1);
            int lower_right = (x + 1) + (y + 1) * (RES_X + 1);
            int center = offset + x + y * (RES_X);

            triangles.push_back({upper_left, upper_right, center});
            triangles.push_back({lower_left, upper_left, center});
            triangles.push_back({lower_right, lower_left, center});
            triangles.push_back({upper_right, lower_right, center});
        }
    }

    return 0;
}

std::vector<double> cross(std::vector<double> vector_a, std::vector<double> vector_b) {
   std::vector<double> cross = {0.0f, 0.0f, 0.0f};
   cross[0] = vector_a[1]*vector_b[2] - vector_a[2]*vector_b[1];
   cross[1] = vector_a[2]*vector_b[0] - vector_a[0]*vector_b[2];
   cross[2] = vector_a[0]*vector_b[1] - vector_a[1]*vector_b[0];
   return cross;
}

std::vector<double> vector_subtract(std::vector<double> vector_a, std::vector<double> vector_b) {
   std::vector<double> difference = {0.0f, 0.0f, 0.0f};
   difference[0] = vector_a[0] - vector_b[0];
   difference[1] = vector_a[1] - vector_b[1];
   difference[2] = vector_a[2] - vector_b[2];
   return difference;
}

double norm(std::vector<double> vector) {
    double squared_norm = 0;
    for (int i=0; i<vector.size(); i++) {
        squared_norm += vector[i]*vector[i];
    }
    return sqrt(squared_norm);
}

double dot(std::vector<double> vector_a, std::vector<double> vector_b) {
    double dot_product = 0;
    if (vector_a.size() == vector_b.size()) {
        for (int i=0; i<vector_a.size(); i++) {
            dot_product += vector_a[i]*vector_b[i];
        }
        return dot_product;
    } else {
        std::cout << "ERROR: vectors need to be equally size for dot product. " << std::endl;
        return 0.0f;
    }
}

int sgn(double x) {
    if (x > 0.0f) {
        return 1;
    } else if (x < 0.0f) {
        return -1;
    } else {
        return 0;
    }
}

void export_obj(std::string path, std::vector<std::vector<double>> vertices, std::vector<std::vector<int>> triangles) {
    std::ofstream outfile;
    outfile.open(path);

    for (std::vector<double> vertex : vertices) {
        outfile << "v " << vertex[0] << " " << vertex[2] << " " << -vertex[1] << "\n"; 
    }

    for (std::vector<int> face : triangles) {
        outfile << "f " << face[0] + 1 << " " << face[1] + 1 << " " << face[2] + 1 << "\n"; 
    }

    std::cout << "worte " << path << std::endl;
}

std::vector<double> p(int x, int y, std::vector<std::vector<double>> vertices) {
    return vertices[x + y * (RES_X + 1)];
}

double h(int x, int y, std::vector<std::vector<double>> vertices){
    return p(x, y, vertices)[2];
}

std::vector<double> p_c(int x, int y, std::vector<std::vector<double>> vertices) {
    int offset = (RES_X + 1) * (RES_Y + 1);
    return vertices[offset + x + y * (RES_X)];
}
    
double h_c(int x, int y, std::vector<std::vector<double>> vertices) {
    return p_c(x, y, vertices)[2];
}

double calculate_L(int x, int y, std::vector<std::vector<double>> vertices, std::vector<double> light_normal) {
    std::vector<double> p_upper_left = p(x, y + 1, vertices);
    std::vector<double> p_upper_right = p(x + 1, y + 1, vertices);
    std::vector<double> p_lower_left = p(x, y, vertices);
    std::vector<double> p_lower_right = p(x + 1, y, vertices);
    std::vector<double> p_center = p_c(x, y, vertices);

    std::vector<double> normal_1 = cross(vector_subtract(p_center, p_lower_left), vector_subtract(p_center, p_lower_right));
    std::vector<double> normal_2 = cross(vector_subtract(p_center, p_lower_right), vector_subtract(p_center, p_upper_right));
    std::vector<double> normal_3 = cross(vector_subtract(p_center, p_upper_right), vector_subtract(p_center, p_upper_left));
    std::vector<double> normal_4 = cross(vector_subtract(p_center, p_upper_left), vector_subtract(p_center, p_lower_left));

    double n_1 = norm(normal_1);
    double n_2 = norm(normal_2);
    double n_3 = norm(normal_3);
    double n_4 = norm(normal_4);

    double L_x = -light_normal[0];
    double L_y = -light_normal[1];
    double L_z = -light_normal[2];

    std::vector<double> A = {
        0.5*(L_x + L_y)*(1.0f/n_1 + 1.0f/n_4), 
        0.5*(-L_x + L_y)*(1.0f/n_1 + 1.0f/n_2), 
        0.5*(-L_x - L_y)*(1.0f/n_2 + 1.0f/n_3), 
        0.5*(L_x - L_y)*(1.0f/n_3 + 1.0f/n_4), 
        L_x*(1.0f/n_2 - 1.0f/n_4) + L_y*(1.0f/n_1 - 1.0f/n_3)
    };

    std::vector<double> B = {
        p_lower_left[2], 
        p_lower_right[2], 
        p_upper_right[2], 
        p_upper_left[2], 
        p_center[2]
    };

    return dot(A, B) + L_z*(1.0f/n_1 + 1.0f/n_2 + 1.0f/n_3 + 1.0f/n_4);
}

void gradient(const std::vector<std::vector<double>> grid, std::vector<std::vector<double>>& gradX, std::vector<std::vector<double>>& gradY) {
    for (int y=0; y<RES_Y; y++) {
        gradX.push_back(std::vector<double>());
        for (int x=0; x<RES_X-1; x++) {
            gradX[y].push_back(grid[y][x + 1] - grid[y][x]);
        }
    }

    for (int y=0; y<RES_Y-1; y++) {
        gradY.push_back(std::vector<double>());
        for (int x=0; x<RES_X; x++) {
            gradY[y].push_back(grid[y + 1][x] - grid[y][x]);
        }
    }
}

void compress_gradients(const std::vector<std::vector<double>> gradient_input, std::vector<std::vector<double>>& gradient_output) {
    double alpha_d = 1.0f;

    for (int y=0; y<RES_Y; y++) {
        gradient_output.push_back(std::vector<double>());
        for (int x=0; x<RES_X; x++) {
            double compressed_grad = sgn(gradient_input[y][x]) * (1 / alpha_d) * log(1 + alpha_d * abs(gradient_input[y][x]));
            gradient_output[y].push_back(compressed_grad);
        }
    }
}

void image_to_grid(cv::Mat image, std::vector<std::vector<double>>& image_grid) {
    for (int i = 0; i < image.rows; ++i) {
        std::vector<double> row;
        for (int j = 0; j < image.cols; ++j) {
            cv::Vec3b intensity = image.at<cv::Vec3b>(i, j); // Retrieve BGR values of the pixel
            double b = intensity[0] / 255.0; // Normalize B component
            double g = intensity[1] / 255.0; // Normalize G component
            double r = intensity[2] / 255.0; // Normalize R component
            double gray = (0.299 * r) + (0.587 * g) + (0.114 * b); // Calculate grayscale value using luminosity method
            row.push_back(gray);
        }
        image_grid.push_back(row);
    }
}

void show_grid(std::vector<std::vector<double>> image_grid) {
    cv::Mat grayImg(RES_Y, RES_X, CV_8U);
    for (int y = 0; y < RES_Y; y++) {
        for (int x = 0; x < RES_X; x++) {
            grayImg.at<uchar>(y, x) = static_cast<uchar>(image_grid[y][x] * 255);
        }
    }

    cv::Mat grayImg8bit;
    cv::convertScaleAbs(grayImg, grayImg8bit);
    cv::imshow("Grayscale Image", grayImg8bit);
    int k = cv::waitKey(0);
}

/*struct CostFunctor {
   template <typename T>
   bool operator()(const T* const x, T* residual) const {
     residual[0] = 10.0 - x[0];
     return true;
   }
};*/

int main(int argc, char** argv) {
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> triangles;

    generate_mesh(vertices, triangles);

    // get image
    std::string image_path = cv::samples::findFile("../../img/lena.png");
    cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    if(img.empty())
    {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return 1;
    }

    cv::Mat scaledImg;
    cv::resize(img, scaledImg, cv::Size(RES_X, RES_Y), cv::INTER_LINEAR);

    // convert image to grayscale values
    std::vector<std::vector<double>> pixels;
    image_to_grid(scaledImg, pixels);

    show_grid(pixels);

    std::vector<std::vector<double>> D0_x;
    std::vector<std::vector<double>> D0_y;
    gradient(pixels, D0_x, D0_y);

    show_grid(D0_x);
    
    //optimization
    /*std::vector<double> light_dir = {1.0f, 0.0f, -0.2f};
    std::vector<std::vector<double>> luminances;
    for (int y=0; y<RES_Y; y++) {
        luminances.push_back(std::vector<double>());
        for (int x=0; x<RES_X; x++) {
            luminances[y].push_back(calculate_L(x, y, vertices, light_dir));
        }
    }

    std::vector<std::vector<double>> L0_x;
    std::vector<std::vector<double>> L0_y;
    gradient(luminances, L0_x, L0_y);

    compress_gradients(L0_x, L0_x);
    compress_gradients(L0_y, L0_y);

    export_obj("test.obj", vertices, triangles);*/

    /*google::InitGoogleLogging(argv[0]);

    // The variable to solve for with its initial value.
    double initial_x = 5.0;
    double x = initial_x;

    // Build the problem.
    Problem problem;

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    CostFunction* cost_function =
        new AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
    problem.AddResidualBlock(cost_function, nullptr, &x);

    // Run the solver!
    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    std::cout << "x : " << initial_x
                << " -> " << x << "\n";*/

    return 0;
}