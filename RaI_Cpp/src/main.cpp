#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>

//#include <eigen3/Eigen/Core>
//#include "ceres/ceres.h"
//#include "glog/logging.h"


#include <iostream>
#include <vector>
#include <fstream>

/*using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;*/

#define RES_X 20
#define RES_Y 20

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

/*struct CostFunctor {
   template <typename T>
   bool operator()(const T* const x, T* residual) const {
     residual[0] = 10.0 - x[0];
     return true;
   }
};*/

int main(int argc, char** argv)
{
    std::ofstream outfile;

    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> triangles;

    generate_mesh(vertices, triangles);

    outfile.open("test.obj");

    for (std::vector<double> vertex : vertices) {
        outfile << "v " << vertex[0] << " " << vertex[2] << " " << -vertex[1] << "\n"; 
    }

    for (std::vector<int> face : triangles) {
        outfile << "f " << face[0] + 1 << " " << face[1] + 1 << " " << face[2] + 1 << "\n"; 
    }

    std::cout << "worte test.obj" << std::endl;

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

    /*std::string image_path = cv::samples::findFile("../../img/lena.png");
    cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    if(img.empty())
    {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return 1;
    }
    cv::imshow("Display window", img);
    int k = cv::waitKey(0); // Wait for a keystroke in the window*/
    return 0;
}
