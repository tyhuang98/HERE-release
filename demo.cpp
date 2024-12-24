#include <iostream>
#include <Eigen/Core>
#include <chrono>
#include <vector>

#include <algorithm>
#include <random>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

//external libraries
#include <pcl/io/io.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>

#include "utils.h"
#include "simu_process.h"
#include "registration.h"


inline double getAngularError(Eigen::Matrix3d R_exp, Eigen::Matrix3d R_est) {
    return std::abs(std::acos(fmin(fmax(((R_exp.transpose() * R_est).trace() - 1) / 2, -1.0), 1.0)) * 180 / M_PI);
}

bool distance[17000][17000];


int main(){

    srand(unsigned(time(NULL)));

    // read original point cloud
    std::string ply_path = "../models/bun_zipper.ply";
//    std::string ply_path = "../models/Armadillo.ply";
    pcl::PointCloud<pcl::PointXYZ>::Ptr src_cloud_ori(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::io::loadPLYFile(ply_path, *src_cloud_ori);

    // max N = 17000 for both bunny and Armadillo
    int N = 10000;                          // correspondence number
    double tho = 0.95;                       // outlier ratio
    double noise_bound = 0.02;              // maximum length of noise vector

    // generate random Transformation
    Eigen::Matrix4d T = Eigen::Matrix4d::Zero();

    T << 0.803058, -0.195693, -0.562852, -0.384886,
            -0.282914, 0.706093, -0.649147, -1.38517,
            0.524459, 0.680542, 0.511669, -2.86804,
            0, 0, 0, 1;

    Eigen::Matrix3d Rotation = T.block<3,3>(0, 0);
    Eigen::Vector3d Translation = T.block<3,1>(0, 3);

    // data generation
    cloudPoints src_points;          // 3*N
    cloudPoints tgt_points;          // 3*N
    std::vector<bool> inlier_mask;   // true: inlier   false: outlier

    simu_process gen_data;
    gen_data.model_data_generation(src_cloud_ori, N, tho, T, noise_bound, src_points, tgt_points, inlier_mask);


    // visualization
    pcl::PointCloud<pcl::PointXYZ>::Ptr src_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    for (size_t i = 0; i < N; ++i) {
        src_cloud->push_back(
                pcl::PointXYZ(static_cast<double>(src_points.col(i)[0]), static_cast<double>(src_points.col(i)[1]),
                              static_cast<double>(src_points.col(i)[2])));
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr tgt_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    for (size_t i = 0; i < N; ++i) {
        tgt_cloud->push_back(
                pcl::PointXYZ(static_cast<double>(tgt_points.col(i)[0]), static_cast<double>(tgt_points.col(i)[1]),
                              static_cast<double>(tgt_points.col(i)[2])));
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr tgt_cloud_inlier(new pcl::PointCloud<pcl::PointXYZ>);
    for (size_t i = 0; i < N; ++i) {
        if(inlier_mask[i]){
            tgt_cloud_inlier->push_back(
                    pcl::PointXYZ(static_cast<double>(tgt_points.col(i)[0]), static_cast<double>(tgt_points.col(i)[1]),
                                  static_cast<double>(tgt_points.col(i)[2])));
        }
    }

    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("test"));
    viewer->setBackgroundColor(255, 255, 255);
    viewer->setCameraPosition(-7.088909,0.489558,6.150008,-0.408311,-1.214921,-2.231068,0.290319,0.956217,0.036946);

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> src_color(src_cloud, 255, 180, 0);
    viewer->addPointCloud(src_cloud, src_color, "src");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3.5, "src");

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> dst_color(tgt_cloud, 0, 166, 237);
    viewer->addPointCloud(tgt_cloud, dst_color, "dst");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2.5, "dst");

    viewer->addPointCloud(tgt_cloud_inlier, dst_color, "dst_inlier");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "dst_inlier");


    // visualize the correspondences
//    for (size_t j = 0; j < src_cloud->size(); ++j)
//    {
//        std::stringstream ss_line;
//        ss_line << "correspondence_line_" << j;
//        pcl::PointXYZ & src_keypoint = src_cloud->points[j];
//        pcl::PointXYZ & tgt_keypoint = tgt_cloud->points[j];
//
//        if(inlier_mask[j]){
//            viewer->addLine<pcl::PointXYZ, pcl::PointXYZ> (src_keypoint, tgt_keypoint, 0, 255, 0, ss_line.str ());
//            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 0.01, ss_line.str ());
//        }
//        else{
//            viewer->addLine<pcl::PointXYZ, pcl::PointXYZ> (src_keypoint, tgt_keypoint, 255, 0, 0, ss_line.str ());
//            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 0.01, ss_line.str ());
//        }
//    }

    viewer->spin();
    while (!viewer->wasStopped()) {
        viewer->spinOnce();
    }



    std::chrono::steady_clock::time_point HERE_begin = std::chrono::steady_clock::now();
    HERE::registration regis_HERE;
    regis_HERE.assign(distance);

    regis_HERE.preserve_ratio = 0.15;
    regis_HERE.num_samples = 20;
    regis_HERE.est6DoF(src_points, tgt_points, noise_bound, true);
    std::chrono::steady_clock::time_point HERE_end = std::chrono::steady_clock::now();


    Eigen::Matrix3d rotation_HERE = regis_HERE.est_rotation;
    Eigen::Vector3d translation_HERE = regis_HERE.est_translation;
    Eigen::Matrix4d T_HERE;
    T_HERE.block<3,3>(0, 0) = rotation_HERE;
    T_HERE.block<3,1>(0, 3) = translation_HERE;

    double R_error_HERE = getAngularError(Rotation, rotation_HERE);
    double t_error_HERE = (Translation - translation_HERE).norm();
    double time_cost_HERE = std::chrono::duration_cast<std::chrono::microseconds>(HERE_end - HERE_begin).count() / 1000000.0;

    std::cout << "Rotation error(HERE): " << R_error_HERE << std::endl;
    std::cout << "Translation error(HERE): " << t_error_HERE << std::endl;
    std::cout << "Time cost(HERE): " << time_cost_HERE << std::endl;


    pcl::visualization::PCLVisualizer::Ptr viewer_final(new pcl::visualization::PCLVisualizer("test"));
    viewer_final->setBackgroundColor(255, 255, 255);
    viewer_final->setCameraPosition(-7.088909,0.489558,6.150008,-0.408311,-1.214921,-2.231068,0.290319,0.956217,0.036946);

    pcl::PointCloud<pcl::PointXYZ>::Ptr src_cloud_trans(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::transformPointCloud(*src_cloud, *src_cloud_trans, T_HERE);

    viewer_final->addPointCloud(src_cloud_trans, src_color, "src_mesh");
    viewer_final->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2.5, "src_mesh");

    viewer_final->addPointCloud(tgt_cloud, dst_color, "dst_mesh");
    viewer_final->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2.5, "dst_mesh");

    viewer_final->addPointCloud(tgt_cloud_inlier, dst_color, "dst_inlier");
    viewer_final->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "dst_inlier");



    while (!viewer_final->wasStopped()) {
        viewer_final->spinOnce();
    }

    return 0;
}

