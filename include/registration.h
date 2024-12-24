#ifndef REGISTRATION
#define REGISTRATION

#include "utils.h"
#include <random>
#include <chrono>
#include <iostream>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand


namespace HERE{

    class registration{

    public:

        registration();
        ~registration();

        bool (*distance)[17000];
        std::vector<size_t> sample_index;
        std::vector<double> consistence;
        double preserve_ratio = 0.15;

        size_t num_samples = 20;
        size_t num_seed = 1;


        std::vector<size_t> inlier_index_after_t_est;
        std::vector<size_t> inlier_index_after_r_est;
        std::vector<size_t> inlier_index_after_theta_est;

        std::vector<size_t> inliers_final_all;



        Eigen::Vector3d t_est;
        Eigen::Vector3d r_est;
        double theta_est;

        Eigen::Matrix3d R_est_trSVD;
        Eigen::Vector3d t_est_trSVD;

        Eigen::Matrix3d R_est_trthetaSVD;
        Eigen::Vector3d t_est_trthetaSVD;

        Eigen::Matrix3d est_rotation;
        Eigen::Vector3d est_translation;

        std::vector<int> consistence_int;


        void assign(bool (*distance_out)[17000]);
        void norm_pointcloud(cloudPoints& src_points, cloudPoints& tgt_points);
        void cal_distance(cloudPoints& src_points, cloudPoints& tgt_points, double xi);

        void est6DoF(cloudPoints& src_points, cloudPoints& tgt_points, double xi, bool if_normalization);
        void est3translation(cloudPoints &src_points, cloudPoints &tgt_points, double xi, cloudNorms &X_norm, cloudNorms &Y_norm, std::vector<std::vector<size_t>> &index_seed, std::vector<Eigen::Vector3d> & t_seed);
        void est2rotationaxis(cloudPoints &src_points, cloudPoints &tgt_points, double xi, std::vector<std::vector<size_t>> &index_seed, std::vector<Eigen::Vector3d> & r_seed);
        void est1rotationangle_wo_ro(cloudPoints &src_points, cloudPoints &tgt_points, double xi);


    protected:

        std::vector<Eigen::Vector3d> centers;
        double GlobalScale;


        std::vector<size_t> est_t3_bnb_bound(const int &index_sample, const cloudPoints &src_points, const cloudPoints& tgt_points,
                                             const Eigen::Matrix<double, 1, Eigen::Dynamic> &X_norm, const Eigen::Matrix<double, 1, Eigen::Dynamic> & Y_norm,
                                             utils::state1Dof state_init, double &phi_stabber, double xi, bool if_upper);

        bool circle_intersection(const Eigen::Vector2d &center2_1, const double dis_center, const double r_1, const double r_2,
                                 Eigen::Vector2d &inter);

        template <typename T>
        std::vector<size_t> sort_indexes(std::vector<T> &v);

        void interval_stabbing(std::vector<double> &intervals_phi, std::vector<size_t> &intervals_index, double &stabber, std::vector<size_t> &index_stabbed);

        void post_SVD(cloudPoints &src, cloudPoints &tgt, double xi, Eigen::Matrix3d &R_method, Eigen::Vector3d &t_method, Eigen::Matrix3d &R_post, Eigen::Vector3d &t_post);

        void post_refinement(cloudPoints &src, cloudPoints &tgt, double xi, Eigen::Matrix3d &R_method, Eigen::Vector3d &t_method, Eigen::Matrix3d &R_post, Eigen::Vector3d &t_post, int num);

    };




}


#endif
