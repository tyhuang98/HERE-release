#include "registration.h"

#include <random>
#include "state_priority_hashtable.h"
#include "state_priority_queue.h"

using namespace HERE;
using namespace HERE::utils;

registration::registration() = default;

registration::~registration() = default;

void registration::assign(bool (*distance_out)[17000]){
    distance = distance_out;
}


void registration::norm_pointcloud(cloudPoints &src_points, cloudPoints &tgt_points) {

    double src_scale_max = src_points.colwise().norm().maxCoeff();
    double tgt_scale_max = tgt_points.colwise().norm().maxCoeff();

    double scale = std::max(src_scale_max, tgt_scale_max);
    GlobalScale = scale;

    src_points = src_points/GlobalScale;
    tgt_points = tgt_points/GlobalScale;

}


void registration::cal_distance(cloudPoints& src_points, cloudPoints& tgt_points, double xi){

    size_t N_corrs = src_points.cols();

    size_t count_max = 17000;
    size_t count = std::min(N_corrs, count_max);
    std::vector<double> consistence_double;

    for(int i = 0; i< count; i++){

        int num = 0;
        for(int j = i + 1; j < count; j++){

            float error = (src_points.col(i) - src_points.col(j)).norm() - (tgt_points.col(i) - tgt_points.col(j)).norm();

            if(abs(error) <= 0.5*xi) {
                num++;
                distance[i][j] = 1;
                distance[j][i] = 1;
            }
        }
        consistence_int.push_back(num);
    }


    for(int n = 0; n <count; n++){
        double score_n = 0.4*consistence_int[n];
        for(int m = 0; m < count ; m++){
            if(m != n && distance[n][m]) {
                score_n += 0.1 * consistence_int[m];
            }
        }

        consistence_double.push_back(score_n);
    }

    consistence = consistence_double;

    sample_index.resize(count);
    iota(sample_index.begin(), sample_index.end(), 0);
    sort(sample_index.begin(), sample_index.end(),
         [&consistence_double](size_t i1, size_t i2) {return consistence_double[i1] > consistence_double[i2]; });


}

void registration::est6DoF(cloudPoints& src_points, cloudPoints& tgt_points, double xi, bool if_normalization){


    if(if_normalization){
        norm_pointcloud(src_points, tgt_points);
        xi = xi/GlobalScale;
    }

    Eigen::Matrix<double, 1, Eigen::Dynamic> X_norm = src_points.colwise().norm();
    Eigen::Matrix<double, 1, Eigen::Dynamic> Y_norm = tgt_points.colwise().norm();

    std::chrono::steady_clock::time_point begin_trthetaSVD = std::chrono::steady_clock::now();
    cal_distance(src_points, tgt_points, xi);

    std::vector<std::vector<size_t>> index_after_t_seed;
    std::vector<Eigen::Vector3d> t_seed;

    est3translation(src_points, tgt_points, xi, X_norm, Y_norm, index_after_t_seed, t_seed);


    Eigen::Matrix3d R_final;
    Eigen::Vector3d t_final;
    std::vector<size_t> consensus_final;
    double error = 100000000;

    size_t num = index_after_t_seed.size();
    for(int seed_index = 0; seed_index < num; seed_index++){

        inlier_index_after_t_est = index_after_t_seed[seed_index];
        t_est = t_seed[seed_index];

        std::vector<std::vector<size_t>> index_after_tr_seed;
        std::vector<Eigen::Vector3d> r_seed;

        est2rotationaxis(src_points, tgt_points, xi, index_after_tr_seed, r_seed);

        size_t num_r = index_after_tr_seed.size();
        for(int seed_index_r = 0; seed_index_r < num_r; seed_index_r++){

            inlier_index_after_r_est = index_after_tr_seed[seed_index_r];
            r_est = r_seed[seed_index_r];

            est1rotationangle_wo_ro(src_points, tgt_points, 1.5*xi);

            size_t k_trtheta = inlier_index_after_theta_est.size();
            if(k_trtheta<3){
                continue;
            }

            std::vector<double> consistence_double = consistence;
            sort(inlier_index_after_theta_est.begin(), inlier_index_after_theta_est.end(),
                 [&consistence_double](size_t i1, size_t i2) {return consistence_double[i1] > consistence_double[i2]; });

            size_t num_all = 0;

            int seed_num = std::max(int(k_trtheta / 2), 3);
            Eigen::Matrix<double, 3, Eigen::Dynamic> inliers_src_trtheta(3, seed_num);
            Eigen::Matrix<double, 3, Eigen::Dynamic> inliers_tgt_trtheta(3, seed_num);
            for(int j = 0;j<seed_num;j++) {

                int index_j = inlier_index_after_theta_est[j];
                inliers_src_trtheta.col(j) = src_points.col(index_j);
                inliers_tgt_trtheta.col(j) = tgt_points.col(index_j);

                num_all++;
            }

            Eigen::Matrix4d transformation_matrix_trtheta;
            transformation_matrix_trtheta = pcl::umeyama(inliers_src_trtheta, inliers_tgt_trtheta);
            Eigen::Matrix3d R_est_trthetaSVD_pre = transformation_matrix_trtheta.block<3,3>(0, 0);
            Eigen::Vector3d t_est_trthetaSVD_pre = transformation_matrix_trtheta.block<3,1>(0, 3);

            std::vector<size_t> inliers_final;
            for(int n = 0; n<src_points.cols(); n++){
                double error = (tgt_points.col(n) - (R_est_trthetaSVD_pre*src_points.col(n) + t_est_trthetaSVD_pre)).norm();
                if(abs(error) < 0.8*xi){
                    inliers_final.push_back(n);
                }
            }

            size_t num_all_final = 0;
            Eigen::Matrix<double, 3, Eigen::Dynamic> inliers_src_final(3, inliers_final.size());
            Eigen::Matrix<double, 3, Eigen::Dynamic> inliers_tgt_final(3, inliers_final.size());
            for(int n_index = 0; n_index<inliers_final.size();n_index++){
                size_t index = inliers_final[n_index];
                inliers_src_final.col(n_index) = src_points.col(index);
                inliers_tgt_final.col(n_index) = tgt_points.col(index);

                num_all_final ++;
            }


            Eigen::Matrix4d transformation_matrix_final;
            transformation_matrix_final = pcl::umeyama(inliers_src_final, inliers_tgt_final);
            R_est_trthetaSVD = transformation_matrix_trtheta.block<3,3>(0, 0);
            t_est_trthetaSVD = transformation_matrix_trtheta.block<3,1>(0, 3);

            Eigen::Matrix<double, 3, Eigen::Dynamic> tgt_solve= (R_est_trthetaSVD*src_points).colwise() + t_est_trthetaSVD;
            auto error_norm = (tgt_solve - src_points).colwise().norm();
            double error_seed = error_norm.sum();

            if(error_seed < error){
                error = error_seed;
                R_final = R_est_trthetaSVD;
                t_final = t_est_trthetaSVD;
                consensus_final = inliers_final;
            }

        }

    }

    inliers_final_all = consensus_final;

    post_refinement(src_points, tgt_points, xi, R_final, t_final, est_rotation, est_translation, 5);

    if(if_normalization){
        est_translation = est_translation*GlobalScale;
    }

}



void registration::est3translation(cloudPoints &src_points, cloudPoints &tgt_points, double xi, cloudNorms &X_norm, cloudNorms &Y_norm, std::vector<std::vector<size_t>> &index_seed, std::vector<Eigen::Vector3d> & t_seed) {

    size_t num_corrs = src_points.cols();

    int num_iter = std::min(num_samples, num_corrs);
    num_samples = num_iter;

    int j;
    Eigen::VectorXd xi_vector(1, 1);
    xi_vector << xi;

    for(int i = -1; i <=1 ; i=i+2){

        std::vector<std::vector<size_t>> index_stabbing_poss;
        std::vector<size_t> index_stabbing_num;
        std::vector<Eigen::Vector3d> t_poss;


        X_norm = X_norm.colwise() + 0.5*i*xi_vector;
#pragma omp parallel for default(none) shared(sample_index, src_points, tgt_points, X_norm, Y_norm, \
                                       index_stabbing_num, index_stabbing_poss, t_poss, xi, num_iter) private(j) num_threads(12)
        for(j = 0; j < num_iter;j++){

            int index_sample = sample_index[j];

            // branch and bound search
            double t_3_up = tgt_points.col(index_sample)[2] + X_norm.col(index_sample).value();
            double t_3_lo = tgt_points.col(index_sample)[2] - X_norm.col(index_sample).value();

            state1Dof state_est;
            state1Dof state_init(t_3_lo, t_3_up);
            double phi_stabber_est;

            size_t lower_bound = 0;
            std::vector<size_t> lower_bound_index;
            size_t state_lower_bound;

            double phi_stabber_ini = 0;

            std::vector<size_t> upper_bound_index = est_t3_bnb_bound(index_sample, src_points, tgt_points, X_norm, Y_norm,
                                                                     state_init, phi_stabber_ini, xi, true);
            size_t upper_bound = upper_bound_index.size();

            if(upper_bound - lower_bound <= 0)
                continue;

            const int buckets = std::max((int)src_points.cols()/10, 10);

            StatePriorityHashtable<state1Dof, int> table(buckets);
            auto *sstate = new SearchState<state1Dof>(state_init, upper_bound);
            auto **state_array = new state1Dof*[2];
            table.push(sstate);

            int iter = 0;
            while (table.size()){

                iter++;

                sstate = table.pop();

                double phi_stabber = 0;
                std::vector<size_t> state_lower_bound_index = est_t3_bnb_bound(index_sample, src_points, tgt_points, X_norm, Y_norm,
                                                                               sstate->ssr, phi_stabber, xi, false);

                state_lower_bound = state_lower_bound_index.size();

                if(state_lower_bound > lower_bound){
                    lower_bound = state_lower_bound;
                    lower_bound_index = state_lower_bound_index;
                    state_est = sstate->ssr;
                    table.prune(lower_bound);
                    phi_stabber_est = phi_stabber;
                }

                if(sstate->bnd <= lower_bound) {
                    delete sstate;
                    break;
                }

                int np = split(sstate->ssr, state_array);
                delete sstate;
                for(int i = 0; i<np;++i){
                    double phi_stabber_upper = 0;
                    std::vector<size_t> upper_bound_in = est_t3_bnb_bound(index_sample, src_points, tgt_points, X_norm, Y_norm,
                                                                          *(state_array[i]), phi_stabber_upper, xi, true);
                    upper_bound = upper_bound_in.size();
                    if(upper_bound > lower_bound){
                        sstate = new SearchState<state1Dof>(*(state_array[i]), upper_bound);
                        table.push(sstate);
                    }
                    delete state_array[i];
                }
            }
            delete []state_array;

            int num_stab = lower_bound_index.size();
            Eigen::Vector3d t;
            t[2] = state_est.mean();
            double r_dis = tgt_points.col(index_sample)[2] - t[2];
            double r_new = sqrt(X_norm.col(index_sample).value() * X_norm.col(index_sample).value() - r_dis*r_dis);
            t[0] = tgt_points.col(index_sample)[0] + r_new * cos(phi_stabber_est);
            t[1] = tgt_points.col(index_sample)[1] + r_new * sin(phi_stabber_est);

            for(int index = 0; index < num_stab; index++){

                double angle_y_t = std::acos(((tgt_points.col(index_sample) - t).dot(tgt_points.col(lower_bound_index[index]) - t)) /
                                             ((tgt_points.col(index_sample) - t).norm() * (tgt_points.col(lower_bound_index[index]) - t).norm()));

                double angle_x = std::acos((src_points.col(index_sample).dot(src_points.col(lower_bound_index[index]))) /
                                           (src_points.col(index_sample).norm() * src_points.col(lower_bound_index[index]).norm()));

                double error = 0.5*std::asin(xi/src_points.col(index_sample).norm()) + 0.5*std::asin(xi/src_points.col(lower_bound_index[index]).norm());

                if(abs(angle_y_t - angle_x) > error){

                    lower_bound_index.erase(std::remove(lower_bound_index.begin(), lower_bound_index.end(), lower_bound_index[index]), lower_bound_index.end());

                }
            }

            lower_bound_index.insert(lower_bound_index.begin(), index_sample);
            lower_bound = lower_bound_index.size();

#pragma omp critical
            {
                index_stabbing_poss.push_back(lower_bound_index);
                t_poss.push_back(t);
                index_stabbing_num.push_back(lower_bound);
            }
        }

        std::vector<size_t> num_stabbing;
        num_stabbing.resize(index_stabbing_poss.size());
        iota(num_stabbing.begin(), num_stabbing.end(), 0);
        sort(num_stabbing.begin(), num_stabbing.end(),
             [&index_stabbing_num](size_t i1, size_t i2) {return index_stabbing_num[i1] > index_stabbing_num[i2]; });


        for(int k = 0; k < num_seed; k++){

            size_t index_k = num_stabbing[k];
            index_seed.push_back(index_stabbing_poss[index_k]);
            t_seed.push_back(t_poss[index_k]);
        }
    }


}


std::vector<size_t> registration::est_t3_bnb_bound(const int &index_sample, const cloudPoints &src_points, const cloudPoints& tgt_points,
                                                   const Eigen::Matrix<double, 1, Eigen::Dynamic> &X_norm, const Eigen::Matrix<double, 1, Eigen::Dynamic> & Y_norm,
                                                   utils::state1Dof state, double &phi_stabber, double xi, bool if_upper){


    std::vector<double> intervals_phi;
    std::vector<size_t> intervals_index;

    size_t num = src_points.cols() - 1;

    double t_3 = state.mean();

    Eigen::Vector2d center_1(tgt_points.col(index_sample)[0], tgt_points.col(index_sample)[1]);
    double r_dis = tgt_points.col(index_sample)[2] - t_3;
    double r_1 = sqrt(X_norm.col(index_sample).value() * X_norm.col(index_sample).value() - r_dis*r_dis);

    double xi_new = xi;
    if(if_upper){

        double y_3 = tgt_points.col(index_sample)[2];
        double r = X_norm.col(index_sample).value();

        double omega_u = abs(state.r_-y_3) <= r ? asin((state.r_ - y_3)/r) : M_PI_2;
        double omega_l = abs(state.l_-y_3) <= r ? asin((state.l_ - y_3)/r) : -M_PI_2;

        double sigma = 2 * r * sin((omega_u - omega_l)/4);
        assert(sigma > state.edge());
        xi_new = xi + sigma;
//          xi_new = xi + state.edge();
    }

    int num_com = int(num*preserve_ratio);
    for(int t =0;t<num_com;t++){
        int i = sample_index[t];
        if(i != index_sample){
            if(distance[index_sample][i]){
                Eigen::Vector2d  center_2(tgt_points.col(i)[0], tgt_points.col(i)[1]);
                Eigen::Vector2d  center_2_1 = center_2 - center_1;
                double dis_center = center_2_1.norm();

                double edge_1 = tgt_points.col(i)[2] + X_norm.col(i).value() + xi_new;
                double edge_2 = tgt_points.col(i)[2] + X_norm.col(i).value() - xi_new;
                double edge_3 = tgt_points.col(i)[2] - X_norm.col(i).value() + xi_new;
                double edge_4 = tgt_points.col(i)[2] - X_norm.col(i).value() - xi_new;

                if(xi_new < X_norm.col(i).value()){
                    if (t_3 > edge_3 && t_3 < edge_2){

                        double r_2_dis = tgt_points.col(i)[2] - t_3;

                        double r_2_in = sqrt((X_norm.col(i).value()-xi_new) * (X_norm.col(i).value()-xi_new) - r_2_dis*r_2_dis);
                        double r_2_out = sqrt((X_norm.col(i).value()+xi_new) * (X_norm.col(i).value()+xi_new) - r_2_dis*r_2_dis);

                        bool inter_in, inter_out;
                        Eigen::Vector2d phi_inter_in;
                        Eigen::Vector2d phi_inter_out;


                        inter_in = circle_intersection(center_2_1, dis_center, r_1, r_2_in, phi_inter_in);
                        inter_out = circle_intersection(center_2_1, dis_center, r_1, r_2_out, phi_inter_out);


                        if(inter_out){
                            if(inter_in){

                                if(phi_inter_in[0] >= phi_inter_out[0] && phi_inter_in[1] <= phi_inter_out[1]){
                                    if(phi_inter_in[0] - phi_inter_out[0] != 0){

                                        if(phi_inter_out[0] < 0){
                                            if(phi_inter_in[0] <= 0){
                                                intervals_phi.push_back(phi_inter_out[0] + 2*M_PI);
                                                intervals_phi.push_back(phi_inter_in[0] + 2*M_PI);
                                                intervals_index.push_back(i);
                                                intervals_index.push_back(i);
                                            }
                                            else{
                                                intervals_phi.push_back(phi_inter_out[0] + 2*M_PI);
                                                intervals_phi.push_back(2*M_PI);
                                                intervals_index.push_back(i);
                                                intervals_index.push_back(i);

                                                intervals_phi.push_back(0);
                                                intervals_phi.push_back(phi_inter_in[0]);
                                                intervals_index.push_back(i);
                                                intervals_index.push_back(i);
                                            }
                                        }
                                        else{
                                            intervals_phi.push_back(phi_inter_out[0]);
                                            intervals_phi.push_back(phi_inter_in[0]);
                                            intervals_index.push_back(i);
                                            intervals_index.push_back(i);
                                        }

                                        if(phi_inter_out[1] > 2*M_PI){
                                            if(phi_inter_in[1] >= 2*M_PI){
                                                intervals_phi.push_back(phi_inter_in[1] - 2*M_PI);
                                                intervals_phi.push_back(phi_inter_out[1] - 2*M_PI);
                                                intervals_index.push_back(i);
                                                intervals_index.push_back(i);
                                            }
                                            else{
                                                intervals_phi.push_back(phi_inter_in[1]);
                                                intervals_phi.push_back(2*M_PI);
                                                intervals_index.push_back(i);
                                                intervals_index.push_back(i);

                                                intervals_phi.push_back(0);
                                                intervals_phi.push_back(phi_inter_out[1] - 2*M_PI);
                                                intervals_index.push_back(i);
                                                intervals_index.push_back(i);
                                            }
                                        }
                                        else{
                                            intervals_phi.push_back(phi_inter_in[1]);
                                            intervals_phi.push_back(phi_inter_out[1]);
                                            intervals_index.push_back(i);
                                            intervals_index.push_back(i);
                                        }

                                    }
                                }
                                else{
                                    if(phi_inter_in[0] < 0){
                                        intervals_phi.push_back(phi_inter_in[1]);
                                        intervals_phi.push_back(phi_inter_in[0] + 2*M_PI);
                                        intervals_index.push_back(i);
                                        intervals_index.push_back(i);
                                    }
                                    else if(phi_inter_in[1] > 2*M_PI){
                                        intervals_phi.push_back(phi_inter_in[1] - 2*M_PI);
                                        intervals_phi.push_back(phi_inter_in[0]);
                                        intervals_index.push_back(i);
                                        intervals_index.push_back(i);
                                    }
                                    else{
                                        intervals_phi.push_back(0);
                                        intervals_phi.push_back(phi_inter_in[0]);
                                        intervals_index.push_back(i);
                                        intervals_index.push_back(i);

                                        intervals_phi.push_back(phi_inter_in[1]);
                                        intervals_phi.push_back(2*M_PI);
                                        intervals_index.push_back(i);
                                        intervals_index.push_back(i);
                                    }


                                }

                            }
                            else{

                                if(phi_inter_out[0] < 0){
                                    intervals_phi.push_back(phi_inter_out[0] + 2*M_PI);
                                    intervals_phi.push_back(2*M_PI);
                                    intervals_index.push_back(i);
                                    intervals_index.push_back(i);

                                    intervals_phi.push_back(0);
                                    intervals_phi.push_back(phi_inter_out[1]);
                                    intervals_index.push_back(i);
                                    intervals_index.push_back(i);
                                }
                                else if(phi_inter_out[1] > 2*M_PI){
                                    intervals_phi.push_back(phi_inter_out[0]);
                                    intervals_phi.push_back(2*M_PI);
                                    intervals_index.push_back(i);
                                    intervals_index.push_back(i);

                                    intervals_phi.push_back(0);
                                    intervals_phi.push_back(phi_inter_out[1] - 2*M_PI);
                                    intervals_index.push_back(i);
                                    intervals_index.push_back(i);
                                }
                                else{
                                    intervals_phi.push_back(phi_inter_out[0]);
                                    intervals_phi.push_back(phi_inter_out[1]);
                                    intervals_index.push_back(i);
                                    intervals_index.push_back(i);
                                }
                            }
                        }

                    }
                    else if((t_3 < edge_1 && t_3 > edge_2) || (t_3 > edge_4 && t_3 < edge_3)){

                        double r_2_dis = tgt_points.col(i)[2] - t_3;
                        double r_2 = sqrt((X_norm.col(i).value()+xi_new) * (X_norm.col(i).value()+xi_new) - r_2_dis*r_2_dis);

                        bool inter;
                        Eigen::Vector2d phi_inter;

                        inter = circle_intersection(center_2_1, dis_center, r_1, r_2, phi_inter);

                        if(inter){

                            if(phi_inter[0] < 0){
                                intervals_phi.push_back(phi_inter[0] + 2*M_PI);
                                intervals_phi.push_back(2*M_PI);
                                intervals_index.push_back(i);
                                intervals_index.push_back(i);

                                intervals_phi.push_back(0);
                                intervals_phi.push_back(phi_inter[1]);
                                intervals_index.push_back(i);
                                intervals_index.push_back(i);
                            }
                            else if(phi_inter[1] > 2*M_PI){
                                intervals_phi.push_back(phi_inter[0]);
                                intervals_phi.push_back(2*M_PI);
                                intervals_index.push_back(i);
                                intervals_index.push_back(i);

                                intervals_phi.push_back(0);
                                intervals_phi.push_back(phi_inter[1] - 2*M_PI);
                                intervals_index.push_back(i);
                                intervals_index.push_back(i);
                            }
                            else{
                                intervals_phi.push_back(phi_inter[0]);
                                intervals_phi.push_back(phi_inter[1]);
                                intervals_index.push_back(i);
                                intervals_index.push_back(i);
                            }

                        }
                    }
                }
                else{

                    double r_2_dis = tgt_points.col(i)[2] - t_3;
                    double r_2 = sqrt((X_norm.col(i).value()+xi_new) * (X_norm.col(i).value()+xi_new) - r_2_dis*r_2_dis);

                    bool inter;
                    Eigen::Vector2d phi_inter;

                    inter = circle_intersection(center_2_1, dis_center, r_1, r_2, phi_inter);

                    if(inter){

                        if(phi_inter[0] < 0){
                            intervals_phi.push_back(phi_inter[0] + 2*M_PI);
                            intervals_phi.push_back(2*M_PI);
                            intervals_index.push_back(i);
                            intervals_index.push_back(i);

                            intervals_phi.push_back(0);
                            intervals_phi.push_back(phi_inter[1]);
                            intervals_index.push_back(i);
                            intervals_index.push_back(i);
                        }
                        else if(phi_inter[1] > 2*M_PI){
                            intervals_phi.push_back(phi_inter[0]);
                            intervals_phi.push_back(2*M_PI);
                            intervals_index.push_back(i);
                            intervals_index.push_back(i);

                            intervals_phi.push_back(0);
                            intervals_phi.push_back(phi_inter[1] - 2*M_PI);
                            intervals_index.push_back(i);
                            intervals_index.push_back(i);
                        }
                        else{
                            intervals_phi.push_back(phi_inter[0]);
                            intervals_phi.push_back(phi_inter[1]);
                            intervals_index.push_back(i);
                            intervals_index.push_back(i);
                        }
                    }
                }
            }

        }
    }

    std::vector<size_t> index_stabbed;
    interval_stabbing(intervals_phi, intervals_index, phi_stabber, index_stabbed);

    return index_stabbed;
}



bool registration::circle_intersection(const Eigen::Vector2d &center2_1, const double dis_center, const double r_1, const double r_2,
                                       Eigen::Vector2d &inter){

    bool if_inter = false;

    if(dis_center > r_1 + r_2){
        return if_inter;
    }
    else if(dis_center > abs(r_1 - r_2)){

        if_inter = true;
        double omega_1 = acos(center2_1[0] / dis_center);
        double omega = center2_1[1]>=0 ? omega_1 : (2*M_PI - omega_1);
        double theta = acos((dis_center*dis_center + r_1*r_1 - r_2*r_2) / (2*dis_center*r_1));

        inter[0] = omega - theta;
        inter[1] = omega + theta;

    }
    else if(r_2 > r_1){
        if_inter = true;
        inter[0] = 0;
        inter[1] = 2*M_PI;
    }

    return if_inter;
}


void registration::est2rotationaxis(cloudPoints &src_points, cloudPoints &tgt_points, double xi, std::vector<std::vector<size_t>> &index_seed, std::vector<Eigen::Vector3d> & r_seed){

    size_t num_corrs = inlier_index_after_t_est.size();

    int num_iter = std::min(num_samples, num_corrs);

    std::vector<double> rotate_all = {0};
    for(int index_m = 0; index_m<rotate_all.size();index_m++){

        double rotate = rotate_all[index_m];

        int k;
        size_t num_stabbing_result = 0;
        std::vector<size_t> index_stabbing_result;
        double theta_result;
        Eigen::Vector3d m_k_result;

#pragma omp parallel for default(none) shared(num_iter, num_corrs, rotate, src_points, tgt_points, num_stabbing_result, index_stabbing_result, theta_result, m_k_result, xi, distance) private(k) num_threads(12)
        for(k = 0; k<num_iter; k++){

            size_t index_k = inlier_index_after_t_est[k];
            Eigen::Vector3d m_k_center = tgt_points.col(index_k) - t_est - src_points.col(index_k);
            Eigen::Vector3d m_k_center_normlize = m_k_center.normalized();
            double gama_cos = std::min(xi / m_k_center.norm(), 1.0);
            double game_sin = std::sqrt(1 - gama_cos*gama_cos);
            Eigen::Vector3d  m_k_center_cos = gama_cos*m_k_center_normlize;

            Eigen::Vector3d rotate_vector;
            rotate_vector[0] = std::cos(0);
            rotate_vector[2] = std::sin(0);
            rotate_vector[1] = (- m_k_center_normlize[0] * rotate_vector[0] - m_k_center_normlize[2] * rotate_vector[2]) / m_k_center_normlize[1];

            Eigen::Vector3d m_k = (m_k_center_cos + rotate*game_sin*rotate_vector.normalized()).normalized();

            double m_k1_k1 = m_k[0] * m_k[0];
            double m_k2_k2 = m_k[1] * m_k[1];
            double m_k3_k3 = 1 - m_k1_k1 - m_k2_k2;
            double m_k1_k2 = m_k[0] * m_k[1];
            double m_k2_k3 = m_k[1] * m_k[2];
            double m_k1_k3 = m_k[0] * m_k[2];

            std::vector<double> intervals_theta;
            std::vector<size_t> intervals_index;

            for(int i = 0; i < num_corrs; i++){
                size_t index_i = inlier_index_after_t_est[i];
                if(i != k && distance[index_i][index_k]){
                    double z = (tgt_points.col(index_i) - t_est - src_points.col(index_i)).norm();
                    if(z >= xi){

                        double xi_i = cos(asin(xi / z));
                        double xi_i_2 = xi_i * xi_i;

                        Eigen::Vector3d  m_i = (tgt_points.col(index_i) - t_est - src_points.col(index_i)).normalized();
                        double m_i1_i1 = m_i[0] * m_i[0];
                        double m_i2_i2 = m_i[1] * m_i[1];
                        double m_i3_i3 = 1 - m_i1_i1 - m_i2_i2;
                        double m_i1_i2 = m_i[0] * m_i[1];
                        double m_i2_i3 = m_i[1] * m_i[2];
                        double m_i1_i3 = m_i[0] * m_i[2];

                        double a_2 = m_k2_k2 * (1 - xi_i_2 - m_i3_i3) + m_k3_k3*(1-xi_i_2 - m_i2_i2) + 2*m_i2_i3*m_k2_k3;
                        double b = (1-xi_i_2-m_i3_i3)*m_k1_k2 - m_i1_i2*m_k3_k3 + m_i1_i3*m_k2_k3 + m_i2_i3*m_k1_k3;
                        double delta = (xi_i_2 - 1)*m_k3_k3*(-xi_i_2+m_i1_i1*m_k1_k1 + m_i2_i2*m_k2_k2 + m_i3_i3*m_k3_k3 + 2*m_i1_i2*m_k1_k2 + 2*m_i1_i3*m_k1_k3 + 2*m_i2_i3*m_k2_k3);

                        if(a_2 >= 0){
                            if(delta >= 0){
                                double theta_l = atan((-b- sqrt(delta))/(a_2));
                                double theta_r = atan((-b+ sqrt(delta))/(a_2));

                                if(theta_l > theta_r){
                                    break;
                                }
                                intervals_theta.push_back(-M_PI_2);
                                intervals_theta.push_back(theta_l);
                                intervals_index.push_back(index_i);
                                intervals_index.push_back(index_i);

                                intervals_theta.push_back(theta_r);
                                intervals_theta.push_back(M_PI_2);
                                intervals_index.push_back(index_i);
                                intervals_index.push_back(index_i);
                            }
                            else{
                                intervals_theta.push_back(-M_PI_2);
                                intervals_theta.push_back(M_PI_2);
                                intervals_index.push_back(index_i);
                                intervals_index.push_back(index_i);
                            }
                        }
                        else{

                            double theta_l = atan((-b- sqrt(delta))/(a_2));
                            double theta_r = atan((-b+ sqrt(delta))/(a_2));

                            if(theta_l < theta_r){
                                break;
                            }
                            intervals_theta.push_back(theta_r);
                            intervals_theta.push_back(theta_l);
                            intervals_index.push_back(index_i);
                            intervals_index.push_back(index_i);
                        }

                    }
                    else{
                        intervals_theta.push_back(-M_PI_2);
                        intervals_theta.push_back(M_PI_2);
                        intervals_index.push_back(index_i);
                        intervals_index.push_back(index_i);
                    }
                }

            }

            double theta_stabber;
            std::vector<size_t> index_stabbed;
            interval_stabbing(intervals_theta, intervals_index, theta_stabber, index_stabbed);
            index_stabbed.insert(index_stabbed.begin(), index_k);
            int num_stabbed = index_stabbed.size();

#pragma omp critical
            {
                if (num_stabbed > num_stabbing_result) {
                    theta_result = theta_stabber;
                    num_stabbing_result = num_stabbed;
                    index_stabbing_result = index_stabbed;
                    m_k_result = m_k;
                }
            }

        }

        double alpha = atan((m_k_result[0]*cos(theta_result) + m_k_result[1]*sin(theta_result))/(-m_k_result[2]));
        Eigen::Vector3d r;
        r << cos(alpha)*cos(theta_result), cos(alpha)*sin(theta_result), sin(alpha);
        r_seed.push_back(r.normalized());

        index_seed.push_back(index_stabbing_result);
    }

}

void registration::est1rotationangle_wo_ro(cloudPoints &src_points, cloudPoints &tgt_points, double xi){

    std::vector<double> intervals_angle;
    std::vector<size_t> intervals_index;
    int num_iter = inlier_index_after_r_est.size();
    int i;
    double xi_2 = xi*xi;

#pragma omp parallel for default(none) shared(num_iter, src_points, tgt_points, intervals_angle, intervals_index, xi, xi_2) private(i) num_threads(12)
    for(i=0; i<num_iter;i++){

        int index_i = inlier_index_after_r_est[i];
        Eigen::Vector3d src_vector = src_points.col(index_i);
        Eigen::Vector3d dst_vector = (tgt_points.col(index_i)-t_est);

        double src_proj_r = src_vector.transpose() * r_est;
        double dst_proj_r = dst_vector.transpose() * r_est;

        if(src_proj_r >= dst_proj_r - xi  && src_proj_r <= dst_proj_r + xi){

            Eigen::Vector3d center_1 = src_proj_r * r_est;
            double r_1 = sqrt(src_vector.squaredNorm() - src_proj_r*src_proj_r);
            Eigen::Vector3d center_2 = dst_vector - (dst_proj_r - src_proj_r)*r_est;
            double r_2 = sqrt(xi_2 - (dst_proj_r - src_proj_r)*(dst_proj_r - src_proj_r));

            double dis_center = (center_1 - center_2).norm();

            Eigen::Vector2d inter;
            bool if_inter = false;
            if(dis_center < r_1 + r_2) {
                if (dis_center > abs(r_1 - r_2)) {
                    if_inter = true;
                    Eigen::Vector3d src_proj_l = src_vector - center_1;
                    Eigen::Vector3d dst_proj_l = center_2 - center_1;

                    double norm_2 = src_proj_l.norm() * dst_proj_l.norm();
                    double omega_1 = acos(src_proj_l.dot(dst_proj_l) / norm_2);
                    double omega = dst_proj_l.transpose() * r_est.cross(src_vector) >= 0 ? omega_1 : (2 * M_PI - omega_1);
                    double theta = acos((dis_center * dis_center + r_1 * r_1 - r_2 * r_2) / (2 * dis_center * r_1));

                    inter[0] = omega - theta;
                    inter[1] = omega + theta;
                } else if (r_2 > r_1) {
                    if_inter = true;
                    inter[0] = 0;
                    inter[1] = 2 * M_PI;
                }
            }

#pragma omp critical
            {
                if (if_inter) {
                    if (inter[0] < 0) {
                        intervals_angle.push_back(inter[0] + 2 * M_PI);
                        intervals_angle.push_back(2 * M_PI);
                        intervals_index.push_back(index_i);
                        intervals_index.push_back(index_i);

                        intervals_angle.push_back(0);
                        intervals_angle.push_back(inter[1]);
                        intervals_index.push_back(index_i);
                        intervals_index.push_back(index_i);
                    } else if (inter[1] > 2 * M_PI) {
                        intervals_angle.push_back(inter[0]);
                        intervals_angle.push_back(2 * M_PI);
                        intervals_index.push_back(index_i);
                        intervals_index.push_back(index_i);

                        intervals_angle.push_back(0);
                        intervals_angle.push_back(inter[1] - 2 * M_PI);
                        intervals_index.push_back(index_i);
                        intervals_index.push_back(index_i);
                    } else {
                        intervals_angle.push_back(inter[0]);
                        intervals_angle.push_back(inter[1]);
                        intervals_index.push_back(index_i);
                        intervals_index.push_back(index_i);
                    }
                }
            }
        }
    }

    double stabber;
    std::vector<size_t> index_stabbed;
    interval_stabbing(intervals_angle, intervals_index, stabber, index_stabbed);

    theta_est = stabber;

    inlier_index_after_theta_est = index_stabbed;

}



template <typename T>
std::vector<size_t> registration::sort_indexes(std::vector<T> &v)
{
    std::vector<size_t> idx(v.size());

    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

    return idx;
}


void registration::interval_stabbing(std::vector<double> &intervals_phi, std::vector<size_t> &intervals_index, double &stabber, std::vector<size_t> &index_stabbed){

    size_t num_endpoints = intervals_phi.size();
    size_t num_intervals = num_endpoints / 2;

    int mask[num_endpoints];
    for(int i =0;i<num_intervals;i++){
        mask[2*i] = 0;
        mask[2*i+1] = 1;
    }

    std::vector<size_t> indexs = sort_indexes(intervals_phi);
    std::sort(intervals_phi.begin(), intervals_phi.end());

    size_t max_count = 0;
    stabber = 0;
    int num_stabbed= 0;
    std::vector<size_t> stabbed;
    for(int j =0;j<num_endpoints;j++){
        if(mask[indexs[j]] == 0) {
            max_count++;
            stabbed.push_back(intervals_index[indexs[j]]);
            if(max_count>num_stabbed){
                num_stabbed = max_count;
//                stabber = (intervals_phi[j] + intervals_phi[j+1])/2;
                stabber = intervals_phi[j];
                index_stabbed = stabbed;
            }
        }
        else{
            max_count = max_count - 1;
            stabbed.erase(std::remove(stabbed.begin(), stabbed.end(), intervals_index[indexs[j]]), stabbed.end());
        }
    }
}


void registration::post_SVD(cloudPoints &src, cloudPoints &tgt, double xi, Eigen::Matrix3d &R_method, Eigen::Vector3d &t_method, Eigen::Matrix3d &R_post, Eigen::Vector3d &t_post){


    int N = src.cols();
    int inliers_num = 0;
    std::vector<int> inliers_index;

    for(int i =0 ;i<N;i++){

        Eigen::Vector3d src_i_vector = src.col(i);
        Eigen::Vector3d tgt_i_vector = tgt.col(i);

        double error = (tgt_i_vector - R_method * src_i_vector - t_method).norm();

        if(error <= xi){
            inliers_index.push_back(i);
            inliers_num++;
        }

    }

    cloudPoints inliers_src(3, inliers_num);
    cloudPoints inliers_tgt(3, inliers_num);

    for(int j = 0; j <inliers_num ;j++){
        inliers_src.col(j) = src.col(inliers_index[j]);
        inliers_tgt.col(j) = tgt.col(inliers_index[j]);
    }

    Eigen::Matrix4d transformation_matrix;
    transformation_matrix = pcl::umeyama(inliers_src, inliers_tgt);
    R_post = transformation_matrix.block<3,3>(0, 0);
    t_post = transformation_matrix.block<3,1>(0, 3);

}


void registration::post_refinement(cloudPoints &src, cloudPoints &tgt, double xi, Eigen::Matrix3d &R_method, Eigen::Vector3d &t_method, Eigen::Matrix3d &R_post, Eigen::Vector3d &t_post, int num){

    Eigen::Matrix3d R_flag;
    Eigen::Vector3d t_flag;

    R_flag = R_method;
    t_flag = t_method;
    for(int i = 0; i< num; i++){

        post_SVD(src, tgt, xi, R_flag, t_flag, R_post, t_post);

        R_flag = R_post;
        t_flag = t_post;
    }

}