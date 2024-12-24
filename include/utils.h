#ifndef DEREG_UTILS_H
#define DEREG_UTILS_H

#include <vector>
#include <math.h>
#include <assert.h>

#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/Core>

#include <pcl/io/io.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_cloud.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/common/centroid.h>
#include <pcl/features/fpfh_omp.h>
#include <pcl/filters/random_sample.h>


#define BRANCH_ACCURACY 1e-3

typedef Eigen::Matrix<double, 3, Eigen::Dynamic> cloudPoints;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> cloudNorms;


namespace HERE{

    namespace utils{

        typedef struct state1Dof{

            double l_,r_;

            state1Dof(): l_(-1), r_(1){
            }

            state1Dof(double l, double r): l_(l), r_(r){

            }

            double mean(){
                return 0.5*(l_+r_);
            }

            double edge(){
                return 0.5 * (r_ - l_ );
            }

        }state1Dof;

        inline
        std::ostream &operator<<(std::ostream& os, const state1Dof &reg)
        {
            os << "[("<< reg.l_ <<" "<< reg.r_ <<")]";
            return os;
        }

        inline
        size_t split(state1Dof father, state1Dof **children, double accuracy = BRANCH_ACCURACY){
            if (father.edge() <= accuracy)
                return 0;
            double c = father.mean();
            children[0] = new state1Dof(father.l_,c);
            children[1] = new state1Dof(c, father.r_);
            return 2;
        }


        template<class SSR, typename Scalar=int>
        class SearchState
        {
        public:
            SSR ssr;
            Scalar bnd;
            SearchState(SSR ssr, Scalar bnd): ssr(ssr), bnd(bnd) {}
            ~SearchState()= default;

            friend std::ostream &operator<<(std::ostream& os, const SearchState &ss)
            {
                os<< "ssr "<<ss.ssr<<" "<<" bnd "<<ss.bnd ;
                return os;
            }
        };
    }

}

#endif